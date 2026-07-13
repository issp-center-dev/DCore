import subprocess
import itertools
import numpy as np
import os
import sys

# T_list, n_iw, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output"

def calc_one_body_green_core_parallel(p_common, max_workers=None, sector_occupancy=None, ne_only=None):
    """
    Return:
        np.ndarray(n_site, n_sigma, n_site, n_sigma, n_T, n_omega)

    sector_occupancy: optional (n_up, n_down, n_site) for the canonical (sector-resolved)
        path. When given, operators with an exactly-zero contribution in that (Ne, 2Sz)
        sector -- cross-spin off-diagonal terms (zero by spin conservation), and
        annihilation/creation on an already empty/full spin -- are skipped instead of run,
        which also avoids HPhi excitations into non-existent sectors. Their contribution is
        filled with zeros. None (default) is the grand-canonical path: every operator runs.

    ne_only: optional total electron number Ne of a 2Sz-free (HubbardNConserved) sector, for the
        spin-orbit bra/ket route (DCORE_HPHI_BRAKET=1 with spin-orbit). Both spins share the
        Ne+-1 excited space, so cross-spin G elements are computed too. Mutually exclusive with
        sector_occupancy; only meaningful on the braket path.
    """

    n_sigma = 2
    n_flg = 2
    n_excitation = 2
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut, *rest = p_common

    if n_site < 1:
        raise ValueError("No excitation tasks were generated (n_site must be >= 1).")

    # Stage-3 direct bra/ket path (opt-in). One ket BiCG solve is projected onto all bras via
    # HPhi's SpectrumNumBra, so off-diagonal G is read directly instead of reconstructed from the
    # c_i + i c_j combination trick (BiCG count n_orb^2 -> n_orb). Default off = combination trick.
    # These guards run before check_eta so the unsupported configurations fail cleanly.
    if os.environ.get("DCORE_HPHI_BRAKET", "0") == "1":
        # The bra/ket path runs on a sector decomposition: per (Ne, 2Sz) for the 2Sz-conserving
        # case (same-spin only), or per Ne (ne_only) for the spin-orbit / HubbardNConserved case
        # (both spins, cross-spin G included). The grand-canonical route (all spins in one excited
        # Fock space without sectoring) is not supported -- reject it rather than silently
        # producing a wrong (zero-norm) Green's function.
        if ne_only is None and sector_occupancy is None:
            raise RuntimeError("DCORE_HPHI_BRAKET=1 requires a sector (sector_occupancy or "
                               "ne_only); the grand-canonical bra/ket route is unsupported.")
        check_eta(p_common)
        return _braket_one_body_green(p_common, max_workers, sector_occupancy, ne_only=ne_only)

    check_eta(p_common)

    def composite(site, sigma):
        return site * n_sigma + sigma

    def gen_p():
        for sitei, sigmai in itertools.product(range(n_site), range(n_sigma)):
            for sitej, sigmaj in itertools.product(range(n_site), range(n_sigma)):
                # Only the upper triangle (a_i <= a_j) is computed; the transposed
                # element G_ji is reconstructed from the same excitation data.
                if composite(sitei, sigmai) > composite(sitej, sigmaj):
                    continue
                for idx, flg in enumerate([True, False]):
                    for ex_state in range(n_excitation):
                        # On the diagonal only the flg=True (idx=0) excitation is used.
                        if composite(sitei, sigmai) == composite(sitej, sigmaj) and idx == 1:
                            continue
                        yield sitei, sigmai, sitej, sigmaj, idx, ex_state, p_common

    tasks = list(gen_p())
    if not tasks:
        raise ValueError("No excitation tasks were generated (n_site must be >= 1).")
    from concurrent.futures import ProcessPoolExecutor
    import shutil

    # In the canonical path, only operators with a non-zero contribution in this sector are
    # run; the rest are filled with zeros below. In the grand-canonical path all tasks run.
    if sector_occupancy is not None:
        run_tasks = [t for t in tasks if _task_valid_in_sector(t, *sector_occupancy)]
        if not run_tasks:
            raise RuntimeError("No valid excitation operator for sector {}.".format(sector_occupancy))
    else:
        run_tasks = tasks

    internal_loop = os.environ.get("DCORE_HPHI_INTERNAL_LOOP", "0") == "1"
    if internal_loop:
        # Group operators by Hilbert sector so each HPhi run reads every eigenvector
        # once and reuses it across all operators of that sector (op-inner). Each
        # batch is one HPhi launch; the pool parallelizes over batches.
        batches = _group_tasks_by_sector(run_tasks)
        cleanup_dirs = ["batch_{}".format(bid) for bid in range(len(batches))]
    else:
        cleanup_dirs = ["{}_{}_{}_{}_{}_{}".format(t[0], t[1], t[2], t[3], t[5], t[4])
                        for t in run_tasks]

    try:
        if internal_loop:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                batch_outputs = list(executor.map(calc_one_body_green_batch, list(enumerate(batches))))
            result_map = {}
            for out in batch_outputs:
                result_map.update(out)
        else:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                run_results = list(executor.map(calc_one_body_green_core, run_tasks))
            result_map = {_task_key(t): r for t, r in zip(run_tasks, run_results)}

        assert len(result_map) == len(run_tasks), \
            "routing covered {} of {} run tasks".format(len(result_map), len(run_tasks))
        n_omega = next(iter(result_map.values())).shape[-1]
        # tasks not run (filtered out as exactly zero in this sector) contribute zeros.
        zero = np.zeros((len(T_list), n_omega), dtype=np.complex128)
        results = [result_map.get(_task_key(t), zero) for t in tasks]

        one_body_green_core = np.zeros(
            (n_site, n_sigma, n_site, n_sigma, n_flg, n_excitation, len(T_list), n_omega),
            dtype=np.complex128)
        for (sitei, sigmai, sitej, sigmaj, idx, ex_state, _), res in zip(tasks, results):
            one_body_green_core[sitei][sigmai][sitej][sigmaj][idx][ex_state] = res
        one_body_green = calc_one_body_green(one_body_green_core)
    finally:
        for dir_path in cleanup_dirs:
            if os.path.isdir(dir_path):
                shutil.rmtree(dir_path)

    return one_body_green


def _task_valid_in_sector(task, n_up, n_down, n_site):
    """Whether a one-body-Gf operator has a non-zero contribution in a canonical (Ne, 2Sz) sector.

    Cross-spin off-diagonal operators (sigma_i != sigma_j) vanish by spin conservation and would
    map to two different Sz sectors, so the c_i + i c_j combination is ill-defined canonically.
    A same-spin operator excites spin sigma: annihilation (ex_state 0) needs that spin occupied
    (n_sigma >= 1); creation (ex_state 1) needs a free slot (n_sigma <= n_site - 1). Otherwise the
    excited state is zero (and HPhi would excite into a non-existent sector).
    """
    sitei, sigmai, sitej, sigmaj, idx_flg, ex_state = task[:6]
    if sigmai != sigmaj:
        return False
    n_sigma = n_up if sigmai == 0 else n_down
    if ex_state == 0:  # annihilation c_sigma
        return n_sigma >= 1
    return n_sigma <= n_site - 1  # creation c_sigma^dagger


def _task_key(task):
    """Hashable identity of a task (drops the unhashable p_common)."""
    return task[:6]


def _group_tasks_by_sector(tasks):
    """Group operators that can share one HPhi run (same Hilbert sector).

    A single-excitation operator changes the electron number by +-1 in spin
    sigma. Diagonal (c_{i,sigma}) and same-spin off-diagonal (c_{i,sigma} +
    i c_{j,sigma}) operators map to the sector keyed by (ex_state, sigma) and are
    batched together. Cross-spin off-diagonal operators (sigma_i != sigma_j) mix
    Sz, so they are sector-inconsistent and run as singletons (one launch each).
    """
    groups = {}
    singletons = []
    for t in tasks:
        sitei, sigmai, sitej, sigmaj, i_flg, ex_state, _ = t
        if sigmai == sigmaj:
            groups.setdefault((ex_state, sigmai), []).append(t)
        else:
            singletons.append([t])
    return list(groups.values()) + singletons

def calc_one_body_green_core(p):
    #unpack parameters
    sitei, sigmai, sitej, sigmaj, i_flg, ex_state, p_common = p
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut, *rest = p_common
    mpi_prefix = rest[0] if rest else ""
    calc_spectrum_core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi=path_to_HPhi, header=header,
                                           output_dir=output_dir, mpi_prefix=mpi_prefix)

    calc_spectrum_core.set_energies()
    flg = True if i_flg == 0 else False
    return calc_spectrum_core.get_one_body_green_core(sitei, sigmai, sitej, sigmaj, ex_state, flg, exct_cut)

def calc_one_body_green_batch(payload):
    """Run one HPhi launch for a batch of same-sector operators (op-inner reuse).

    payload = (batch_id, batch); batch is a list of tasks sharing a Hilbert sector.
    Returns {task_key: one_body_green} for every task in the batch.
    """
    batch_id, batch = payload
    p_common = batch[0][6]
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut, *rest = p_common
    mpi_prefix = rest[0] if rest else ""
    core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi=path_to_HPhi, header=header,
                            output_dir=output_dir, mpi_prefix=mpi_prefix)
    core.set_energies()
    # (sitei, sigmai, sitej, sigmaj, ex_state, flg); flg = (i_flg == 0)
    ops = [(t[0], t[1], t[2], t[3], t[5], t[4] == 0) for t in batch]
    res_list = core.get_one_body_green_batch(ops, exct_cut, batch_id=batch_id)
    assert len(res_list) == len(batch), \
        "batch {}: got {} results for {} operators".format(batch_id, len(res_list), len(batch))
    return {_task_key(t): r for t, r in zip(batch, res_list)}

def calc_one_body_green_braket_job(payload):
    """Run one HPhi braket launch (all kets x all bras) for one ex_state sector.

    payload = (job_id, ex_state, ket_ops, bra_ops, p_common); ops are (site, sigma).
    Returns {(bra_op, ket_op): one_body_green} (the sign for this ex_state is already applied).
    """
    job_id, ex_state, ket_ops, bra_ops, p_common = payload
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut, *rest = p_common
    mpi_prefix = rest[0] if rest else ""
    core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi=path_to_HPhi, header=header,
                            output_dir=output_dir, mpi_prefix=mpi_prefix)
    core.set_energies()
    return core.get_one_body_green_braket_batch(ex_state, ket_ops, bra_ops, exct_cut, batch_id=job_id)


def _braket_sector_jobs(n_site, sector_occupancy):
    """Enumerate the (ex_state, ket_ops, bra_ops) launches of the direct bra/ket path.

    One launch per (ex_state, spin): ket_ops = bra_ops = [(0, sigma), .., (n_site-1, sigma)] so the
    full same-spin G block is obtained from n_site solves. sector_occupancy = (n_up, n_down, n_site)
    drops (ex_state, spin) sectors whose excited state vanishes (annihilation on an empty spin or
    creation on a full one); None (grand canonical) keeps every (ex_state, spin). Cross-spin
    elements are never generated (zero by spin conservation; dropped from the final Gimp anyway).
    """
    n_sigma = 2
    n_excitation = 2
    jobs = []
    for ex_state in range(n_excitation):
        for sigma in range(n_sigma):
            if sector_occupancy is not None:
                n_up, n_down, nsite = sector_occupancy
                n_s = n_up if sigma == 0 else n_down
                # annihilation needs the spin occupied; creation needs a free slot.
                valid = (n_s >= 1) if ex_state == 0 else (n_s <= nsite - 1)
                if not valid:
                    continue
            ops = [(site, sigma) for site in range(n_site)]
            jobs.append((ex_state, ops, ops))
    return jobs


def _braket_ne_sector_jobs(n_site, ne):
    """Bra/ket launches for ONE Ne-only sector (2Sz-free, HubbardNConserved) -- the spin-orbit
    capable route. With only Ne fixed, both spins share the Ne+-1 excited space, so each launch's
    kets = bras = ALL 2*n_site spin-orbitals (both spins); one ket solve then yields the same-spin
    AND cross-spin G_{i sigma_i, j sigma_j} elements. One launch per ex_state: annihilation (0)
    needs an electron to remove (Ne >= 1); creation (1) needs a free slot (Ne <= 2*n_site - 1).
    """
    n_sigma = 2
    n_so = n_sigma * n_site
    if not 0 <= ne <= n_so:
        raise ValueError("ne must be in [0, 2*n_site] = [0, {}], got {}".format(n_so, ne))
    jobs = []
    for ex_state in range(2):
        valid = (ne >= 1) if ex_state == 0 else (ne <= n_so - 1)
        if not valid:
            continue
        ops = [(site, sigma) for site in range(n_site) for sigma in range(n_sigma)]
        jobs.append((ex_state, ops, ops))
    return jobs


def _braket_store_op(bra_op, ket_op, ex_state):
    """Where the bra/ket element g = G_{bra,ket} = <c_bra phi|R|c_ket phi> is stored in
    one_body_green, as a (row_op, col_op) pair of (site, sigma) tuples.

    The annihilation channel (ex_state 0) stores the physical element at one_body_green[ket][bra]
    -- the convention the validated combination-trick path uses. The creation channel (ex_state 1)
    builds bra and ket from c^dag operators, whose ordering reverses the two indices, so it is the
    transpose, stored at [bra][ket]. (For a diagonal element bra == ket the two coincide.) This is
    the SAME convention for same-spin and cross-spin elements; it was verified element-wise (to
    1e-10) against the combination trick / scipy for the same-spin blocks, and the cross-spin
    blocks (Ne-only route) rely on the identical operator-ordering argument.
    """
    return (ket_op, bra_op) if ex_state == 0 else (bra_op, ket_op)


def _braket_assemble(results, n_site):
    """Accumulate the per-launch braket results into the one_body_green tensor.

    ``results`` is a list of dicts ``{(row_op, col_op): g}`` (one per HPhi launch), where row_op /
    col_op are (site, sigma) storage indices already produced by _braket_store_op (so the
    creation-channel transpose is baked in). Returns one_body_green of shape
    (n_site, 2, n_site, 2, n_T, n_omega); a given (row_op, col_op) accumulates across ex_state
    launches (annihilation + creation) and, in the Ne-only route, across the both-spin job.
    """
    n_sigma = 2
    sample = next(iter(results[0].values()))
    n_T, n_omega = sample.shape
    one_body_green = np.zeros((n_site, n_sigma, n_site, n_sigma, n_T, n_omega), dtype=np.complex128)
    for res in results:
        for (row_op, col_op), g in res.items():
            (ri, rsg) = row_op
            (ci, csg) = col_op
            one_body_green[ri][rsg][ci][csg] += g
    return one_body_green


def _braket_one_body_green(p_common, max_workers, sector_occupancy, ne_only=None):
    """Stage-3 direct bra/ket one-body Green's function (DCORE_HPHI_BRAKET=1).

    Replaces the c_i + i c_j combination trick: each ket solve is projected onto every bra,
    so G_{i,j} is read directly (BiCG count n_orb^2 -> n_orb).

    Two sector modes:
    - (Ne, 2Sz) [default, ``ne_only=None``]: each launch handles one (ex_state, spin) sector,
      kets = bras = the same-spin sites, so the same-spin G block comes from n_site solves.
      Cross-spin elements vanish by spin conservation and are not computed.
    - Ne-only [``ne_only`` = the sector's Ne; for the spin-orbit / HubbardNConserved route]:
      both spins share the Ne+-1 excited space, so each launch's kets = bras = ALL 2*n_site
      spin-orbitals and one ket solve yields the same-spin AND cross-spin G elements.
    """
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut, *rest = p_common

    if ne_only is not None:
        # Ne-only mode is mutually exclusive with the (Ne, 2Sz) occupancy: the two describe
        # different sector schemes, so a caller passing both is a bug, not a silent fallback.
        if sector_occupancy is not None:
            raise ValueError("ne_only and sector_occupancy are mutually exclusive "
                             "(Ne-only vs (Ne, 2Sz) sectoring); got both.")
        spec = _braket_ne_sector_jobs(n_site, ne_only)
    else:
        spec = _braket_sector_jobs(n_site, sector_occupancy)
    if not spec:
        raise RuntimeError("No valid braket excitation sector for occupancy={}, ne_only={}."
                           .format(sector_occupancy, ne_only))
    jobs = [(ex_state, ket_ops, bra_ops, p_common) for (ex_state, ket_ops, bra_ops) in spec]

    from concurrent.futures import ProcessPoolExecutor
    import shutil
    payloads = [(jid,) + job for jid, job in enumerate(jobs)]
    cleanup_dirs = ["braket_{}".format(jid) for jid in range(len(jobs))]
    try:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(calc_one_body_green_braket_job, payloads))
        one_body_green = _braket_assemble(results, n_site)
    finally:
        for d in cleanup_dirs:
            if os.path.isdir(d):
                shutil.rmtree(d)
    return one_body_green


def check_eta(p_common):
    _, T_list, exct, eta, path_to_HPhi, header, output_dir, _, *rest = p_common
    mpi_prefix = rest[0] if rest else ""
    calc_spectrum_core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi=path_to_HPhi, header=header,
                                           output_dir=output_dir, mpi_prefix=mpi_prefix)
    calc_spectrum_core.set_energies(check_eta=True)

def calc_one_body_green(one_body_green_core):
    n_site, n_sigma, n_site, n_sigma, n_excitation, n_flg, n_T, n_omega = one_body_green_core.shape
    n_excitation = 2
    n_flg = 2
    n_sigma = 2
    one_body_green = np.zeros((n_site, n_sigma, n_site, n_sigma, n_T, n_omega), dtype=np.complex128)
    # Diagonal
    for sitei, sigmai, ex_state in itertools.product(range(n_site), range(n_sigma), range(n_excitation)):
        one_body_green[sitei][sigmai][sitei][sigmai] += one_body_green_core[sitei][sigmai][sitei][sigmai][0][ex_state]

    # Off diagonal: only the upper triangle (a_i < a_j) is reconstructed; each
    # iteration fills both G_ij and the transposed G_ji from the same data.
    for sitei, sigmai, sitej, sigmaj  in itertools.product(range(n_site), range(n_sigma), range(n_site), range(n_sigma)):
        if sitei * n_sigma + sigmai >= sitej * n_sigma + sigmaj:
            continue  # diagonal is set above; lower triangle is set via its transpose
        one_body_green_tmp = np.zeros((n_flg, n_T, n_omega), dtype=np.complex128)
        for idx in range(n_flg):
            for ex_state in range(n_excitation):
                one_body_green_tmp[idx] += one_body_green_core[sitei][sigmai][sitej][sigmaj][idx][ex_state]
        # Subtract the diagonal contribution from both combinations
        #   B = c_i + c_j   -> tmp[1] - (G_ii + G_jj) = G_ij + G_ji
        #   A = c_i + i c_j -> tmp[0] - (G_ii + G_jj) = i (G_ji - G_ij)
        diag = one_body_green[sitei][sigmai][sitei][sigmai] + \
               one_body_green[sitej][sigmaj][sitej][sigmaj]
        one_body_green_tmp[0] -= diag
        one_body_green_tmp[1] -= diag
        one_body_green[sitei][sigmai][sitej][sigmaj] = (one_body_green_tmp[1] + 1J * one_body_green_tmp[0]) / 2.0
        one_body_green[sitej][sigmaj][sitei][sigmai] = (one_body_green_tmp[1] - 1J * one_body_green_tmp[0]) / 2.0
    return one_body_green


class CalcSpectrumCore:
    def __init__(self, T_list, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output", mpi_prefix=""):
        self.T_list = T_list
        self.exct = exct
        self.eta = eta
        self.header = header
        self.output_dir = output_dir
        self.nomega = 0
        self.parent_dir = os.getcwd()
        # self.path_to_HPhi = os.path.join(self.parent_dir, path_to_HPhi)
        self.path_to_HPhi = os.path.abspath(path_to_HPhi)  # converted to full path in DCore
        # MPI launcher for each HPhi run of the Gf step (e.g. "mpirun -np 4");
        # empty string runs HPhi serially (one rank, no mpirun).
        self.mpi_prefix = mpi_prefix
        # When enabled, HPhi loops over the exct_cut eigenstates internally
        # (SpectrumLoopExct in modpara): one launch per excitation operator instead
        # of one per (eigenstate x operator), eliminating the per-state process-startup
        # overhead. Requires an HPhi build with the SpectrumLoopExct feature. The
        # per-idx path (default) stays as the validated fallback / reference.
        self.internal_loop = os.environ.get("DCORE_HPHI_INTERNAL_LOOP", "0") == "1"

    def Make_Spectrum_Input(self, calc_dir="./", spectrum_type="single", bra=False):

        rel_path_org = os.path.relpath(self.parent_dir, calc_dir)
        rel_path = os.path.relpath(os.path.join(self.parent_dir, "output"), calc_dir)
        spvec_base = os.path.join("../", rel_path, self.header) + "_eigenvec"

        def _write_calcmod():
            with open(os.path.join(self.parent_dir, "calcmod.def")) as f:
                lines = f.readlines()
            with open(os.path.join(calc_dir, "calcmod_ex.def"), "w") as fex:
                for line in lines:
                    words = line.split()
                    if words[0] in ("CalcSpec", "OutputExVec", "OutputEigenVec"):
                        continue
                    fex.write(line)
                fex.write("CalcSpec    1\n")

        def _write_namelist(fname, spectrum_vec):
            with open(os.path.join(self.parent_dir, "namelist.def")) as f:
                lines = f.readlines()
            with open(os.path.join(calc_dir, fname), "w") as fex:
                for line in lines:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    if words[0] in ["CalcMod", "SpectrumVec", "ModPara",
                                    "SingleExcitation", "PairExcitation",
                                    "SingleExcitationBra", "PairExcitationBra"]:
                        continue
                    fex.write("{} {}\n".format(words[0], os.path.join(rel_path_org, words[1])))
                fex.write("ModPara modpara_ex.def\n")
                fex.write("CalcMod calcmod_ex.def\n")
                fex.write("SpectrumVec    {}\n".format(spectrum_vec))
                if spectrum_type == "single":
                    fex.write("SingleExcitation single_ex.def\n")
                    if bra:
                        # Stage-3 bra/ket reuse: bra set 0 is the namelist SingleExcitationBra;
                        # bra sets 1.. come from single_ex_bra_<b>.def (SpectrumNumBra).
                        fex.write("SingleExcitationBra single_ex_bra_0.def\n")
                elif spectrum_type == "pair":
                    fex.write("PairExcitation pair_ex.def\n")

        _write_calcmod()
        if self.internal_loop:
            # One launch: HPhi loops eigenstates internally and appends
            # _<idx>_rank_<r>.dat to this common SpectrumVec base per eigenstate.
            _write_namelist("namelist_ex.def", spvec_base)
            # HPhi (loop mode) reads E_idx from <calc_dir>/output/<header>_energy.dat
            # (childfopenMPI prepends 'output/'); stage it from the eigenvalue run.
            import shutil
            os.makedirs(os.path.join(calc_dir, "output"), exist_ok=True)
            shutil.copy(os.path.join(self.parent_dir, self.output_dir, "{}_energy.dat".format(self.header)),
                        os.path.join(calc_dir, "output", "{}_energy.dat".format(self.header)))
        else:
            for idx in range(self.exct):
                _write_namelist("namelist_ex_{}.def".format(idx),
                                "{}_{}".format(spvec_base, idx))

        with open(os.path.join(self.parent_dir,"modpara.def"), "r") as fr:
            lines = fr.readlines()
            for line in lines:
                words = line.split()
                if words[0] == "NOmega":
                    self.nomega = int(words[1])

        if self.nomega == 0:
            print("Error: Please set NOmega in modpara file")
            sys.exit(1)

    def _read_spectrum(self, calc_dir="./"):
        spectrum_dict={}
        frequencies =[]
        for idx in range(self.exct):
            path_to_spectrum_dir = os.path.join(calc_dir, self.output_dir)
            path_to_DG = os.path.join(path_to_spectrum_dir, "{}_DynamicalGreen_{}.dat".format(self.header,idx))
            spectrum = np.loadtxt(path_to_DG)
            spectrum_dict[idx] = spectrum[:,2] + 1J*spectrum[:,3]
            if idx == 0 :
                frequencies = spectrum[:, 0] + 1J*spectrum[:, 1]
        frequencies = frequencies
        return frequencies, spectrum_dict

    def set_energies(self, check_eta=False):
        energy_list = []
        with open(os.path.join(self.output_dir, "{}_energy.dat".format(self.header))) as f:
            lines = f.readlines()
            for line in lines:
                words = line.split()
                if len(words) != 0 and words[0] == "Energy":
                    energy_list.append(float(words[1]))
        self.energy_list = energy_list
        # Use min/max (not [0]/[-1]) so the sector-local Boltzmann normalization here matches the
        # cross-sector recombination in the canonical solver, which uses energies.min(), even if
        # zvo_energy.dat were ever written out of energy order.
        self.ene_min = min(energy_list)
        self.ene_max = max(energy_list)

        if check_eta:
            print(f"\n  Check eta:=exp[-beta(ene_max-ene_mix)] < {self.eta:.1e}")
            for T in self.T_list:
                eta_ene = np.exp(-(self.ene_max-self.ene_min)/T)
                print(f"    T = {T}: eta = {eta_ene:.2e}")
                if eta_ene > self.eta:
                    print(f"Warning: At T = {T}, exp[-beta(ene_max-ene_mix)]={eta_ene:.2e} is larger than eta={self.eta}.", file=sys.stderr)

    def _calc_Z(self, T):
        Z = 0
        for ene in self.energy_list:
            ene_diff = ene-self.ene_min
            Z += np.exp(-ene_diff/T)
        return Z

    def get_finite_T_spectrum(self, calc_dir="./"):
        frequencies, self.spectrums_dict = self._read_spectrum(calc_dir)
        finite_T_spectrum_dict ={}
        for T in self.T_list:
            Z = self._calc_Z(T)
            spectrum = np.zeros_like(self.spectrums_dict[0])
            for idx in range(self.exct):
                spectrum += np.exp(-(self.energy_list[idx]-self.ene_min)/T)*self.spectrums_dict[idx]
            spectrum /= Z
            finite_T_spectrum_dict[T]=spectrum
        self.finite_T_spectrum_dict = finite_T_spectrum_dict
        return frequencies, finite_T_spectrum_dict

    def print_finite_T_spectrum(self, file_name = "Dynamical_Green"):
        for key, spectrum in self.finite_T_spectrum_dict.items():
            file_name_T = self.header + "_" + file_name + "_T_{}.dat".format(key)
            with open(os.path.join(self.output_dir, file_name_T), "w") as fw:
                for idx, value in enumerate(spectrum):
                    fw.write("{} {} {} {}\n".format(self.frequencies[idx].real, self.frequencies[idx].imag, value.real, value.imag))

    def _update_modpara(self, exct, ex_state=0, calc_dir="./", num_op=1, num_bra=1):
        dict_mod={}
        header = []
        with open(os.path.join(self.parent_dir, "modpara.def"), "r") as fr:
            lines = fr.readlines()
            header = lines[:8]
            for line in lines[8:]:
                words = line.split()
                dict_mod[words[0]] = words[1:]
            if self.internal_loop:
                # 'exct' is the loop COUNT (= exct_cut); HPhi reads each E_idx from
                # the energy file and uses it as the per-state shift internally.
                dict_mod["SpectrumLoopExct"] = [exct]
                if num_op > 1:
                    # HPhi evaluates num_op operator sets (op-inner) reading each
                    # eigenvector once; set 0 is single_ex.def, sets 1.. are single_ex_<op>.def.
                    dict_mod["SpectrumNumOp"] = [num_op]
                if num_bra > 1:
                    # Stage-3 bra/ket reuse: project each ket solve onto num_bra bras in one
                    # BiCG run; set 0 is single_ex_bra_0.def, sets 1.. are single_ex_bra_<b>.def.
                    dict_mod["SpectrumNumBra"] = [num_bra]
            else:
                dict_mod["OmegaOrg"] = [self.energy_list[exct], 0]
            if ex_state == 0:
                omega_max = dict_mod["OmegaMax"]
                dict_mod["OmegaMax"] = [-1.0*float(omega_max[0]), -1.0*float(omega_max[1])]
                omega_min = dict_mod["OmegaMin"]
                dict_mod["OmegaMin"] = [-1.0*float(omega_min[0]), -1.0*float(omega_min[1])]
            with open(os.path.join(calc_dir, "modpara_ex.def"), "w") as fw:
                for line in header:
                    fw.write(line)
                for key, value in dict_mod.items():
                    if len(value) == 1:
                        fw.write("{} {}\n".format(key, value[0]))
                    else:
                        fw.write("{} {} {}\n".format(key, value[0], value[1]))

    def _make_single_excitation(self, site_i, sigma_i, site_j, sigma_j, file_name = "single_ex.def", ex_state=0, flg_complex = True, calc_dir="./"):
        # c_{i sigma_i} or c_{i sigma_i} + i c_{j sigma_j}
        nsingle = 2
        if (2 * site_i + sigma_i) == ( 2 * site_j + sigma_j):
            nsingle = 1
        with open(os.path.join(calc_dir, file_name), "w") as fw:
            fw.write("===============================\n")
            fw.write("NSingle {}\n".format(nsingle))
            fw.write("===============================\n")
            fw.write("===============================\n")
            fw.write("===============================\n")
            if nsingle == 1:
                fw.write("{} {} {} 1.0 0.0\n".format(site_i, sigma_i, ex_state))
            else:
                if flg_complex is True:
                    # A = c_i + i c_j (annihilation, ex_state=0), whose conjugate
                    # excitation A^dag = c_i^dag - i c_j^dag (creation, ex_state=1).
                    # The imaginary coefficient of c_j must be conjugated for the
                    # creation channel, otherwise the anti-symmetric (imaginary)
                    # part of the off-diagonal Green's function gets the wrong
                    # high-frequency tail. (ex_state: 0 = c, 1 = c^dag)
                    imag_coeff = -1.0 if ex_state == 1 else 1.0
                    fw.write("{} {} {} 1.0 0.0\n".format(site_i, sigma_i, ex_state))
                    fw.write("{} {} {} 0.0 {}\n".format(site_j, sigma_j, ex_state, imag_coeff))
                else:
                    fw.write("{} {} {} 1.0 0.0\n".format(site_i, sigma_i, ex_state))
                    fw.write("{} {} {} 1.0 0.0\n".format(site_j, sigma_j, ex_state))

    def _run_HPhi(self, exct_cut, ex_state=0, calc_dir="./", num_op=1, num_bra=1):
        os.chdir(calc_dir)
        exec_path = self.path_to_HPhi
        if self.internal_loop:
            # One launch covers all exct_cut eigenstates (and all num_op x num_bra operator
            # sets); HPhi writes output/<header>_DynamicalGreen_<idx>[_<op>[_<bra>]].dat directly.
            self._update_modpara(exct_cut, ex_state, calc_dir, num_op=num_op, num_bra=num_bra)
            input_path = os.path.join(calc_dir, "namelist_ex.def")
            cmd = "{} {} -e {} > std.log".format(self.mpi_prefix, exec_path, input_path).strip()
            ret = subprocess.call(cmd, shell=True)
            if ret != 0:
                raise RuntimeError(
                    "HPhi spectrum run failed (exit {}) in {}; see std.log".format(ret, calc_dir))
        else:
            for idx in range(exct_cut):
                self._update_modpara(idx, ex_state, calc_dir)
                input_path = os.path.join(calc_dir, "namelist_ex_{}.def".format(idx))
                cmd = "{} {} -e {} > std_{}.log".format(self.mpi_prefix, exec_path, input_path, idx).strip()
                subprocess.call(cmd, shell=True)
                cmd = "mv ./output/{0}_DynamicalGreen.dat ./output/{0}_DynamicalGreen_{1}.dat".format(self.header, idx)
                subprocess.call(cmd, shell=True)
        os.chdir(self.parent_dir)

    def get_one_body_green_core(self, sitei, sigmai, sitej, sigmaj, ex_state, flg, exct_cut):
        calc_dir = os.path.join(self.parent_dir, "{}_{}_{}_{}_{}_{}".format(sitei,sigmai,sitej,sigmaj, ex_state, 0 if flg is True else 1))
        os.makedirs(calc_dir, exist_ok=True)
        self.Make_Spectrum_Input(calc_dir)
        one_body_green = np.zeros((len(self.T_list), self.nomega), dtype=np.complex128)
        # print("Calculate G[{},{}][{},{}]".format(sitei, "u" if sigmai == 0 else "d", sitej, "u" if sigmaj == 0 else "d"))
        self._make_single_excitation(sitei, sigmai, sitej, sigmaj, ex_state=ex_state, flg_complex=flg, calc_dir=calc_dir)
        # Run HPhi
        self._run_HPhi(exct_cut, ex_state, calc_dir)
        # Get Finite-T Green
        frequencies, finite_spectrum_list = self.get_finite_T_spectrum(calc_dir)
        if ex_state == 1:
            self.frequencies = frequencies
        sign = 1.0 if ex_state == 1 else -1.0
        for idx, T in enumerate(self.T_list):
            one_body_green[idx] = sign * finite_spectrum_list[T]
        return one_body_green

    def _finite_T_spectrum_multiop(self, calc_dir, num_op):
        """Boltzmann-sum the per-(idx, op) spectra written by a SpectrumNumOp run.

        Returns (frequencies, {op: {T: spectrum}}). Reads
        <header>_DynamicalGreen_<idx>_<op>.dat for op in 0..num_op-1 and the same
        self.exct eigenstates as get_finite_T_spectrum.
        """
        spec_dir = os.path.join(calc_dir, self.output_dir)
        frequencies = None
        raw = {}  # (op, idx) -> complex spectrum
        for op in range(num_op):
            for idx in range(self.exct):
                # HPhi appends the _<op> suffix only when there is more than one
                # operator set (useOp = nop>1); a single-operator batch writes _<idx>.dat.
                if num_op > 1:
                    name = "{}_DynamicalGreen_{}_{}.dat".format(self.header, idx, op)
                else:
                    name = "{}_DynamicalGreen_{}.dat".format(self.header, idx)
                path = os.path.join(spec_dir, name)
                d = np.loadtxt(path)
                raw[(op, idx)] = d[:, 2] + 1J * d[:, 3]
                if op == 0 and idx == 0:
                    frequencies = d[:, 0] + 1J * d[:, 1]
        finite = {op: {} for op in range(num_op)}
        for T in self.T_list:
            Z = self._calc_Z(T)
            weights = [np.exp(-(self.energy_list[idx] - self.ene_min) / T) for idx in range(self.exct)]
            for op in range(num_op):
                spectrum = np.zeros_like(raw[(op, 0)])
                for idx in range(self.exct):
                    spectrum += weights[idx] * raw[(op, idx)]
                finite[op][T] = spectrum / Z
        return frequencies, finite

    def get_one_body_green_batch(self, ops, exct_cut, batch_id=0):
        """Evaluate a batch of single-excitation operators in ONE HPhi run.

        ops: list of (sitei, sigmai, sitej, sigmaj, ex_state, flg), all sharing the
        same ex_state and the same Hilbert sector (verified by HPhi at run time).
        The eigenvector for each eigenstate is read once and reused across all
        operators (op-inner). Returns a list of one_body_green arrays aligned with ops.
        """
        ex_state = ops[0][4]
        num_op = len(ops)
        calc_dir = os.path.join(self.parent_dir, "batch_{}".format(batch_id))
        os.makedirs(calc_dir, exist_ok=True)
        self.Make_Spectrum_Input(calc_dir)
        # operator set 0 -> single_ex.def (the namelist SingleExcitation), sets 1.. -> single_ex_<op>.def
        for op, (si, sgi, sj, sgj, exs, flg) in enumerate(ops):
            fname = "single_ex.def" if op == 0 else "single_ex_{}.def".format(op)
            self._make_single_excitation(si, sgi, sj, sgj, file_name=fname,
                                         ex_state=exs, flg_complex=flg, calc_dir=calc_dir)
        self._run_HPhi(exct_cut, ex_state, calc_dir, num_op=num_op)
        frequencies, finite = self._finite_T_spectrum_multiop(calc_dir, num_op)
        if ex_state == 1:
            self.frequencies = frequencies
        sign = 1.0 if ex_state == 1 else -1.0
        out = []
        for op in range(num_op):
            one_body_green = np.zeros((len(self.T_list), self.nomega), dtype=np.complex128)
            for idx, T in enumerate(self.T_list):
                one_body_green[idx] = sign * finite[op][T]
            out.append(one_body_green)
        return out

    def _finite_T_spectrum_braket(self, calc_dir, num_op, num_bra):
        """Boltzmann-sum the per-(idx, op, bra) spectra written by a SpectrumNumBra run.

        Returns (frequencies, {(op, bra, T): spectrum}). The file name matches HPhi's naming,
        which depends on what is active: with num_bra > 1 it is
        <header>_DynamicalGreen_<idx>_<op>_<bra>.dat (the _<op> field is forced on); with a single
        operator AND single bra (num_op == num_bra == 1, e.g. a single-orbital impurity) HPhi
        falls back to the unsuffixed <header>_DynamicalGreen_<idx>.dat.
        """
        spec_dir = os.path.join(calc_dir, self.output_dir)
        use_bra = num_bra > 1
        use_op = use_bra or num_op > 1
        frequencies = None
        raw = {}  # (op, bra, idx) -> complex spectrum
        for op in range(num_op):
            for b in range(num_bra):
                for idx in range(self.exct):
                    if use_bra:
                        name = "{}_DynamicalGreen_{}_{}_{}.dat".format(self.header, idx, op, b)
                    elif use_op:
                        name = "{}_DynamicalGreen_{}_{}.dat".format(self.header, idx, op)
                    else:
                        name = "{}_DynamicalGreen_{}.dat".format(self.header, idx)
                    d = np.loadtxt(os.path.join(spec_dir, name))
                    raw[(op, b, idx)] = d[:, 2] + 1J * d[:, 3]
                    if op == 0 and b == 0 and idx == 0:
                        frequencies = d[:, 0] + 1J * d[:, 1]
        finite = {}
        for T in self.T_list:
            Z = self._calc_Z(T)
            weights = [np.exp(-(self.energy_list[idx] - self.ene_min) / T) for idx in range(self.exct)]
            for op in range(num_op):
                for b in range(num_bra):
                    spectrum = np.zeros_like(raw[(op, b, 0)])
                    for idx in range(self.exct):
                        spectrum += weights[idx] * raw[(op, b, idx)]
                    finite[(op, b, T)] = spectrum / Z
        return frequencies, finite

    def get_one_body_green_braket_batch(self, ex_state, ket_ops, bra_ops, exct_cut, batch_id=0):
        """Stage-3 bra/ket reuse: ONE HPhi run that solves the resolvent for every ket
        c_{j,sigma_j} (SpectrumNumOp) and projects each solve onto every bra c_{i,sigma_i}
        (SpectrumNumBra), giving G_{i,j} = <c_i phi|(z-(H-E))^{-1}|c_j phi> directly -- no
        c_i + i c_j combination trick. This cuts the BiCG count from n_orb^2 to n_orb.

        ket_ops / bra_ops: lists of (site, sigma). For a canonical (Sz-resolved) sector both
        lists are the sites of ONE spin; for grand canonical they are all 2*n_site spin-orbitals
        (kets and bras of both spins share the Ne+-1 excited space). All ket/bra operators must
        map to the same excited sector, which HPhi verifies (skipped for grand-canonical models).

        Returns {(row_op, col_op): one_body_green[(n_T, n_omega)]} where (row_op, col_op) is the
        storage index of g = G_{bra_op, ket_op} given by _braket_store_op (which encodes the
        creation-channel transpose); the caller accumulates one_body_green[row_op][col_op] += g.
        """
        num_op = len(ket_ops)
        num_bra = len(bra_ops)
        calc_dir = os.path.join(self.parent_dir, "braket_{}".format(batch_id))
        os.makedirs(calc_dir, exist_ok=True)
        # The bra/ket projection rides on the SpectrumLoopExct internal eigenstate loop.
        self.internal_loop = True
        self.Make_Spectrum_Input(calc_dir, bra=True)
        # kets: op 0 -> single_ex.def (the namelist SingleExcitation), ops 1.. -> single_ex_<k>.def
        for k, (sj, sgj) in enumerate(ket_ops):
            fname = "single_ex.def" if k == 0 else "single_ex_{}.def".format(k)
            self._make_single_excitation(sj, sgj, sj, sgj, file_name=fname,
                                         ex_state=ex_state, flg_complex=True, calc_dir=calc_dir)
        # bras: bra b -> single_ex_bra_<b>.def (set 0 is the namelist SingleExcitationBra)
        for b, (si, sgi) in enumerate(bra_ops):
            self._make_single_excitation(si, sgi, si, sgi, file_name="single_ex_bra_{}.def".format(b),
                                         ex_state=ex_state, flg_complex=True, calc_dir=calc_dir)
        self._run_HPhi(exct_cut, ex_state, calc_dir, num_op=num_op, num_bra=num_bra)
        frequencies, finite = self._finite_T_spectrum_braket(calc_dir, num_op, num_bra)
        if ex_state == 1:
            self.frequencies = frequencies
        sign = 1.0 if ex_state == 1 else -1.0
        out = {}
        for k, ket_op in enumerate(ket_ops):
            for b, bra_op in enumerate(bra_ops):
                g = np.zeros((len(self.T_list), self.nomega), dtype=np.complex128)
                for t_idx, T in enumerate(self.T_list):
                    g[t_idx] = sign * finite[(k, b, T)]
                # Key by the storage index (which encodes the creation-channel transpose).
                out[_braket_store_op(bra_op, ket_op, ex_state)] = g
        return out


def test_main():
    args = sys.argv
    if len(args) != 2:
        print("Error: Wrong argument.")
        print("Usage: python hphi_spectrum.py filename")
        exit(1)

    file_name = sys.argv[1]
    import toml
    dict_toml = toml.load(open(file_name))
    NOmega = 200
    T_list = dict_toml.get("T_list", [1.0])
    exct = dict_toml.get("exct", 10)
    eta = dict_toml.get("eta", 1e-4)
    path_to_HPhi = dict_toml.get("path_to_HPhi", "./HPhi")
    header = dict_toml.get("header", "zvo")
    output_dir = dict_toml.get("output_dir", "./output")
    n_site = dict_toml.get("n_site", 2)
    max_workers=4

    #Calculate one body Green's functions using parallel
    n_sigma = 2
    n_flg = 2
    n_excitation = 2
    p_common = (n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct)
    one_body_green = calc_one_body_green_core_parallel(p_common)
    np.save("test_g", one_body_green)

if __name__ == "__main__":
    test_main()
