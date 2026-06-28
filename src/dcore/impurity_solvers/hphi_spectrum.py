import subprocess
import itertools
import numpy as np
import os
import sys

# T_list, n_iw, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output"

def calc_one_body_green_core_parallel(p_common, max_workers=None, sector_occupancy=None):
    """
    Return:
        np.ndarray(n_site, n_sigma, n_site, n_sigma, n_T, n_omega)

    sector_occupancy: optional (n_up, n_down, n_site) for the canonical (sector-resolved)
        path. When given, operators with an exactly-zero contribution in that (Ne, 2Sz)
        sector -- cross-spin off-diagonal terms (zero by spin conservation), and
        annihilation/creation on an already empty/full spin -- are skipped instead of run,
        which also avoids HPhi excitations into non-existent sectors. Their contribution is
        filled with zeros. None (default) is the grand-canonical path: every operator runs.
    """

    n_sigma = 2
    n_flg = 2
    n_excitation = 2
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut, *rest = p_common

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

    def Make_Spectrum_Input(self, calc_dir="./", spectrum_type="single"):

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
                                    "SingleExcitation", "PairExcitation"]:
                        continue
                    fex.write("{} {}\n".format(words[0], os.path.join(rel_path_org, words[1])))
                fex.write("ModPara modpara_ex.def\n")
                fex.write("CalcMod calcmod_ex.def\n")
                fex.write("SpectrumVec    {}\n".format(spectrum_vec))
                if spectrum_type == "single":
                    fex.write("SingleExcitation single_ex.def\n")
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

    def _update_modpara(self, exct, ex_state=0, calc_dir="./", num_op=1):
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

    def _run_HPhi(self, exct_cut, ex_state=0, calc_dir="./", num_op=1):
        os.chdir(calc_dir)
        exec_path = self.path_to_HPhi
        if self.internal_loop:
            # One launch covers all exct_cut eigenstates (and all num_op operator sets);
            # HPhi writes output/<header>_DynamicalGreen_<idx>[_<op>].dat directly.
            self._update_modpara(exct_cut, ex_state, calc_dir, num_op=num_op)
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
