#
# DCore -- Integrated DMFT software for correlated electrons
# Copyright (C) 2017 The University of Tokyo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import numpy
from itertools import product
import os
import sys
from collections import namedtuple
import shlex
import math

from dcore._dispatcher import *

from ..tools import make_block_gf, launch_mpi_subprocesses, extract_H0, extract_bath_params, expand_path
from .base import SolverBase
from .hphi_spectrum import calc_one_body_green_core_parallel
from .pomerol import assign_from_numpy_array


namelist_def = """\
ModPara  modpara.def
CalcMod  calcmod.def
LocSpin  locspn.def
Trans  trans.def
InterAll  interall.def
"""


modpara_def = """\
--------------------
Model_Parameters   0
--------------------
HPhi_Cal_Parameters
--------------------
CDataFileHead  zvo
CParaFileHead  zqp
--------------------
Nsite          {0}
Lanczos_max    2000
initial_iv     -1
exct           {1}
LanczosEps     14
LanczosTarget  2
LargeValue     4.000000000000000e+00
NumAve         5
ExpecInterval  20
NOmega         {2}
OmegaMax       0.0     {3}
OmegaMin       0.0     {4}
OmegaOrg       0.0     0.0
"""


calcmod_def = """\
#CalcType = 0:Lanczos, 1:TPQCalc, 2:FullDiag, 3:CG, 4:Time-evolution
#CalcModel = 0:Hubbard, 1:Spin, 2:Kondo, 3:HubbardGC, 4:SpinGC, 5:KondoGC
#Restart = 0:None, 1:Save, 2:Restart&Save, 3:Restart
#CalcSpec = 0:None, 1:Normal, 2:No H*Phi, 3:Save, 4:Restart, 5:Restart&Save
CalcType   3
CalcModel   3
ReStart   0
CalcSpec   0
CalcEigenVec   0
InitialVecType   0
InputEigenVec   0
OutputEigenVec   1
InputHam   0
OutputHam   0
OutputExVec   0
"""


locspn_def = """\
================================
NlocalSpin     0
================================
========i_0LocSpn_1IteElc ======
================================
"""


trans_def = """\
========================
NTransfer      {0}
========================
========i_j_s_tijs======
========================
"""


interall_def = """\
======================
NInterAll      {0}
======================
========zInterAll=====
======================
"""


def read_eigenenergies(energy_file):
    """Read the eigenenergies (one per computed state) from an HPhi zvo_energy.dat file."""
    energies = []
    with open(energy_file) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) >= 2 and tokens[0] == 'Energy':
                energies.append(float(tokens[1]))
    return numpy.array(energies)


def warn_if_exct_truncates_thermal_trace(energy_file, beta, exct, exct_max, weight_threshold=1e-3):
    """
    Warn if too few eigenstates were computed to span the thermally relevant
    multiplet at temperature T = 1/beta.

    HPhi builds the finite-T Green's function from the lowest ``exct`` eigenstates,
    weighting state n by the Boltzmann factor exp(-beta*(E_n - E_0)).  If the
    highest computed state still carries a non-negligible weight, states just
    above the cutoff are missing from the thermal trace.  This typically happens
    when the ground state is degenerate (common in multi-orbital models) and the
    default ``exct = 1`` keeps only one member of the multiplet -- silently
    breaking orbital symmetry and producing a wrong (e.g. spurious off-diagonal)
    self-energy.

    A warning is printed (it does not abort the solver).  It never fires when the
    full Hilbert space is already covered (exct == exct_max).
    """
    if exct >= exct_max:
        return  # full Hilbert space is covered; nothing is truncated

    try:
        energies = read_eigenenergies(energy_file)
    except (OSError, ValueError):
        return  # tolerate a missing/odd energy file rather than aborting the solver
    if energies.size == 0:
        return

    e0 = energies.min()
    e_last = energies.max()
    # eigenenergies are printed to ~1e-9 precision, so a loose tol identifies degeneracy
    degeneracy = int(numpy.count_nonzero(numpy.abs(energies - e0) < 1e-6))

    # All omitted states have energy >= e_last (we kept the lowest exct states),
    # so w_last = exp(-beta*(e_last - e0)) is an UPPER BOUND on the Boltzmann
    # weight of each omitted state. If w_last <= threshold the omitted tail is
    # provably negligible; otherwise truncation cannot be ruled out (we cannot
    # see whether the next state sits just above e_last or far above it), so we
    # warn conservatively.
    w_last = math.exp(-beta * (e_last - e0))

    if w_last > weight_threshold:
        if degeneracy >= exct:
            deg_msg = (f"  All {exct} computed states are degenerate with the ground state, "
                       f"so the cutoff falls inside the ground multiplet.\n")
        else:
            deg_msg = f"  The ground state is {degeneracy}-fold degenerate.\n"
        print(
            "\n*** WARNING (HPhi solver): 'exct' may be too small ***\n"
            f"  exct = {exct} eigenstates were computed (full space = {exct_max}).\n"
            + deg_msg +
            f"  The highest computed state has Boltzmann weight {w_last:.2e} at "
            f"T = {1.0 / beta:.4g}; this is an upper bound on the weight of every\n"
            f"  omitted (higher) state, and it exceeds {weight_threshold:.0e}. A "
            "thermally-relevant state above the cutoff therefore cannot be ruled\n"
            "  out: the finite-T trace may be truncated, which can break orbital "
            "symmetry and yield a wrong (e.g. spurious off-diagonal) self-energy.\n"
            "  => Increase 'exct'; if the result does not change, it was already converged.\n",
            file=sys.stderr,
        )


def enumerate_particle_sectors(n_site):
    """Enumerate the (Ne, 2Sz) sectors of a system of ``n_site`` sites (2*n_site spin-orbitals).

    A sector is fixed in HPhi by ``Ncond = Ne`` and ``2Sz``; here ``n_up = (Ne + 2Sz)//2`` and
    ``n_down = (Ne - 2Sz)//2``. The grand-canonical Hilbert space (dimension ``4**n_site``) is
    the disjoint union of these sectors as a counting fact (every Fock state has a definite
    (Ne, 2Sz)), which the unit tests check via ``sum(dim) == 4**n_site``.

    **Validity of the (Ne, 2Sz) block decomposition for the spectrum.** Replacing one
    grand-canonical run with per-sector canonical runs is correct only when the generated HPhi
    Hamiltonian is block-diagonal in (Ne, 2Sz) -- i.e. it conserves BOTH the total electron
    number Ne and 2Sz. That holds for density-density interactions plus number-conserving,
    spin-diagonal one-body terms. It is BROKEN by anomalous/pairing terms (Ne not conserved) or
    by spin-mixing one-body terms such as spin-orbit coupling (the transfer block can carry
    ``s1 != s2``), which conserve Ne but not 2Sz. When only Ne is conserved one must sector by
    Ne alone (no 2Sz constraint); when neither is conserved the canonical-by-sector path does
    not apply and the grand-canonical run must be used. The caller (per-sector solver loop) is
    responsible for checking these conservation laws on the actual Trans/InterAll before using
    this enumeration. When the decomposition is valid, a sector's eigenvalues are exactly the
    grand-canonical eigenvalues living in it, so the finite-T trace equals the sum of the
    per-sector contributions (with a single global E_min and partition function Z).

    Returns a list of dicts ``{'Ne', 'two_Sz', 'n_up', 'n_down', 'dim'}`` sorted by Ne then 2Sz.
    """
    from math import comb
    if n_site < 0:
        raise ValueError("n_site must be non-negative, got {}".format(n_site))
    sectors = []
    for n_up in range(n_site + 1):
        for n_down in range(n_site + 1):
            sectors.append({
                'Ne': n_up + n_down,
                'two_Sz': n_up - n_down,
                'n_up': n_up,
                'n_down': n_down,
                'dim': comb(n_site, n_up) * comb(n_site, n_down),
            })
    sectors.sort(key=lambda s: (s['Ne'], s['two_Sz']))
    return sectors


def gf_parallel_layout(mpirun_command, np_total, n_procs_per_hphi):
    """
    Decide the two-level parallel layout for the Gf (one-body Green's function) step.

    The Gf step runs one HPhi calculation per (excitation) task; the tasks are
    independent, so they can be spread over ``n_outer`` concurrent runs, each of
    which may itself use ``n_inner`` MPI ranks.

    Parameters
    ----------
    mpirun_command : str
        The MPI launcher, e.g. ``"mpirun -np 16"`` (last token = number of ranks).
    np_total : int or None
        Total number of processes (``None`` if it could not be parsed).
    n_procs_per_hphi : int
        Requested MPI ranks per HPhi run (``n_inner``). ``1`` keeps the previous
        behaviour (serial HPhi, ``np_total`` concurrent runs).

    Returns
    -------
    (n_inner, n_outer, hphi_mpi_command, gf_mpi_prefix)
        ``n_inner`` is clamped to a power of four not exceeding ``np_total``;
        ``n_outer = np_total // n_inner`` is the ProcessPool size.
        ``hphi_mpi_command`` is the launcher for a *single* HPhi run with
        ``n_inner`` ranks -- it is used for BOTH the eigenvalue step and each
        Green's-function run. They must use the same rank count because the
        eigenvectors are written MPI-distributed (one file per rank) by the
        eigenvalue step and read back by the Green's-function step.
        ``gf_mpi_prefix`` is that same launcher prepended to each Gf shell
        command (``""`` when ``n_inner == 1``, i.e. a bare/serial HPhi run).
    """
    if np_total is None:
        # The rank count is not the launcher's last token, so we cannot rewrite
        # it to split the Gf step. Preserve correctness instead of speed: drive
        # BOTH phases with the *same* launcher (so the eigenvector rank counts
        # match) and do not run Gf runs concurrently (n_outer = 1), since we
        # cannot tell how many ranks each run would take.
        print("Note (HPhi solver): could not parse the process count from the MPI "
              "command; the Gf step is not parallelized (put the rank count as the "
              "last token of [mpi] command to enable it).", file=sys.stderr)
        return 1, 1, mpirun_command, mpirun_command
    n_inner = max(1, min(int(n_procs_per_hphi), np_total))
    if not math.log(n_inner, 4).is_integer():  # HPhi requires a power of four
        n_inner = 4 ** int(math.log(n_inner, 4))
        print(f"Warning: n_procs_per_hphi must be a power of four in HPhi. "
              f"It is set to {n_inner} (ranks per HPhi run).", file=sys.stderr)
    n_outer = max(1, np_total // n_inner)
    cmds = shlex.split(mpirun_command)
    cmds[-1] = str(n_inner)
    hphi_mpi_command = " ".join(cmds)
    # The Gf runs may go bare when single-rank, but the eigenvalue step always
    # uses the launcher so its rank count matches what the Gf step expects.
    gf_mpi_prefix = "" if n_inner == 1 else hphi_mpi_command
    return n_inner, n_outer, hphi_mpi_command, gf_mpi_prefix


class HPhiSolver(SolverBase):

    def __init__(self, beta, gf_struct, u_mat, n_iw=1025):
        """
        Initialize the solver.
        """

        super(HPhiSolver, self).__init__(beta, gf_struct, u_mat, n_iw)

    def solve(self, rot, mpirun_command, params_kw):
        """
        In addition to the parameters described in the docstring of SolverBase,
        one can pass solver-dependent parameters using params_kw. For example,
          exec_path : str, path to an executable, mandatory
          dry_run   : bool, actual computation is not performed if dry_run is True, optional
        """

        # (1) Set configuration for the impurity solver
        # input:
        #   self.beta
        #   self.set_G0_iw
        #   self.u_mat
        #
        # Additionally, the following variables may be used:
        #   self.n_orb
        #   self.n_flavor
        #   self.gf_struct
        #   self.use_spin_orbit

        exec_path = expand_path(params_kw['exec_path'])

        # Parse the total number of processes from the launcher's last token.
        commands = shlex.split(mpirun_command)
        try:
            np = int(commands[-1])
        except ValueError:
            np = None
            print("A check of np is skipped.")

        # Two-level parallelism for the Gf step. CRUCIAL: the eigenvalue step and
        # every Green's-function run must use the SAME number of MPI ranks
        # (n_inner), because the eigenvalue step writes the eigenvectors
        # MPI-distributed (one file per rank) and the Gf step reads them back -- a
        # rank-count mismatch makes HPhi stop while inputting the eigenvector.
        #   n_procs_per_hphi (= n_inner): MPI ranks per HPhi run (power of four;
        #   1 = serial, the default). n_outer = np // n_inner runs run concurrently.
        n_inner, n_outer, mpirun_command_eigen, gf_mpi_prefix = gf_parallel_layout(
            mpirun_command, np, params_kw.get('n_procs_per_hphi', 1))

        # Matsubara frequencies omega_n = (2*n+1)*pi*T
        omega_min = numpy.pi / self.beta  # n=0
        omega_max = (2*self.n_iw + 1) * numpy.pi / self.beta  # n=n_iw
        # NOTE: omega_max is NOT included in the omega mesh.
        #           omega_n = (omega_max - omega_min) / n_iw * n
        #       for n=[0:n_iw)

        # bath fitting
        n_bath = params_kw.get('n_bath', 0)  # 0 for Hubbard-I approximation
        exct = params_kw.get('exct', 1)  # number of states to be computed
        # Boltzmann-weight threshold above which an under-sized exct triggers a warning
        exct_weight_threshold = params_kw.get('exct_weight_threshold', 1e-3)

        fit_params = {}
        for key in ['fit_gtol',]:
            if key in params_kw:
                fit_params[key] = params_kw[key]

        n_site = self.n_orb + n_bath

        exct_max = 4**n_site
        if exct > exct_max:
            print(f"Warning: exct={exct} is larger than {exct_max}. exct is set to {exct_max}", file=sys.stderr)
            exct = exct_max

        # Output namelist.def
        with open('./namelist.def', 'w') as f:
            print(namelist_def, end="", file=f)

        # Output modpara.def
        with open('./modpara.def', 'w') as f:
            print(modpara_def.format(n_site, exct, self.n_iw, omega_max, omega_min), end="", file=f)

        # Output calcmod.def
        with open('./calcmod.def', 'w') as f:
            print(calcmod_def, end="", file=f)

        # Output locspn.def
        with open('./locspn.def', 'w') as f:
            print(locspn_def, end="", file=f)

            for i in range(n_site):
                print("{0} 0".format(i), file=f)

        # -------------------------------------------------------------------------

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
        h0_mat = extract_H0(self._G0_iw, self.block_names)
        assert h0_mat.shape == (self.n_flavors, self.n_flavors)

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(self._G0_iw)

        bath_levels, bath_hyb = extract_bath_params(self._Delta_iw, self.beta, self.block_names, n_bath, **fit_params)
        assert bath_levels.shape == (2*n_bath,)
        assert bath_hyb.shape == (self.n_flavors, 2*n_bath)

        # make hopping matrix
        Transfer = namedtuple('Transfer', ('i1', 's1', 'i2', 's2', 't'))
        transfer = []

        # A. correlated sites
        h0_isjs = h0_mat.reshape((2, self.n_orb, 2, self.n_orb))
        for s1, s2 in product(range(2), repeat=2):
            for i1, i2 in product(range(self.n_orb), repeat=2):
                t = -h0_isjs[s1, i1, s2, i2]
                # print(i1, s1, i2, s2, t.real, t.imag, file=f)
                if t != 0:
                    transfer.append(Transfer(i1, s1, i2, s2, t))

        # B. bath levels
        bath_levels_is = bath_levels.reshape((2, n_bath))
        for s1 in range(2):
            for i1 in range(n_bath):
                j1 = self.n_orb + i1
                eps = -bath_levels_is[s1, i1]
                # print(j1, s1, j1, s1, eps.real, eps.imag, file=f)
                if eps != 0:
                    transfer.append(Transfer(j1, s1, j1, s1, eps))

        # C. hopping between correlated sites and bath sites
        bath_hyb_is = bath_hyb.reshape((2, self.n_orb, 2, n_bath))
        for s1 in range(2):
            for i1 in range(self.n_orb):
                for s2 in range(2):
                    for i2 in range(n_bath):
                        j1 = i1
                        j2 = self.n_orb + i2
                        v = -bath_hyb_is[s1, i1, s2, i2]
                        # print(j1, s1, j2, s2, v.real, v.imag, file=f)
                        # print(j2, s2, j1, s1, v.real, -v.imag, file=f)
                        if v != 0:
                            transfer.append(Transfer(j1, s1, j2, s2, v))
                            transfer.append(Transfer(j2, s2, j1, s1, numpy.conj(v)))

        # Output trans.def
        with open('./trans.def', 'w') as f:
            print(trans_def.format(len(transfer)), end="", file=f)

            for t in transfer:
                print(t.i1, t.s1, t.i2, t.s2, t.t.real, t.t.imag, file=f)

        # -------------------------------------------------------------------------

        # (1c) Set U_{ijkl} for the solver
        # for i, j, k, l in product(range(self.n_flavors), repeat=4):
        #     self.u_mat[i, j, k, l]

        # make U matrix
        InterAll = namedtuple('InterAll', ('i1', 's1', 'i2', 's2', 'i3', 's3', 'i4', 's4', 'U'))
        interall = []

        # (1/2) U_{1234} c_1^+ c_2^+ c_4 c_3  # Dcore
        # = I_{1324} c_1^+ c_3 c_2^+ c_4      # Hphi
        u_1234 = self.u_mat.reshape((2, self.n_orb, 2, self.n_orb, 2, self.n_orb, 2, self.n_orb))
        for s1, s2, s3, s4 in product(range(2), repeat=4):
            for o1, o2, o3, o4 in product(range(self.n_orb), repeat=4):
                u = u_1234[s1, o1, s2, o2, s3, o3, s4, o4] / 2.
                if s1==s2==s3==s4 and o1==o2==o3==o4:
                    continue
                    # u = 0.0
                # print(o1, s1, o3, s3, o2, s2, o4, s4, u.real, u.imag, file=f)
                if numpy.abs(u) > 1e-10:
                    interall.append(InterAll(o1, s1, o3, s3, o2, s2, o4, s4, u))

        # Output interall.def
        with open('./interall.def', 'w') as f:
            print(interall_def.format(len(interall)), end="", file=f)

            for u in interall:
                print(u.i1, u.s1, u.i2, u.s2, u.i3, u.s3, u.i4, u.s4, u.U.real, u.U.imag, file=f)

        # (2) Run a working horse
        # Sector-resolved (canonical) path is valid only when H conserves both Ne and 2Sz.
        # Ne is always conserved (only c^dag c hopping + density interactions); 2Sz is broken
        # by spin-orbit coupling (spin-mixing one-body terms), so guard on use_spin_orbit.
        use_canonical = (os.environ.get("DCORE_HPHI_CANONICAL_SECTORS", "0") == "1"
                         and not self.use_spin_orbit)

        if use_canonical:
            # H is block-diagonal in (Ne, 2Sz), so the grand-canonical finite-T trace equals the
            # weighted sum of per-sector canonical contributions. Each sector spans a much smaller
            # Hilbert space than the full 4**n_site grand-canonical space, which is where the
            # dominant Green's-function (BiCG) cost shrinks. Each sector run returns its own
            # finite-T-normalized Gf G_s and we recombine exactly with
            #   G = sum_s w_s G_s / sum_s w_s,  w_s = Z_s * exp(-beta (Emin_s - Emin_global)),
            # which reproduces the single global Boltzmann sum.
            print("\nComputing eigenenergies + Gf per (Ne, 2Sz) sector ...")
            canonical_calcmod = calcmod_def.replace("CalcModel   3", "CalcModel   0")
            T_list = [1. / self.beta]
            eta = 1e-4
            # Drop sectors whose (bounded) thermal contribution is below this, relative to the
            # global ground state: they contribute negligibly to the finite-T trace.
            # 0 disables the selection (compute every sector -- exact but slower).
            sector_weight_threshold = params_kw.get('canonical_sector_weight_threshold', 1e-8)
            if not (0.0 <= sector_weight_threshold < 1.0):
                raise ValueError("canonical_sector_weight_threshold must be in [0, 1), got {}"
                                 .format(sector_weight_threshold))

            def run_sector_eigenvalues(Ne, two_Sz, exct_use):
                """Write the canonical (Ncond, 2Sz) modpara/calcmod and run the eigenvalue step;
                return the eigenenergies (empty array if the sector produced nothing)."""
                sector_modpara = modpara_def.format(n_site, exct_use, self.n_iw, omega_max, omega_min)
                # the convergence target cannot exceed the number of states in this sector
                sector_modpara = sector_modpara.replace("LanczosTarget  2",
                                                        "LanczosTarget  {}".format(min(2, exct_use)))
                sector_modpara += "Ncond          {}\n2Sz            {}\n".format(Ne, two_Sz)
                with open('./modpara.def', 'w') as f:
                    f.write(sector_modpara)
                with open('./calcmod.def', 'w') as f:
                    f.write(canonical_calcmod)
                with open('./stdout.log', 'w') as output_f:
                    launch_mpi_subprocesses(mpirun_command_eigen, [exec_path, '-e', 'namelist.def'], output_f)
                return read_eigenenergies(os.path.join('output', 'zvo_energy.dat'))

            all_sectors = [s for s in enumerate_particle_sectors(n_site) if min(exct, s['dim']) >= 1]

            # Phase-2 thermal selection: a cheap ground-state-only pre-pass per sector picks the
            # thermally relevant (Ne, 2Sz) sectors, so the dominant per-sector Gf step runs only
            # for those. The spectrum reaches the (Ne+-1) excited sectors internally, so only the
            # thermally OCCUPIED sectors need an eigenvalue/Gf run here.
            if sector_weight_threshold > 0.0:
                gs = []
                for sec in all_sectors:
                    e = run_sector_eigenvalues(sec['Ne'], sec['two_Sz'], 1)
                    if e.size > 0:
                        gs.append((sec, float(e.min())))
                if not gs:
                    raise RuntimeError("No (Ne, 2Sz) sector produced eigenstates.")
                e_gs_global = min(e for _, e in gs)
                # A sector's contribution to the un-normalized trace is
                # Z_sec * exp(-beta(Egs_sec - Egs_global)) with Z_sec = sum_i exp(-beta(E_i-Egs_sec))
                # <= the number of states summed (<= min(exct, dim)). Bounding Z_sec by that count
                # (rather than filtering on the ground-state weight alone) makes the cut safe even
                # for sectors whose many low-lying / nearly degenerate states give a large Z_sec.
                sectors = [sec for sec, e in gs
                           if min(exct, sec['dim']) * numpy.exp(-self.beta * (e - e_gs_global))
                           > sector_weight_threshold]
                print("  thermal selection: {} of {} sectors kept (weight > {:.1e})".format(
                    len(sectors), len(all_sectors), sector_weight_threshold), flush=True)
            else:
                sectors = all_sectors

            contributions = []  # list of (G_sector, E_min_sector, Z_sector)
            for sec in sectors:
                Ne, two_Sz, dim = sec['Ne'], sec['two_Sz'], sec['dim']
                exct_sec = min(exct, dim)
                energies = run_sector_eigenvalues(Ne, two_Sz, exct_sec)
                if energies.size == 0:
                    continue
                # Same thermal-truncation guard as the grand-canonical path, per sector: if
                # exct_sec < dim and the highest retained state still carries weight, this
                # sector's trace (and Z_sec) is truncated, which would bias the recombined Gf.
                warn_if_exct_truncates_thermal_trace(
                    os.path.join('output', 'zvo_energy.dat'), self.beta, exct_sec, dim,
                    weight_threshold=exct_weight_threshold)
                E_min_sec = float(energies.min())
                Z_sec = float(numpy.sum(numpy.exp(-self.beta * (energies - E_min_sec))))
                p_common = (self.n_orb, T_list, exct_sec, eta, exec_path, "zvo", "./output",
                            exct_sec, gf_mpi_prefix)
                G_sec = calc_one_body_green_core_parallel(
                    p_common, max_workers=n_outer,
                    sector_occupancy=(sec['n_up'], sec['n_down'], n_site))
                contributions.append((G_sec, E_min_sec, Z_sec))
                print("  sector Ne={:2d} 2Sz={:+d} dim={:5d} exct={:3d} E_min={:.6g} Z={:.3g}".format(
                    Ne, two_Sz, dim, exct_sec, E_min_sec, Z_sec), flush=True)
            if not contributions:
                raise RuntimeError("No (Ne, 2Sz) sector produced eigenstates.")
            E_min_global = min(c[1] for c in contributions)
            weights = [c[2] * numpy.exp(-self.beta * (c[1] - E_min_global)) for c in contributions]
            Z_global = sum(weights)
            one_body_g = sum(w * c[0] for w, c in zip(weights, contributions)) / Z_global
            print("\nFinish Gf calc ({} sectors).".format(len(contributions)))
        else:
            print("\nComputing eigeneneries ...")
            with open('./stdout.log', 'w') as output_f:
                launch_mpi_subprocesses(mpirun_command_eigen, [exec_path, '-e', 'namelist.def'], output_f)

            # Warn if too few eigenstates were computed to span the thermally relevant
            # multiplet (e.g. a degenerate ground state with the default exct=1).
            warn_if_exct_truncates_thermal_trace(
                os.path.join('output', 'zvo_energy.dat'), self.beta, exct, exct_max,
                weight_threshold=exct_weight_threshold)

            print("\nComputing Gf ...")
            if n_inner > 1:
                print(f"  Gf parallel layout: {n_outer} concurrent HPhi run(s) x {n_inner} MPI rank(s) each")
            header = "zvo"
            T_list = [1./self.beta]
            eta = 1e-4
            output_dir = "./output"
            p_common = (self.n_orb, T_list, exct, eta, exec_path, header, output_dir, exct, gf_mpi_prefix)
            one_body_g = calc_one_body_green_core_parallel(p_common, max_workers=n_outer)

            print("\nFinish Gf calc.")

        # print(one_body_g.shape)
        assert isinstance(one_body_g, numpy.ndarray)
        assert one_body_g.shape == (self.n_orb, 2, self.n_orb, 2, 1, self.n_iw)

        gf = one_body_g[..., 0, :]
        assert gf.shape == (self.n_orb, 2, self.n_orb, 2, self.n_iw)

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

        # Change data structure of gf from [o1, s1, o2, s2, iw] to ...
        if self.use_spin_orbit:
            # [1, (s1,o1), (s2,o2), iw]
            gf = gf.transpose((1, 0, 3, 2, 4)).reshape((1, 2*self.n_orb, 2*self.n_orb, self.n_iw))
            assert gf.shape == (1, 2*self.n_orb, 2*self.n_orb, self.n_iw)
        else:
            # [s, o1, o2, iw]
            gf = numpy.einsum("isjsw->sijw", gf)
            assert gf.shape == (2, self.n_orb, self.n_orb, self.n_iw)

        assign_from_numpy_array(self._Gimp_iw, gf, self.block_names)

        # if triqs_major_version == 1:
        #     set_tail(self._Gimp_iw)

        if self.use_spin_orbit:
            print("Sigma is not implemented for SOC")
            raise NotImplementedError

        # Make H0 matrix
        h0_full = numpy.zeros((2, n_site, 2, n_site), dtype=complex)
        for t in transfer:
            h0_full[t.s1, t.i1, t.s2, t.i2] = -t.t
        h0_full = h0_full.reshape((2*n_site, 2*n_site))

        # TODO: move into a function -- begin
        # Cut H0 into block structure
        n_block = len(self.gf_struct)
        n_inner = h0_full.shape[0] // n_block
        h0_block = [h0_full[s*n_inner:(s+1)*n_inner, s*n_inner:(s+1)*n_inner] for s in range(n_block)]

        # Construct G0 including bath sites
        bath_names = ["bath" + str(i_bath) for i_bath in range(n_bath)]
        bath_names = ["bath" + str(i_bath) for i_bath in range(n_bath)]
        gf_struct_full = {block: list(inner_names) + bath_names for block, inner_names in self.gf_struct.items()}
        g0_full = make_block_gf(GfImFreq, gf_struct_full, self.beta, self.n_iw)
        g0_full << iOmega_n
        for i, block in enumerate(self.block_names):
            g0_full[block] -= h0_block[i]
        g0_full.invert()

        # Project G0 onto impurity site
        g0_imp = make_block_gf(GfImFreq, self.gf_struct, self.beta, self.n_iw)
        for block in self.block_names:
            for o1, o2 in product(self.gf_struct[block], repeat=2):
                g0_imp[block].data[:, o1, o2] = g0_full[block].data[:, o1, o2]
        # TODO: move into a function -- end

        self._Sigma_iw << inverse(g0_imp) - inverse(self._Gimp_iw)

    def name(self):
        return "HPhi"
