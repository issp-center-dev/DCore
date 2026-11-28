import subprocess
import itertools
import numpy as np
import os
import sys

# T_list, n_iw, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output"


def available_cores():
    """Number of CPU cores available to this process"""
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:
        return os.cpu_count() or 1


def omp_env(n_parallel, omp_num_threads=None):
    """
    Environment for HPhi subprocesses with OMP_NUM_THREADS set.

    HPhi is parallelized by OpenMP. When several HPhi processes run at the same time
    (process pool or MPI) without a limit on the number of threads, the cores are
    oversubscribed and the OpenMP synchronization can slow down the calculation by
    orders of magnitude.

    OMP_NUM_THREADS: omp_num_threads (argument) > OMP_NUM_THREADS in the environment
                     > available cores // n_parallel.
    OMP_WAIT_POLICY: PASSIVE unless set in the environment. Idle threads then sleep instead of
                     spin-waiting, which avoids the slowdown when the cores are shared with
                     other jobs (observed: 22 s -> 0.07 s for one spectrum run of a 256-dim. model).
    """
    env = dict(os.environ)
    if omp_num_threads is not None:
        env['OMP_NUM_THREADS'] = str(int(omp_num_threads))
    elif 'OMP_NUM_THREADS' not in env:
        env['OMP_NUM_THREADS'] = str(max(1, available_cores() // max(1, int(n_parallel))))
    env.setdefault('OMP_WAIT_POLICY', 'PASSIVE')
    return env


def calc_one_body_green_core_parallel(p_common, max_workers=None, omp_num_threads=None):
    """
    Return:
        np.ndarray(n_site, n_sigma, n_site, n_sigma, n_T, n_omega)

    max_workers: number of HPhi processes run in parallel (default: number of available cores)
    omp_num_threads: number of OpenMP threads of each HPhi process (see omp_env)
    """

    n_sigma = 2
    n_flg = 2
    n_excitation = 2
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut = p_common

    check_eta(p_common)

    n_workers = max_workers if max_workers is not None else available_cores()
    env = omp_env(n_workers, omp_num_threads)
    print(f"  HPhi spectrum runs: {n_workers} process(es) in parallel, "
          f"OMP_NUM_THREADS={env['OMP_NUM_THREADS']}, OMP_WAIT_POLICY={env['OMP_WAIT_POLICY']}")

    def gen_p():
        for sitei, sigmai in itertools.product(range(n_site), range(n_sigma)):
            for sitej, sigmaj in itertools.product(range(n_site), range(n_sigma)):
                    for idx, flg in enumerate([True, False]):
                        for ex_state in range(n_excitation):
                            yield sitei, sigmai, sitej, sigmaj, idx, ex_state, p_common, env

    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        one_body_g_tmp = np.array(list(executor.map(calc_one_body_green_core, gen_p())))
    n_omega = one_body_g_tmp.shape[2]
    one_body_green_core = one_body_g_tmp.reshape((n_site, n_sigma, n_site, n_sigma, n_flg, n_excitation, len(T_list), n_omega))
    one_body_green = calc_one_body_green(one_body_green_core)

    import shutil
    for sitei, sigmai in itertools.product(range(n_site), range(n_sigma)):
        for sitej, sigmaj in itertools.product(range(n_site), range(n_sigma)):
            for idx in range(n_flg):
                for ex_state in range(n_excitation):
                    dir_path = "{}_{}_{}_{}_{}_{}".format(sitei, sigmai, sitej, sigmaj, ex_state, idx)
                    shutil.rmtree(dir_path)

    return one_body_green

def calc_one_body_green_core(p):
    #unpack parameters
    sitei, sigmai, sitej, sigmaj, i_flg, ex_state, p_common, env = p
    n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct_cut = p_common
    calc_spectrum_core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi=path_to_HPhi, header=header,
                                           output_dir=output_dir, env=env)

    calc_spectrum_core.set_energies()
    flg = True if i_flg == 0 else False
    return calc_spectrum_core.get_one_body_green_core(sitei, sigmai, sitej, sigmaj, ex_state, flg, exct_cut)

def check_eta(p_common):
    _, T_list, exct, eta, path_to_HPhi, header, output_dir, _ = p_common
    calc_spectrum_core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi=path_to_HPhi, header=header,
                                           output_dir=output_dir)
    calc_spectrum_core.set_energies(check_eta=True)

def calc_one_body_green(one_body_green_core):
    """
    Reconstruct G_{i,j}(iw) from the spectra computed by HPhi.

    Parameters
    ----------
    one_body_green_core : ndarray(n_site, n_sigma, n_site, n_sigma, n_flg, n_excitation, n_T, n_omega)
        one_body_green_core[i, si, j, sj, flg, ex] is the spectrum
            G_X = <X^+ R X>
        of the excitation operator X (HPhi applies the coefficients to the ket),
            X = c_{i si} + lam c_{j sj}  for ex = 0 (hole part, R = (iw + H - E_n)^{-1}),
            X = c^+_{i si} + lam c^+_{j sj}  for ex = 1 (particle part, R = (iw - H + E_n)^{-1}),
        with lam = 1j for flg = 0 and lam = 1 for flg = 1 (X = c_{i si} or c^+_{i si} if (i, si) = (j, sj)).
        The sign of the hole part is already included (see CalcSpectrumCore.get_one_body_green_core).

    Returns
    -------
    one_body_green : ndarray(n_site, n_sigma, n_site, n_sigma, n_T, n_omega)
        G_{ij}(iw) = P_ij + H_ij with the particle part P_ij = <c_i R c^+_j>
        and the hole part H_ij = <c^+_j R c_i>.

    Notes
    -----
    With D_re = G_X(lam=1) - G_ii - G_jj and D_im = G_X(lam=1j) - G_ii - G_jj,
        particle: D_re = P_ij + P_ji,  D_im = 1j (P_ij - P_ji)  ->  P_ij = (D_re - 1j D_im) / 2
        hole:     D_re = H_ij + H_ji,  D_im = 1j (H_ji - H_ij)  ->  H_ij = (D_re + 1j D_im) / 2
    The particle and hole parts must therefore be reconstructed separately.
    """
    n_site, n_sigma, n_site2, n_sigma2, n_flg, n_excitation, n_T, n_omega = one_body_green_core.shape
    assert (n_site, n_sigma) == (n_site2, n_sigma2)
    assert n_flg == 2 and n_excitation == 2
    core = one_body_green_core
    one_body_green = np.zeros((n_site, n_sigma, n_site, n_sigma, n_T, n_omega), dtype=np.complex128)
    for sitei, sigmai, sitej, sigmaj in itertools.product(range(n_site), range(n_sigma), range(n_site), range(n_sigma)):
        for ex_state in range(n_excitation):
            g_ii = core[sitei, sigmai, sitei, sigmai, 0, ex_state]
            if (sitei, sigmai) == (sitej, sigmaj):
                one_body_green[sitei, sigmai, sitej, sigmaj] += g_ii
                continue
            g_jj = core[sitej, sigmaj, sitej, sigmaj, 0, ex_state]
            d_re = core[sitei, sigmai, sitej, sigmaj, 1, ex_state] - g_ii - g_jj
            d_im = core[sitei, sigmai, sitej, sigmaj, 0, ex_state] - g_ii - g_jj
            sign = -1.0 if ex_state == 1 else 1.0
            one_body_green[sitei, sigmai, sitej, sigmaj] += (d_re + sign * 1J * d_im) / 2.0
    return one_body_green


class CalcSpectrum:
    def __init__(self, T_list, n_iw, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output"):
        self.T_list = T_list
        self.exct = exct
        self.eta = eta
        self.header = header
        self.output_dir = output_dir
        self.path_to_HPhi = os.path.abspath(path_to_HPhi)
        self.calc_spectrum_core = CalcSpectrumCore(T_list, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output")
        self.nomega = n_iw

    def get_one_body_green_core(self, n_site, exct_cut):
        self.calc_spectrum_core.set_energies()
        n_excitation = 2 # type of excitation operator
        n_flg = 2
        n_sigma = 2
        one_body_green = np.zeros((n_site, n_sigma, n_site, n_sigma, len(self.T_list), self.nomega), dtype=np.complex128)
        one_body_green_core = np.zeros((n_site, n_sigma, n_site, n_sigma, n_flg, n_excitation, len(self.T_list), self.nomega), dtype=np.complex128)

        for sitei, sigmai, sitej, sigmaj  in itertools.product(range(n_site), range(n_sigma), range(n_site), range(n_sigma)):
            for idx, flg in enumerate([True, False]):
                for ex_state in range(n_excitation):
                    one_body_green_core[sitei][sigmai][sitej][sigmaj][idx][ex_state] = self.calc_spectrum_core.get_one_body_green_core(sitei, sigmai, sitej, sigmaj, ex_state, flg, exct_cut)

        one_body_green = calc_one_body_green(one_body_green_core)
        return one_body_green

class CalcSpectrumCore:
    def __init__(self, T_list, exct, eta, path_to_HPhi="./HPhi", header="zvo", output_dir="./output", env=None):
        self.T_list = T_list
        self.env = env  # environment of the HPhi subprocesses (None: inherit)
        self.exct = exct
        self.eta = eta
        self.header = header
        self.output_dir = output_dir
        self.nomega = 0
        self.parent_dir = os.getcwd()
        # self.path_to_HPhi = os.path.join(self.parent_dir, path_to_HPhi)
        self.path_to_HPhi = os.path.abspath(path_to_HPhi)  # converted to full path in DCore

    def Make_Spectrum_Input(self, calc_dir="./", spectrum_type="single"):

        rel_path_org = os.path.relpath(self.parent_dir, calc_dir)
        for idx in range(self.exct):
            with open(os.path.join(self.parent_dir, "calcmod.def")) as f:
                lines = f.readlines()
            with open(os.path.join(calc_dir, "calcmod_ex.def"), "w") as fex:
                for line in lines:
                   words = line.split()
                   if words[0] == "CalcSpec" or words[0] == "OutputExVec" or words[0] == "OutputEigenVec":
                    continue
                   fex.write(line)
                fex.write("CalcSpec    1\n")
            with open(os.path.join(self.parent_dir, "namelist.def")) as f:
                lines = f.readlines()
            with open(os.path.join(calc_dir, "namelist_ex_{}.def".format(idx)), "w") as fex:
                for line in lines:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    if words[0] in ["CalcMod", "SpectrumVec", "ModPara"]:
                        continue
                    if words[0] == "SingleExcitation":
                        continue
                    if words[0] == "PairExcitation":
                        continue
                    fex.write("{} {}\n".format(words[0], os.path.join(rel_path_org, words[1])))
                fex.write("ModPara modpara_ex.def\n")
                fex.write("CalcMod calcmod_ex.def\n")
                rel_path = os.path.relpath(os.path.join(self.parent_dir, "output"), calc_dir)
                fex.write("SpectrumVec    {}_eigenvec_{}\n".format(os.path.join("../", rel_path, self.header), idx))
                if spectrum_type == "single":
                    fex.write("SingleExcitation single_ex.def\n")
                elif spectrum_type == "pair":
                    fex.write("PairExcitation pair_ex.def\n")

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
        self.ene_min = energy_list[0]
        self.ene_max = energy_list[len(energy_list)-1]

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

    def _update_modpara(self, exct, ex_state=0, calc_dir="./"):
        dict_mod={}
        header = []
        with open(os.path.join(self.parent_dir, "modpara.def"), "r") as fr:
            lines = fr.readlines()
            header = lines[:8]
            for line in lines[8:]:
                words = line.split()
                dict_mod[words[0]] = words[1:]
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
                    fw.write("{} {} {} 1.0 0.0\n".format(site_i, sigma_i, ex_state))
                    fw.write("{} {} {} 0.0 1.0\n".format(site_j, sigma_j, ex_state))
                else:
                    fw.write("{} {} {} 1.0 0.0\n".format(site_i, sigma_i, ex_state))
                    fw.write("{} {} {} 1.0 0.0\n".format(site_j, sigma_j, ex_state))

    def _run_HPhi(self, exct_cut, ex_state=0, calc_dir="./"):
        os.chdir(calc_dir)
        for idx in range(exct_cut):
            self._update_modpara(idx, ex_state, calc_dir)
            input_path = os.path.join(calc_dir, "namelist_ex_{}.def".format(idx))
            exec_path = self.path_to_HPhi
            cmd = "{} -e {} > std_{}.log".format(exec_path, input_path, idx)
            subprocess.call(cmd, shell=True, env=self.env)
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

    #Calculate one body Green's functions directly
    # calcg = CalcSpectrum(T_list, NOmega, exct, eta, path_to_HPhi, header, output_dir)
    # one_body_g_direct = calcg.get_one_body_green(n_site, exct)
    # print(one_body_g_direct)
    # exit(0)

    #Calculate one body Green's functions using parallel
    n_sigma = 2
    n_flg = 2
    n_excitation = 2
    p_common = (n_site, T_list, exct, eta, path_to_HPhi, header, output_dir, exct)
    one_body_green = calc_one_body_green_core_parallel(p_common)
    np.save("test_g", one_body_green)

if __name__ == "__main__":
    test_main()
