import subprocess
import itertools
import numpy as np
import os
import sys
from ..tools import launch_mpi_subprocesses


class CalcSpectrum:
    def __init__(self, filename, T_list, exct, eta, mpirun_command="", path_to_HPhi="./HPhi", header="zvo", output_dir="./output"):
        self.filename = filename
        self.T_list = T_list
        self.exct = exct
        self.eta = eta
        self.header = header
        self.output_dir = output_dir
        self.nomega = 0
        self.path_to_HPhi = path_to_HPhi
        self.mpirun_command = mpirun_command

    def Make_Spectrum_Input(self, spectrum_type="single"):
        for idx in range(self.exct):
            with open("calcmod.def") as f:
                lines = f.readlines()
            with open("calcmod_ex.def", "w") as fex:
                for line in lines:
                   words = line.split()
                   if words[0] == "CalcSpec" or words[0] == "OutputExVec" or words[0] == "OutputEigenVec":
                    continue
                   fex.write(line)
                fex.write("CalcSpec    1\n")
            with open("namelist.def") as f:
                lines = f.readlines()
            with open("namelist_ex_{}.def".format(idx), "w") as fex:
                for line in lines:
                    words = line.split()
                    if len(words) == 0:
                        continue
                    if words[0] == "CalcMod" or words[0] == "SpectrumVec":
                        continue
                    if words[0] == "SingleExcitation":
                        continue
                    if words[0] == "PairExcitation":
                        continue
                    fex.write(line)

                fex.write("CalcMod calcmod_ex.def\n")
                fex.write("SpectrumVec    {}_eigenvec_{}\n".format(self.header, idx))
                if spectrum_type == "single":
                    fex.write("SingleExcitation single_ex.def\n")
                elif spectrum_type == "pair":
                    fex.write("PairExcitation pair_ex.def\n")

        with open("modpara.def", "r") as fr:
            lines = fr.readlines()
            for line in lines:
                words = line.split()
                if words[0] == "NOmega":
                    self.nomega = int(words[1])
        if self.nomega == 0:
            print("Error: Please set NOmega in modpara file")
            sys.exit(1)

    def read_spectrum(self):
        spectrum_dict={}
        frequencies =[]
        for idx in range(self.exct):
            path_to_DG = os.path.join(self.output_dir, "{}_DynamicalGreen_{}.dat".format(self.header,idx))
            spectrum = np.loadtxt(path_to_DG)
            spectrum_dict[idx] = spectrum[:,2] + 1J*spectrum[:,3]
            if idx == 0 :
                frequencies = spectrum[:, 0] + 1J*spectrum[:, 1]
        self.spectrums_dict = spectrum_dict
        self.frequencies = frequencies
        return spectrum_dict

    def get_energies(self):
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
        for T in self.T_list:
            eta_ene = np.exp(-(self.ene_max-self.ene_min)/T)
            print("T = {}: exp[-beta(ene_max-ene_mix)] = {}".format(T, eta_ene))
            if eta_ene > self.eta:
                print("Warning: At T = {}, eta_ene is larger than eta.".format(T))
        return energy_list

    def _calc_Z(self, T):
        Z = 0
        for ene in self.energy_list:
            ene_diff = ene-self.ene_min
            Z += np.exp(-ene_diff/T)
        return Z

    def get_finite_T_spectrum(self):
        self.read_spectrum()
        finite_T_spectrum_dict ={}
        for T in self.T_list:
            Z = self._calc_Z(T)
            spectrum = np.zeros_like(self.spectrums_dict[0])
            for idx in range(self.exct):
                spectrum += np.exp(-(self.energy_list[idx]-self.ene_min)/T)*self.spectrums_dict[idx]
            spectrum /= Z
            finite_T_spectrum_dict[T]=spectrum
        self.finite_T_spectrum_dict = finite_T_spectrum_dict
        return self.frequencies, finite_T_spectrum_dict

    def print_finite_T_spectrum(self, file_name = "Dynamical_Green"):
        for key, spectrum in self.finite_T_spectrum_dict.items():
            file_name_T = self.header + "_" + file_name + "_T_{}.dat".format(key)
            with open(os.path.join(self.output_dir, file_name_T), "w") as fw:
                for idx, value in enumerate(spectrum):
                    fw.write("{} {} {} {}\n".format(self.frequencies[idx].real, self.frequencies[idx].imag, value.real, value.imag))

    def _make_single_excitation(self, site_i, sigma_i, site_j, sigma_j, file_name = "single_ex.def", ex_state=0, flg_complex = True):
        # c_{i sigma_i} or c_{i sigma_i} + i c_{j sigma_j}
        nsingle = 2 
        if (2 * site_i + sigma_i) == ( 2 * site_j + sigma_j):
            nsingle = 1
        with open(file_name, "w") as fw:
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

    def _run_HPhi(self, exct_cut):
        for idx in range(exct_cut):
            #print("Process: {}/{}".format(idx, exct_cut))
            input_path = "namelist_ex_{}.def".format(idx)
            exec_path = self.path_to_HPhi
            with open('./stdout_{}.log'.format(idx), 'w') as output_f:
                launch_mpi_subprocesses(self.mpirun_command, [exec_path, '-e', input_path], output_f)
            cmd = "mv ./output/{0}_DynamicalGreen.dat ./output/{0}_DynamicalGreen_{1}.dat".format(self.header, idx)
            subprocess.call(cmd, shell=True)

    def get_one_body_green(self, n_site, exct_cut):
        self.Make_Spectrum_Input()
        one_body_green = {}
        for T in self.T_list:
            one_body_green[T] = np.zeros((n_site, 2, n_site, 2, self.nomega), dtype=np.complex)
        #diagonal 
        print("Calculate Diagonal Green function")
        for sitei, sigmai in itertools.product(range(n_site), range(2)):
            print("G[{},{}][{},{}]".format(sitei, "u" if sigmai == 0 else "d", sitei, "u" if sigmai == 0 else "d"))
            for ex_state in range(2):
                self._make_single_excitation(sitei, sigmai, sitei, sigmai, ex_state=ex_state)
                #Run HPhi
                self._run_HPhi(exct_cut)
                #Get Finite-T Green
                frequencies, finite_spectrum_list = self.get_finite_T_spectrum()
                for T in self.T_list:
                    one_body_green[T][sitei][sigmai][sitei][sigmai] += finite_spectrum_list[T]
                
        #off diagonal
        print("Calculate Off-Diagonal Green function")

        one_body_green_tmp = [{}, {}]
        for T in self.T_list:
            one_body_green_tmp[0][T] = np.zeros((self.nomega), dtype=np.complex)
            one_body_green_tmp[1][T] = np.zeros((self.nomega), dtype=np.complex)

        for sitei, sigmai in itertools.product(range(n_site), range(2)):
            sitei_idx = 2 * sitei +  sigmai   
            for sitej, sigmaj in itertools.product(range(n_site), range(2)): 
                sitej_idx = 2 * sitej + sigmaj
                if sitei_idx >= sitej_idx:
                    continue
                print("G[{},{}][{},{}]".format(sitei, "u" if sigmai == 0 else "d", sitej, "u" if sigmaj == 0 else "d"))
                for ex_state in range(2):
                    #True c_i + c_j, False: c_i + i c_j
                    for idx, flg in enumerate([True, False]):
                        self._make_single_excitation(sitei, sigmai, sitej, sigmaj, ex_state=ex_state, flg_complex=flg)
                        #Run HPhi
                        self._run_HPhi(exct_cut)
                        #Get Finite-T Green
                        frequencies, finite_spectrum_list = self.get_finite_T_spectrum()
                        for T in self.T_list:
                            one_body_green_tmp[idx][T] += finite_spectrum_list[T]
                for T in self.T_list:
                    one_body_green_tmp[idx][T] -= one_body_green[T][sitei][sigmai][sitei][sigmai]+one_body_green[T][sitej][sigmaj][sitej][sigmaj]
                #Get Offdiagonal Green
                for T in self.T_list:
                    one_body_green[T][sitei][sigmai][sitej][sigmaj] = (one_body_green_tmp[0][T] + 1J * one_body_green_tmp[1][T] )/2.0
                    one_body_green[T][sitej][sigmaj][sitei][sigmai] = (one_body_green_tmp[0][T] - 1J * one_body_green_tmp[1][T] )/2.0
        return one_body_green
