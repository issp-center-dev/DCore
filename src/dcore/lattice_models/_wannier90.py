import numpy
from itertools import product

class Wannier90(object):
    def __init__(self, seedname, verbose=0):
        file_name = seedname + "_hr.dat"
        f = open(file_name)
        f.readline()
        self.Nwann = int(f.readline())
        self.nrpts = int(f.readline())

        if verbose > 0:
            print("Num of Wannier functions = ", self.Nwann)
            print("Num of R points = ", self.nrpts)

        num_lines = self.nrpts // 15
        if self.nrpts%15 != 0:
            num_lines += 1
        ndgen = []
        for iline in range(num_lines):
            ndgen.extend(f.readline().split())
        self.ndgen = numpy.array(ndgen, dtype=int)

        self.HamR = numpy.zeros((self.nrpts, self.Nwann, self.Nwann), dtype=complex)
        self.irvec = numpy.zeros((self.nrpts, 3), dtype=int)
        for ir in range(self.nrpts):
            for j in range(self.Nwann):
                for i in range(self.Nwann):
                    i1, i2, i3, i_in, j_in, hr_real,hr_imag = f.readline().split()
                    if i==0 and j==0:
                        self.irvec[ir,0] = i1
                        self.irvec[ir,1] = i2
                        self.irvec[ir,2] = i3
                    assert i == int(i_in)-1
                    assert j == int(j_in)-1
                    self.HamR[ir, i, j] = complex(float(hr_real), float(hr_imag))

    def get_Hk(self, kvec):
        """
        Compute H(k)
        :param kvec: (float, float, float). Fraction coordinates in k space.
        :return: matrix of H(k)
        """

        Hk = numpy.zeros((2*self.norb, 2*self.norb),dtype=complex)
        for iR in range(self.nrpts):
            factor = numpy.exp(2J*numpy.pi*(self.irvec[iR,0]*kvec[0]+self.irvec[iR,1]*kvec[1]+self.irvec[iR,2]*kvec[2]))
            Hk += self.HamR_full[iR,:,:] * factor / self.ndgen[iR]
        return Hk
