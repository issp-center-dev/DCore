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

        lines = f.readlines()

        # Cound how many R points are defined
        rmap = {} # (int,int,int) -> index of R point
        rvec = []
        ir = 0
        for line in lines:
            i1, i2, i3, _, _, _, _ = line.split()
            i1, i2, i3 = map(int, (i1, i2, i3))
            r = (i1, i2, i3)
            if r in rmap:
                continue
            rvec.append(r)
            rmap[r] = ir
            ir += 1
        assert ir == self.nrpts

        self.irvec = numpy.array(rvec, dtype=numpy.int64)
        self.HamR = numpy.zeros((self.nrpts, self.Nwann, self.Nwann), dtype=complex)

        for line in lines:
            i1, i2, i3, i, j, hr_real, hr_imag = line.split()
            i1, i2, i3, i, j = map(int, (i1, i2, i3, i, j))
            ir = rmap[(i1,i2,i3)]
            self.HamR[ir, i-1, j-1] = complex(float(hr_real), float(hr_imag))


    def get_Hk(self, kvec):
        """
        Compute H(k)
        :param kvec: (float, float, float). Fraction coordinates in k space.
        :return: matrix of H(k)
        """

        Hk = numpy.zeros((self.Nwann, self.Nwann),dtype=numpy.complex128)
        for iR in range(self.nrpts):
            factor = numpy.exp(2J*numpy.pi*(self.irvec[iR,0]*kvec[0]+self.irvec[iR,1]*kvec[1]+self.irvec[iR,2]*kvec[2]))
            Hk += self.HamR[iR,:,:] * factor / self.ndgen[iR]
        return Hk
