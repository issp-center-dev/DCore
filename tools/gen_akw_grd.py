

# a b c alpha beta gamma
# Na Nb Nc
# (((z) y) x)

import numpy as np
import math
import linecache
import argparse
import scipy.interpolate as interp


def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('path_input_file',
                        action='store',
                        default=None,
                        type=str,
                        help="input file name.")
    
    parser.add_argument('path_output_file',
                        action='store',
                        default=None,
                        type=str,
                        help="output file name.")
    parser.add_argument('--omega', help='real_frequency')
    parser.add_argument('--kz', help='wave vector along the z-direction')
    args = parser.parse_args()
    if args.omega:
        fomega = float(args.omega)
    if args.kz:
        Nkz_new = int(args.kz)
    
    print("Opening {}...".format(args.path_input_file))
    with open(args.path_input_file, 'r') as fin:
        l = fin.readline()
        N_kx, N_ky, N_kz, N_omega = list(map(int, l.replace('#', '').split()))
        N_k = N_kx * N_ky * N_kz
    
        print("N_kx= ", N_kx)
        print("N_ky= ", N_ky)
        print("N_kz= ", N_kz)
        print("N_omega= ", N_omega)
    
        bvec = []
        for i in range(3):
            bvec.append(np.array(list(map(float, fin.readline().replace('#', '').split()))))
            print("bvec{}: {}".format(i, bvec[-1]))
    
    fout = open(args.path_output_file, 'w')
    
    # Find a, b, c, alpha, beta, gamma from vec
    # get norm of each vector -> a b c
    a = np.zeros(3)
    angle = np.zeros(3)  # order: alpha beta gamma
    print("#", file=fout)
    for i in range(3):
        a[i] = np.linalg.norm(bvec[i])
        print(a[i], end=" ", file=fout)
    
    # Get angle between each vector -> alpha beta gamma
    for i in range(3):
        temp = math.acos(np.dot(bvec[(i+1) % 3], bvec[(i+2) % 3])/(a[(i+1) % 3]*a[(i+2) % 3]))
        angle[i] = math.degrees(temp)
        print(angle[i], end=" ", file=fout)
    print(" ", file=fout)
    
    
    # read mesh data
    data = np.loadtxt(args.path_input_file).reshape((N_kx, N_ky, N_kz, N_omega, 5))
    
    if args.omega:
        print(N_kx, N_ky, N_kz, file=fout)
        omega_mesh = data[0, 0, 0, :, 3]
        omega_min = omega_mesh[0]
        omega_max = omega_mesh[-1]
        print("omega_min={}, omega_max={}".format(omega_min, omega_max))
        if fomega < omega_min or omega_max < fomega:
            raise RuntimeError("omega={} is out of interpolation range!".format(fomega))
    
        akw = data[:, :, :, :, 4].reshape((N_k, N_omega))
    
        # interpolate
        """
        for k in range(N_k):
            x = np.linspace(omega_min, omega_max, N_omega)
            yup = akw[0][k][:]
            ydw = akw[1][k][:]
            xx = np.linspace(omega_min, omega_max, N_omega_new)
            yyup = interp.spline(x, yup, xx)
            yydw = interp.spline(x, ydw, xx)
            for w in range(N_omega_new):
               if abs(xx[w] - fomega) < ((omega_max - omega_min)/N_omega_new/2):
                    akwout[0][k] = yyup[w]
                    akwout[1][k] = yydw[w]
                    flag = 1
        """
    
        print("Interpolating...")
        for k in range(N_k):
            f = interp.interp1d(omega_mesh, akw[k, :], kind='cubic')
            print(f(fomega), file=fout)
        print("done")
    
    if args.kz:
        print(N_kx, N_ky, N_kz * Nkz_new, file=fout)
        kz_mesh = data[0, 0, :, 0, 2]
        kz_mesh_new = np.linspace(kz_mesh[0], kz_mesh[-1], N_kz * Nkz_new)
        print(kz_mesh.size, kz_mesh_new.size)
        Nxy = N_kx * N_ky
        if N_omega % 2 != 1:
            raise RuntimeError("Nomega should be an odd number", format(N_omega))
        akw = data[:, :, :, N_omega/2, 4].reshape((N_kx * N_ky, N_kz))
    
        # interpolate
        print("Interpolating...")
        for kxy in range(Nxy):
            fxy = interp.interp1d(kz_mesh, akw[kxy, :], kind='cubic')
            for kz in range(kz_mesh_new.size):
                print(fxy(kz_mesh_new[kz]), file=fout)
        print("done")