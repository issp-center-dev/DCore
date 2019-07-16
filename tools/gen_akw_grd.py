from __future__ import print_function
# a b c alpha beta gamma
# Na Nb Nc
# (((z) y) x)
import numpy as np
import math
import linecache
import argparse
import scipy.interpolate as interp
parser = argparse.ArgumentParser()
parser.add_argument('--omega', help='real_frequency')
args = parser.parse_args()
# print('omega='+args.omega)
fomega = float(args.omega)

N_omega = 100
N_kx = 10
N_ky = 10
N_kz = 10
N_k = N_kx * N_ky * N_kz
Flavor = 2
N_omega_new = N_omega * 10  # #datas in the interpolated data
fup = open("akw_up.grd", "w")
fdw = open("akw_dw.grd", "w")

vec = np.array([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
# find a, b, c, alpha, beta, gamma from vec
# get norm of each vector -> a b c
a = np.zeros(3)
angle = np.zeros(3)  # order: alpha beta gamma
print("#", file=fup)
print("#", file=fdw)
for i in range(3):
    a[i] = np.linalg.norm(vec[i])
    print(a[i], end=" ", file=fup)
    print(a[i], end=" ", file=fdw)

# get angle between each vector -> alpha beta gamma
for i in range(3):
    temp = math.acos(np.dot(vec[(i+1) % 3], vec[(i+2) % 3])/a[(i+1) % 3]/a[(i+2) % 3])
    angle[i] = math.degrees(temp)
    print (angle[i], end=" ", file=fup)
    print (angle[i], end=" ", file=fdw)
print(" ", file=fup)
print(" ", file=fdw)

print(N_kx, N_ky, N_kz, file=fup)
print(N_kx, N_ky, N_kz, file=fdw)


# read from case_akw_mesh.dat
infile = open("square_akw_mesh.dat", 'r')
lines = infile.readlines()
infile.close()

line_num = 0
akw = np.zeros((Flavor, N_k, N_omega))
akwout = np.zeros((Flavor, N_k))
# omega_list = np.zeros(N_omega)
omega_min = 0
omega_max = 0
flag = 0
for line in lines:
    if line_num > 0:
        column = line.strip().split()
        if column[0] == 'up':
            spin = 0
        else:
            spin = 1
        if line_num == 1:
            omega_min = float(column[4])
            print("omega_min:", column[4])
        elif line_num == N_omega:
            omega_max = float(column[4])
            print("omega_max:", column[4])
        akw[spin][(line_num-1) / N_omega - spin * N_k][(line_num-1) % N_omega] = column[5]
        # omega_list = column[4]

    line_num += 1


# for spin in range(2):
#    for k in range(N_k):
#        print(akwout[spin][k],akw[spin][k][omega_index])

# interpolate
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

if flag == 0:
    print("something wrong")

for spin in range(2):
    for k in range(N_k):
        print(akwout[0][k], file=fup)
        print(akwout[1][k], file=fdw)
    print(" ")
