import numpy as np
from itertools import product


t = np.zeros((3, 3, 2, 2))
so = 2  # 1 for no SO, 2 for SO

print(f"square lattice, two orbitals")
print(f"{2*so}")  # number of orbitals
print(f"9")  # number of sites
print(f"1 1 1 1 1 1 1 1 1")

# hopping (orbital diagonal)
t_orb = [-1, -0.5]
for o, tval in enumerate(t_orb):
    t[0, 1, o, o] = tval
    t[1, 0, o, o] = tval
    t[0, -1, o, o] = tval
    t[-1, 0, o, o] = tval

# crystal field splitting
cf = [-0.2, 0.2]
for l, cfval in enumerate(cf):
    t[0, 0, l, l] = cfval

# print in wannier90 format

# spin : outer
for i, j in product(range(-1, 2), repeat=2):
    for s1 in range(so):
        for o1 in range(2):
            for s2 in range(so):
                for o2 in range(2):
                    m = 2*s1 + o1 + 1
                    n = 2*s2 + o2 + 1
                    tval = t[i, j, o1, o2] if s1==s2 else 0
                    print(f"{i:2d} {j:2d} 0  {m} {n}  {tval:4.1f} 0.0")

# spin : inner
# for i, j in product(range(-1, 2), repeat=2):
#     for o1 in range(2):
#         for s1 in range(so):
#             for o2 in range(2):
#                 for s2 in range(so):
#                     m = 2*o1 + s1 + 1
#                     n = 2*o2 + s2 + 1
#                     tval = t[i, j, o1, o2] if s1==s2 else 0
#                     print(f"{i:2d} {j:2d} 0  {m} {n}  {tval:4.1f} 0.0")

