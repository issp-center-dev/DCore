import numpy as np
import scipy.sparse as sp
import scipy.linalg
import argparse
import json
import sys
import time
from itertools import product
from collections import namedtuple
from pprint import pp
from dcore.impurity_solvers.lanczos import LanczosEigenSolver

# Sparse solvers for eigenvalue problems
eigsolver = (
    # scipy.sparse.linalg.eigsh()
    'eigsh',

    # Lanczos method (self-implemented)
    'lanczos',
)


# Sparse solvers for linear equations
# Choose by 'gf_solver' parameter
spsolver = {
    # Direct methods
    # No parameters
    'spsolve': sp.linalg.spsolve,

    # Iterative methods
    # Parameters:
    #   gf_atol, gf_rtol: absolute and relative tolerance
    'bicg': sp.linalg.bicg,
    'bicgstab': sp.linalg.bicgstab,
    'cg': sp.linalg.cg,
    'cgs': sp.linalg.cgs,
    'gmres': sp.linalg.gmres,
    'lgmres': sp.linalg.lgmres,
    'minres': sp.linalg.minres,
    'qmr': sp.linalg.qmr,
    'gcrotmk': sp.linalg.gcrotmk,
    'tfqmr': sp.linalg.tfqmr,

    # Lanczos method with continued fraction (self-implemented)
    'lanczos': None,
}


class Timer:
    def __init__(self, prefix="Elapsed time: "):
        """Start the timer."""
        self.prefix = prefix
        self.start_time = time.time()
        self.end_time = None

    def restart(self):
        """Restart the timer."""
        self.start_time = time.time()

    def print(self):
        """Print the elapsed time in 'XmXs' format."""
        if self.start_time is None:
            raise RuntimeError("Timer has not been started.")

        self.end_time = time.time()
        elapsed_time = self.end_time - self.start_time

        # Convert elapsed time to minutes and seconds
        minutes = int(elapsed_time // 60)
        seconds = elapsed_time % 60

        # Print in the desired format
        print(f"{self.prefix}{minutes}m{seconds:.3f}s", flush=True)


def make_local_ops():
    ops = {}
    ops['c^+'] = sp.coo_matrix([[0, 0], [1, 0]], dtype=int)
    ops['c'] = sp.coo_matrix([[0, 1], [0, 0]], dtype=int)
    ops['I'] = sp.coo_matrix(np.identity(2, dtype=int))
    ops['F'] = sp.coo_matrix([[1, 0], [0, -1]], dtype=int)  # to represent fermionic anticommutation
    return ops


# kronecker product for an arbitrary number of operators
def kron(*ops):
    r = 1
    for op in ops:
        r = sp.kron(r, op, format='coo')
    return r


def slice_spmatrix(matrix, indices, N1, N2):
    if 0 <= N1 < indices.size and 0 <= N2 < indices.size:
        return sp.csr_matrix(matrix[np.ix_(indices[N1], indices[N2])])
    else:
        return None


def sort_eigs(eigs):
    # Sort eigenvalues and eigenvectors in ascending order
    eigs.sort(key=lambda x: x.val)
    return eigs


def convert_to_real_dtype(array, name="array"):
    # Check if array are real, and if so, convert to real dtype.
    if np.any(np.iscomplexobj(array)):  # check dtype
        if np.any(np.iscomplex(array)):  # check values
            print(f"{name} is complex. Keep dtype={array.dtype}.")
            return array
        else:
            print(f"{name} is real. -> Convert to real dtype.")
            return array.real
    else:
        print(f"{name} is real dtype.")
        return array


def sort_eigenvalues(eigvals, eigvecs):
    # Sort eigenvalues and eigenvectors in ascending order
    idx = np.argsort(eigvals)
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:, idx]
    return eigvals, eigvecs


def orthogonalize(v1, v2):
    # Orthogonalize v1 with respect to v2, assuming v2 is normalized
    return v1 - v2 * np.vdot(v2, v1)
    # return v1 - v2 * (np.vdot(v2, v1) / np.vdot(v2, v2))


def normalize(v):
    # Normalize v
    return v / np.linalg.norm(v)


# Orthonormalize degenerated eigenvectors
def orthonormalize_eigvecs(eigvals, eigvecs, tol=1e-6):
    dim = eigvals.size

    # sort first
    eigvals, eigvecs = sort_eigenvalues(eigvals, eigvecs)

    i_start = 0
    eigvecs[:, 0] = normalize(eigvecs[:, 0])

    for i in range(1, dim):
        if abs(eigvals[i] - eigvals[i-1]) < tol:  # degenerate
            # orthogonalize
            for j in range(i_start, i):
                eigvecs[:, i] = orthogonalize(eigvecs[:, i], eigvecs[:, j])
        else:
            i_start = i
        eigvecs[:, i] = normalize(eigvecs[:, i])

    return eigvals, eigvecs


def func_check_orthonormality(eigvecs, check_orthonormality=True):
    # Check orthonormality of eigenvectors
    overlap = eigvecs.conj().T @ eigvecs  # <m|n>
    if not np.allclose(overlap, np.identity(eigvecs.shape[1])):
        sys.stdout.flush()
        print("Warning: Eigenvectors are not orthonormal!", file=sys.stderr)

        if check_orthonormality:
            print("  -> Exit. Set check_orthonormality{bool}=False to continue calculation.", file=sys.stderr)
            sys.exit(1)
        else:
            print("  -> Continue. Set check_orthonormality{bool}=True to abort calculation.", file=sys.stderr)


def print_ndarray_info(array, prefix=""):
    assert isinstance(array, np.ndarray)
    print(f"{prefix}type = {type(array)}")
    print(f"{prefix}shape = {array.shape}")
    print(f"{prefix}dtype = {array.dtype}")


def print_sparse_matrix_info(matrix, prefix=""):
    assert isinstance(matrix, sp.spmatrix)
    print(f"{prefix}type = {type(matrix)}")
    print(f"{prefix}shape = {matrix.shape}")
    print(f"{prefix}dtype = {matrix.dtype}")
    print(f"{prefix}number of non-zero elements = {matrix.nnz}")
    density = matrix.nnz / np.prod(matrix.shape)
    print(f"{prefix}rate of non-zero elements = {density:.6f}")


# Compute
#   <n| c_i (iw - H + E_n) c_j^+ |n>
# using the eigenvalues of H (Lehmann representation)
def calc_gf_Lehmann(iws, Cdag, spin_conserve, eigvec, E_n, eigvals_ex, eigvecs_ex, pm):
    n_flavors = Cdag.size
    n_iw = iws.size
    dim_ex = eigvals_ex.size
    assert eigvecs_ex.shape == (dim_ex, dim_ex)

    gf = np.zeros((n_flavors, n_flavors, n_iw), dtype=complex)

    # cdag_im[i, m] = <m|c_i^+|n>  for a given n
    cdag_im = np.empty((n_flavors, dim_ex), dtype=complex)
    for i in range(n_flavors):
        cdag_im[i] = eigvecs_ex.conj().T @ Cdag[i] @ eigvec

    for l, iw in enumerate(iws):
        if pm == +1:
            # ene_denom[m] = 1 / (iw - E_m + E_n)
            ene_denom = 1 / (iw - eigvals_ex + E_n)
        else:
            # ene_denom[m] = 1 / (iw + E_m - E_n)
            ene_denom = 1 / (iw + eigvals_ex - E_n)
        assert ene_denom.shape == (dim_ex,)

        for i, j in np.ndindex(n_flavors, n_flavors):
            if spin_conserve:
                n_orb = n_flavors // 2
                if i // n_orb != j // n_orb:  # skip different spins
                    continue

            gf[i, j, l] = np.einsum("m, m, m", cdag_im[i].conj().T, ene_denom, cdag_im[j])

    return gf


# Compute
#   <n| c_i (iw - H + E_n) c_j^+ |n>
# by solving linear equations
def calc_gf_iterative(iws, Cdag, spin_conserve, eigvec, E_n, hamil_ex, pm, solver):
    n_flavors = Cdag.size
    n_iw = iws.size
    dim_ex = hamil_ex.shape[0]
    assert hamil_ex.shape == (dim_ex, dim_ex)

    spsolve = spsolver[solver]

    gf = np.zeros((n_flavors, n_flavors, n_iw), dtype=complex)

    for j in range(n_flavors):
        # cdag_j = c_j^+ |n>  for a given (j, n)
        cdag_j = Cdag[j] @ eigvec

        x0 = None

        for l, iw in enumerate(iws):
            if pm == +1:
                # A = iw - H + E_n
                A = (iw + E_n) * sp.identity(dim_ex) - hamil_ex
            else:
                # A = iw + H - E_n
                A = (iw - E_n) * sp.identity(dim_ex) + hamil_ex

            assert isinstance(A, sp.spmatrix)
            assert A.shape == (dim_ex, dim_ex)

            # Solve A |x> = c_j^+ |n>
            if solver == 'spsolve':
                x = spsolve(A, cdag_j)
            else:
                x, _ = spsolve(A, cdag_j, x0=x0)
                # x, _ = spsolve(A, cdag_j, x0=x0, rtol=rtol, atol=atol)
                x0 = x  # for the next iteration

            assert x.shape == (dim_ex,)

            for i in range(n_flavors):
                if spin_conserve:
                    if i // 2 != j // 2:  # skip different spins
                        continue

                if i == j:
                    gf[i, j, l] = cdag_j.conj().T @ x
                else:
                    cdag_i = Cdag[i] @ eigvec
                    gf[i, j, l] = cdag_i.conj().T @ x
                    del cdag_i

    return gf


# Compute
#   <n| c_i (iw - H + E_n) c_j^+ |n>
# by self-implemented Lanczos code with continued fraction
def calc_gf_lanczos(iws, Cdag, spin_conserve, eigvec, E_n, hamil_ex, pm, ncv):
    n_flavors = Cdag.size
    n_iw = iws.size
    dim_ex = hamil_ex.shape[0]
    assert hamil_ex.shape == (dim_ex, dim_ex)

    gf = np.zeros((n_flavors, n_flavors, n_iw), dtype=complex)

    # cdag_i = c_i^+ |n>  for a given (j, n)
    cdag_i = np.empty((n_flavors, dim_ex), dtype=complex)
    for i in range(n_flavors):
        # cdag_i = c_i^+ |n>  for a given (j, n)
        cdag_i[i, :] = Cdag[i] @ eigvec

    # <n| c_i (iw - H + E_n) c_i^+ |n> for particle excitation
    # <n| c_i^+ (iw + H - E_n) c_i |n> for hole excitation
    def _calc_gf(cdag_i_n):
        # cdag_i_n = c_i^+ |n>
        lan = LanczosEigenSolver(hamil_ex)
        _gf = lan.calc_gf(pm * iws + E_n, cdag_i_n, ncv=ncv)
        del lan
        return pm * _gf

        # lan = LanczosEigenSolver(pm * hamil_ex)
        # _gf = lan.calc_gf(iws + pm * E_n, cdag_i_n, ncv=ncv)
        # del lan
        # return _gf

    for i in range(n_flavors):
        gf[i, i, :] = _calc_gf(cdag_i[i])

    for i, j in np.ndindex(n_flavors, n_flavors):
        if i <= j:  # do only i > j
            continue

        if spin_conserve:
            if i // 2 != j // 2:  # skip different spins
                continue

        # (c_i^+ + c_j^+) |n>  for particle excitation
        # (c_i + c_j) |n>      for hole excitation
        F_ij_1 = _calc_gf(cdag_i[i] + cdag_i[j])

        # (c_i^+ + i c_j^+) |n>  for particle excitation
        # (c_i + i c_j) |n>      for hole excitation
        F_ij_2 = _calc_gf(cdag_i[i] + 1j * cdag_i[j])
        # F_ij_2 = _calc_gf(cdag_i[i] + 1j * pm * cdag_i[j])

        gf_diag = gf[i, i, :] + gf[j, j, :]
        F_ij_1 -= gf_diag
        F_ij_2 -= gf_diag

        gf[i, j] = (F_ij_1 - 1j * F_ij_2) / 2   # [j, i] for hole
        gf[j, i] = (F_ij_1 + 1j * F_ij_2) / 2   # [i, j] for hole

    return gf


def main():
    # ----------------------------------------------------------------
    # Parse input parameters

    parser = argparse.ArgumentParser(description='Process a string argument.')
    parser.add_argument('filename', type=str, help='Input file in json format')
    args = parser.parse_args()

    # Load input parameters
    with open(args.filename, 'r') as f:
        params = json.load(f)

    assert isinstance(params, dict)
    print("Input parameters:")
    pp(params)

    n_flavors = params['n_flavors']
    n_sites = params['n_sites']
    beta = params['beta']
    n_eigen = params['n_eigen']
    n_iw = params['n_iw']
    flag_spin_conserve = params['flag_spin_conserve']
    dim_full_diag = params['dim_full_diag']
    particle_numbers = params['particle_numbers']
    weight_threshold = params['weight_threshold']
    ncv = params['ncv']
    eigen_solver = params['eigen_solver']
    gf_solver = params['gf_solver']
    # gf_atol = params['gf_atol']
    # gf_rtol = params['gf_rtol']
    check_n_eigen = params['check_n_eigen']
    check_orthonormality = params['check_orthonormality']

    if eigen_solver not in eigsolver:
        raise ValueError(f"Invalid eigen_solver: {eigen_solver}")

    if gf_solver not in spsolver:
        raise ValueError(f"Invalid gf_solver: {gf_solver}")

    # ----------------------------------------------------------------
    # Load H0 and U_ijkl

    # Load H0
    print("\nLoading H0 and U_ijkl", flush=True)
    h0 = np.load(params['file_h0'])
    print("\nH0 info:")
    print_ndarray_info(h0, prefix=" | ")
    # assert h0.shape == (n_flavors, n_flavors)
    assert h0.shape == (2*n_sites, 2*n_sites)

    # Load U_ijkl
    umat = np.load(params['file_umat'])
    print("\nU_ijkl info:", flush=True)
    print_ndarray_info(umat, prefix=" | ")
    assert umat.shape == (n_flavors, n_flavors, n_flavors, n_flavors)

    # Check if H0 and U_ijkl are real. If so, convert to real dtype.
    print("\nConvert H0 and U_ijkl to real dtype if possible", flush=True)
    h0 = convert_to_real_dtype(h0, name="H0")
    umat = convert_to_real_dtype(umat, name="U_ijkl")

    dim = 2 ** (2*n_sites)
    print("\nDimension of Hilbert space:", dim, flush=True)

    # ----------------------------------------------------------------
    # Make creation/annihilation operators

    # building blocks for creation/anihilation operators
    local_ops = make_local_ops()
    cdag = local_ops['c^+']
    I = local_ops['I']
    F = local_ops['F']
    assert cdag.shape == I.shape == F.shape == (2, 2)

    # Make creation operators
    # ex. for i=3:
    #   kron(I, I, I, cdag, F, F)
    # Cdag = []
    print("\nMaking creation/annihilation operators...", flush=True)
    timer = Timer(prefix="Time: ")
    Cdag = np.empty((2*n_sites,), dtype=object)
    C = np.empty((2*n_sites,), dtype=object)
    for i in range(2*n_sites):
        n_I = i
        n_F = 2*n_sites - i - 1
        # print(n_I, n_F)
        args = [I] * n_I + [cdag] + [F] * n_F
        Cdag[i] = kron(*args)
        C[i] = Cdag[i].T  # real (omit conj)
        assert isinstance(Cdag[i], sp.spmatrix)
        assert Cdag[i].shape == (dim, dim)
    timer.print()

    print("\nCreation/annihilation operators info:", flush=True)
    print_sparse_matrix_info(Cdag[0], prefix=" | ")

    # ----------------------------------------------------------------
    # Make many-body Hamiltonian matrix

    # TODO: define density matrix first

    print("\nMaking many-body Hamiltonian matrix...", flush=True)
    timer.restart()
    hamil = 0
    for i, j in np.ndindex(h0.shape):
        if h0[i, j] != 0:
            hamil += h0[i, j] * Cdag[i] @ C[j]
    # print(hamil.shape)
    assert hamil.shape == (dim, dim)

    # U_ijkl
    for i, j, k, l in np.ndindex(umat.shape):
        if umat[i, j, k, l] != 0:
            hamil += 0.5 * umat[i, j, k, l] * Cdag[i] @ Cdag[j] @ C[l] @ C[k]

    assert isinstance(hamil, sp.spmatrix)
    assert hamil.shape == (dim, dim)
    timer.print()

    print("\nHamiltonian matrix info:", flush=True)
    print_sparse_matrix_info(hamil, prefix=" | ")

    # ----------------------------------------------------------------
    # Particle number conservation

    particle_numbers_all = np.arange(2*n_sites+1)

    print("\nParticle-number conservation:", flush=True)
    print("  N       dim")
    print(" -------------")
    dims = np.zeros(2*n_sites+1, dtype=int)
    indices = np.empty(2*n_sites+1, dtype=object)
    for N in particle_numbers_all:
        indices[N] = [i for i, state in enumerate(product([0, 1], repeat=2*n_sites)) if sum(state)==N]
        dims[N] = len(indices[N])
        print(f" {N:2d} {dims[N]:9,d}")
    assert np.sum(dims) == dim

    # Split the Hamiltonian, Cdag, and C matrices into blocks
    hamil = sp.csr_matrix(hamil)  # for slicing
    hamils = np.empty(2*n_sites+1, dtype=object)
    for N in particle_numbers_all:
        hamils[N] = slice_spmatrix(hamil, indices, N, N)
    del hamil

    # retain only Cdag and C at impurity sites
    Cdags = np.empty((2*n_sites+1, n_flavors), dtype=object)
    Cs = np.empty((2*n_sites+1, n_flavors), dtype=object)
    for i in range(n_flavors):
        cdag_i = sp.csr_matrix(Cdag[i])
        c_i = sp.csr_matrix(C[i])
        for N in particle_numbers_all:
            Cdags[N, i] = slice_spmatrix(cdag_i, indices, N+1, N)
            Cs[N, i] = slice_spmatrix(c_i, indices, N-1, N)
        del cdag_i, c_i
    del Cdag, C

    # particle numbers to be considered
    flag_prescreening = False
    if particle_numbers == 'all':
        particle_numbers = particle_numbers_all
    elif particle_numbers == 'auto':
        # particle_numbers is determined by the prescreening
        flag_prescreening = True
    elif isinstance(particle_numbers, (list, tuple)):
        particle_numbers = np.array(particle_numbers, dtype=int)
    else:
        raise ValueError(f"Invalid particle_numbers: {particle_numbers}")

    # ----------------------------------------------------------------
    # Prescreening

    if flag_prescreening:
        print("\nPrescreening -- Solving the eigenvalue problem...", flush=True)
        timer.restart()
        eigvals = np.zeros(2*n_sites+1, dtype=float)
        for N in particle_numbers_all:
            print(f"\nN = {N}  (dim = {dims[N]})", flush=True)

            _timer = Timer(prefix=" Time: ")
            if dims[N] <= 2:
                # print(" full diagonalization", flush=True)
                # n_eigen = dim[N]
                _eigvals = scipy.linalg.eigh(hamils[N].toarray(), eigvals_only=True)
                eigvals[N] = _eigvals[0]
            else:
                # print(f" Iterative solver: The smallest eigenvalue is computed.", flush=True)
                _eigvals = sp.linalg.eigsh(hamils[N], k=1, which='SA', ncv=ncv, return_eigenvectors=False)
                eigvals[N] = _eigvals[0]
                # ‘SA’ : Smallest (algebraic) eigenvalues.

            _timer.print()

        print("\nFinish the eigenvalue problem", flush=True)
        timer.print()

        Emin = np.min(eigvals)
        weights_rel = np.exp(-beta * (eigvals - Emin))  # relative Boltzmann weights

        print("\nSummary of lowest energy state in each N block:", flush=True)
        print("  N       dim    E_min     weight_rel")
        print(" ------------------------------------")
        for N in particle_numbers_all:
            print(f" {N:2d} {dims[N]:>9,d}   {eigvals[N]:9.2e}  {weights_rel[N]:<8.1e}")

        # Select particle numbers of thermally occupied states
        bools_thermal = weights_rel > weight_threshold
        particle_numbers_thermal = np.flatnonzero(bools_thermal)

        print(f"\nParticle numbers of thermally occupied states:\n {particle_numbers_thermal}", flush=True)

        def expand_right_and_left(array):
            expanded = array.copy()
            expanded[1:] |= array[:-1]  # propagate True to the right neighbor
            expanded[:-1] |= array[1:]  # propagate True to the left neighbor
            return expanded

        # N+1 and N-1 are also considered for Green's function calc
        bools_solve = expand_right_and_left(bools_thermal)
        particle_numbers = np.flatnonzero(bools_solve)

    print(f"\nParticle numbers to be considered:\n {particle_numbers}", flush=True)

    # 0 <= n <= 2*n_sites
    assert np.all((particle_numbers >= 0) & (particle_numbers <= 2*n_sites))

    # ----------------------------------------------------------------
    # Solve the eigenvalue problem

    print("\nSolving the eigenvalue problem...", flush=True)
    timer.restart()
    eigvals = np.empty(2*n_sites+1, dtype=object)
    eigvecs = np.empty(2*n_sites+1, dtype=object)
    full_diagonalization = np.full(2*n_sites+1, False)
    for N in particle_numbers:
        print(f"\nN = {N}  (dim = {dims[N]})", flush=True)
        full_diagonalization[N] = (dims[N] <= dim_full_diag) or (n_eigen >= dims[N] - 1)

        _timer = Timer(prefix=" Time: ")
        if full_diagonalization[N]:
            print(" full diagonalization", flush=True)
            # n_eigen = dim[N]
            eigvals[N], eigvecs[N] = scipy.linalg.eigh(hamils[N].toarray())
        else:
            print(f" Iterative solver: n_eigen={n_eigen} eigenvalues are computed.", flush=True)
            if eigen_solver == 'lanczos':
                lanczos = LanczosEigenSolver(hamils[N])
                eigvals[N], eigvecs[N] = lanczos.solve(k=n_eigen, ncv=ncv)
                del lanczos
            else:
                eigvals[N], eigvecs[N] = sp.linalg.eigsh(hamils[N], k=n_eigen, which='SA', ncv=ncv)
                # ‘SA’ : Smallest (algebraic) eigenvalues.

        _timer.print()

        # Orthonormalize eigenvectors
        #   Degenerated eigenstates are not orthogonal with each other
        #   in scipy.sparse.linalg.eigsh() for complex Hermitian matrices
        #   since eigs() is called internally.
        eigvals[N], eigvecs[N] = orthonormalize_eigvecs(eigvals[N], eigvecs[N])

        func_check_orthonormality(eigvecs[N], check_orthonormality=check_orthonormality)

    print("\nFinish the eigenvalue problem", flush=True)
    timer.print()

    # One-dimensionalize
    eigs = []
    EIG = namedtuple('EIG', ['val', 'N', 'i'])
    for N in particle_numbers:
        eigs += [EIG(eigval, N, i) for i, eigval in enumerate(eigvals[N])]
    print("\nTotal eigenvalues computed: ", len(eigs))

    # Sort eigenvalues
    eigs = sort_eigs(eigs)

    E = np.array([eig.val for eig in eigs])

    # Boltzmann factors
    weights = np.exp(-beta * (E - E[0]))
    weights /= np.sum(weights)

    # Save eigenvalues
    print("\nSave eigenvalues and eigenvectors in\n 'eigenvalues.dat'\n 'eigenvectors.dat'", flush=True)
    with open("eigenvalues.dat", "w") as f:
        f.write(f"# dim = {dim}\n")
        f.write(f"# n_eigen = {n_eigen} (for each n)\n")
        f.write(f"# N  E_i  Boltzmann_weight\n")
        for i, eig in enumerate(eigs):
            f.write(f"{eig.N}  {eig.val:.8e}  {weights[i]:.5e}\n")

    # Save eigenvectors
    # print("\nSave eigenvectors in 'eigenvectors.dat'", flush=True)
    with open("eigenvectors.dat", "w") as f:
        length = 2*n_sites
        f.write(f"# The i-th digit from the left of the {length}-dimensional Fock state\n")
        f.write(f"# indicates whether (spin, site) state is occupied (1) or empty (0).\n")
        f.write(f"# 'site' includes bath sites (b) after Wannier orbitals (w).\n")
        f.write(f"#\n")
        f.write(f"#   w1+ w2+ ... b1+ b2+ ... w1- w2- ... b1- b2- ...\n")
        f.write(f"#\n")
        for eig in eigs:
            # print(eig)
            f.write(f"# E={eig.val:.8e}, N={eig.N}\n")
            eigvec = eigvecs[eig.N][:, eig.i]
            for j in range(eigvec.size):
                if abs(eigvec[j]) > 1e-6:  # print only non-zero elements
                    state = indices[eig.N][j]
                    f.write(f" |{state:0{length}b}> {eigvec[j]:.8e}\n")

    n_initial_states = np.count_nonzero(weights > weight_threshold)
    print("\nNumber of initial states:", n_initial_states, flush=True)

    # Check if n_eigen is enough
    # if not np.all(full_diagonalization==True):  # do no check if only full diag is used
    #     # TODO: improve this condition
    #     if n_eigen < len(weights) and weights[n_eigen] > weight_threshold:
    #         sys.stdout.flush()
    #         print(f"\nWarning: n_eigen={n_eigen} may be too small: The weight for the highest-energy state computed is {weights[-1]}.", file=sys.stderr)

    #         if check_n_eigen:
    #             print("  -> Exit. Set check_n_eigen{bool}=False to continue calculation.", file=sys.stderr)
    #             sys.exit(1)
    #         else:
    #             print("  -> Continue. Set check_n_eigen{bool}=True to abort calculation.", file=sys.stderr)

    # number of N-particle states in initial states
    n_states_in_initial = np.zeros(N+1, dtype=int)
    Ns_initial = np.array([eig.N for eig in eigs[0:n_initial_states]])
    for N in particle_numbers:
        n_states_in_initial[N] = np.count_nonzero(Ns_initial == N)
        print(f"  N={N}: {n_states_in_initial[N]} / {dims[N]}")

    # Check if n_eigen is enough
    for N in particle_numbers:
        if full_diagonalization[N]:
            continue
        if n_eigen == n_states_in_initial[N]:
            sys.stdout.flush()
            print(f"\nWarning: n_eigen={n_eigen} is too small for {N}-particle states. Increase n_eigen!", file=sys.stderr)

            if check_n_eigen:
                print("  -> Exit. Set check_n_eigen{bool}=False to continue calculation.", file=sys.stderr)
                sys.exit(1)
            else:
                print("  -> Continue. Set check_n_eigen{bool}=True to abort calculation.", file=sys.stderr)

    # ----------------------------------------------------------------
    # Calculate impurity Green's function

    # Matsubara frequencies
    iws = 1j * (2 * np.arange(n_iw) + 1) * np.pi / beta

    print("\nCalculating impurity Green's function...", flush=True)
    timer.restart()
    gf = np.zeros((n_flavors, n_flavors, n_iw), dtype=complex)
    for n in range(n_initial_states):
        N = eigs[n].N
        E_n = eigs[n].val
        eigvec = eigvecs[N][:, eigs[n].i]

        print(f"\nInitial state {n+1}/{n_initial_states}  (N = {N})", flush=True)

        # loop for [particle excitation, hole excitation]
        for _Cdag, pm in [[Cdags[N], +1], [Cs[N], -1]]:
            if pm == +1:
                N_ex = N + 1
                print(f"\n particle excitation: N + 1 = {N_ex}", flush=True)
            else:
                N_ex = N - 1
                print(f"\n hole excitation: N - 1 = {N_ex}", flush=True)

            if 0 <= N_ex <= 2*n_sites:
                if N_ex not in particle_numbers:
                    print(f"ERROR: N={N_ex} is not in particle_numbers.", file=sys.stderr)
                    sys.exit(1)

                # parameters for gf_solver
                params_gf = dict(
                    iws = iws,
                    Cdag = _Cdag,
                    spin_conserve = flag_spin_conserve,
                    eigvec = eigvec,
                    E_n = E_n,
                    pm = pm,
                )
                # Compute
                #   <n| c_i (iw - H + E_n) c_j^+ |n> for particle excitation
                #   <n| c_j^+ (iw + H - E_n) c_i |n> for hole excitation
                if full_diagonalization[N_ex]:
                    print("  Use the Lehmann representation", flush=True)
                    params_gf.update(
                        eigvals_ex = eigvals[N_ex],
                        eigvecs_ex = eigvecs[N_ex],
                    )
                    gf_1 = calc_gf_Lehmann(**params_gf)

                else:
                    print("  Solve linear equations", flush=True)
                    _timer = Timer(prefix="  Time: ")
                    if gf_solver == 'lanczos':
                        params_gf.update(
                            hamil_ex = hamils[N_ex],
                            ncv = ncv,
                        )
                        gf_1 = calc_gf_lanczos(**params_gf)
                    else:
                        params_gf.update(
                            hamil_ex = hamils[N_ex],
                            solver = gf_solver,
                        )
                        gf_1 = calc_gf_iterative(**params_gf)
                    _timer.print()

                if pm == +1:
                    gf += gf_1 * weights[n]
                else:
                    gf += gf_1.transpose([1, 0, 2]) * weights[n]  # [i, j, l] -> [j, i, l]

    print("\nFinish impurity Green's function", flush=True)
    timer.print()

    # Save Green's function
    np.save("gf", gf)

    # ----------------------------------------------------------------


if __name__ == '__main__':
    main()

