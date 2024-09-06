import numpy as np
import scipy.sparse as sp
import scipy.linalg
import argparse
import json
import sys
import time
from itertools import product
from collections import namedtuple


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
        print(f"{self.prefix}{minutes}m{seconds:.3f}s")


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


def check_orthonormality(eigvecs, ignore_orthonormality=False):
    # Check orthonormality of eigenvectors
    overlap = eigvecs.conj().T @ eigvecs  # <m|n>
    if not np.allclose(overlap, np.identity(eigvecs.shape[1])):
        sys.stdout.flush()
        print("Warning: Eigenvectors are not orthonormal!", file=sys.stderr)
        if ignore_orthonormality:
            print("  -> Continue. Set ignore_orthonormality=False to abort calculation.", file=sys.stderr)
        else:
            print("  -> Exit. Set ignore_orthonormality=True to ignore this warning.", file=sys.stderr)
            sys.exit(1)


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
                if i // 2 != j // 2:  # skip different spins
                    continue

            gf[i, j, l] = np.einsum("m, m, m", cdag_im[i].conj().T, ene_denom, cdag_im[j])

    return gf


# Compute
#   <n| c_i (iw - H + E_n) c_j^+ |n>
# by solving linear equations
def calc_gf_iterative(iws, Cdag, spin_conserve, eigvec, E_n, hamil_ex, pm):
    n_flavors = Cdag.size
    n_iw = iws.size
    dim_ex = hamil_ex.shape[0]
    assert hamil_ex.shape == (dim_ex, dim_ex)

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
            # x = sp.linalg.spsolve(A, cdag_j)
            # x = np.linalg.solve(A.toarray(), cdag_j)

            # Use BiCG algorithm
            # x, _ = sp.linalg.bicg(A, cdag_j, x0=x0)
            x, _ = sp.linalg.bicgstab(A, cdag_j, x0=x0)

            # Use LGMRES algorithm
            # x, _ = sp.linalg.lgmres(A, cdag_j, x0=x0)
            # x, _ = sp.linalg.lgmres(A, cdag_j, x0=x0, outer_v=outer_v)

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
    print(params)

    weight_threshold = 1e-6
    ignore_orthonormality = True
    particle_numbers = None

    n_flavors = params['n_flavors']
    n_sites = params['n_sites']
    beta = params['beta']
    n_eigen = params['n_eigen']
    n_iw = params['n_iw']
    flag_spin_conserve = params['flag_spin_conserve']
    dim_full_diag = params['dim_full_diag']
    ncv = params['ncv']

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

    print("\nCreation/annihilation operators info:", flush=True)
    print_sparse_matrix_info(Cdag[0], prefix=" | ")

    # ----------------------------------------------------------------
    # Make many-body Hamiltonian matrix

    # TODO: define density matrix first

    print("\nMaking many-body Hamiltonian matrix...", flush=True)
    hamil = 0
    for i, j in np.ndindex(h0.shape):
        hamil += h0[i, j] * Cdag[i] @ C[j]
    # print(hamil.shape)
    assert hamil.shape == (dim, dim)

    # U_ijkl
    for i, j, k, l in np.ndindex(umat.shape):
        hamil += 0.5 * umat[i, j, k, l] * Cdag[i] @ Cdag[j] @ C[l] @ C[k]

    assert isinstance(hamil, sp.spmatrix)
    assert hamil.shape == (dim, dim)

    print("\nHamiltonian matrix info:", flush=True)
    print_sparse_matrix_info(hamil, prefix=" | ")

    # ----------------------------------------------------------------
    # Particle number conservation

    # particle numbers to be considered
    if particle_numbers is None:
        particle_numbers = np.arange(2*n_sites+1)

    # 0 <= n <= 2*n_sites
    assert np.all((particle_numbers >= 0) & (particle_numbers <= 2*n_sites))

    print("\nParticle number conservation:", flush=True)
    print("  N dim[N]")
    dims = np.zeros(2*n_sites+1, dtype=int)
    indices = np.empty(2*n_sites+1, dtype=object)
    for N in particle_numbers:
        indices[N] = [i for i, state in enumerate(product([0, 1], repeat=2*n_sites)) if sum(state)==N]
        dims[N] = len(indices[N])
        print(f" {N:2d} {dims[N]}")
    assert np.sum(dims) == dim

    # Split the Hamiltonian, Cdag, and C matrices into blocks
    hamil = sp.csr_matrix(hamil)  # for slicing
    hamils = np.empty(2*n_sites+1, dtype=object)
    for N in particle_numbers:
        hamils[N] = slice_spmatrix(hamil, indices, N, N)
    del hamil

    # retain only Cdag and C at impurity sites
    Cdags = np.empty((2*n_sites+1, n_flavors), dtype=object)
    Cs = np.empty((2*n_sites+1, n_flavors), dtype=object)
    for i in range(n_flavors):
        cdag_i = sp.csr_matrix(Cdag[i])
        c_i = sp.csr_matrix(C[i])
        for N in particle_numbers:
            Cdags[N, i] = slice_spmatrix(cdag_i, indices, N+1, N)
            Cs[N, i] = slice_spmatrix(c_i, indices, N-1, N)
        del cdag_i, c_i
    del Cdag, C

    # ----------------------------------------------------------------
    # Solve the eigenvalue problem

    print("\nSolving the eigenvalue problem...", flush=True)
    timer = Timer(prefix="Time: ")
    eigvals = np.empty(2*n_sites+1, dtype=object)
    eigvecs = np.empty(2*n_sites+1, dtype=object)
    full_diagonalization = np.full(2*n_sites+1, False)
    for N in particle_numbers:
        print(f"\nN = {N}  (dim[N] = {dims[N]})")
        full_diagonalization[N] = (dims[N] <= dim_full_diag) or (n_eigen >= dims[N] - 1)

        _timer = Timer(prefix=" Time: ")
        if full_diagonalization[N]:
            print(" full diagonalization")
            # n_eigen = dim
            eigvals[N], eigvecs[N] = scipy.linalg.eigh(hamils[N].toarray())
        else:
            print(" Lanczos method")
            eigvals[N], eigvecs[N] = sp.linalg.eigsh(hamils[N], k=n_eigen, which='SA', ncv=ncv)
            # ‘SA’ : Smallest (algebraic) eigenvalues.
        _timer.print()

        # TODO: orthogonalize eigenvectors
        check_orthonormality(eigvecs[N], ignore_orthonormality=ignore_orthonormality)

    print("\nFinish the eigenvalue problem", flush=True)
    timer.print()

    # One-dimensionalize
    eigs = []
    EIG = namedtuple('EIG', ['val', 'N', 'i'])
    for N in particle_numbers:
        eigs += [EIG(eigval, N, i) for i, eigval in enumerate(eigvals[N])]
    print("\nTotal eigenvalues computed: ", len(eigs))

    # Sort eigenvalues and eigenvectors
    #   The eigenvalues are not ordered in scipy.sparse.linalg.eigsh()
    #   for complex Hermitian matrices since eigs() is called internally.
    eigs = sort_eigs(eigs)

    print("\nEigenvalues:", flush=True)
    E = np.array([eig.val for eig in eigs])
    print(E)

    # Boltzmann factors
    weights = np.exp(-beta * (E - E[0]))
    weights /= np.sum(weights)
    print("\nWeights (Boltzmann factors / Z):", flush=True)
    print(weights)

    n_initial_states = np.count_nonzero(weights > weight_threshold)
    print("\nNumber of initial states:", n_initial_states, flush=True)

    # Check if n_eigen is enough
    # TODO: improve this condition
    if n_eigen < len(weights) and weights[n_eigen] > weight_threshold:
        sys.stdout.flush()
        print(f"\nWarning: n_eigen={n_eigen} may be too small: The weight for the highest-energy state computed is {weights[-1]}.", file=sys.stderr)

    # Save eigenvalues
    with open("eigenvalues.dat", "w") as f:
        f.write(f"# dim = {dim}\n")
        f.write(f"# n_eigen = {n_eigen} (for each n)\n")
        f.write(f"# N  E_i  Boltzmann_weight\n")
        for i, eig in enumerate(eigs):
            f.write(f"{eig.N}  {eig.val:.8e}  {weights[i]:.5e}\n")

    # ----------------------------------------------------------------
    # Calculate impurity Green's function

    # Matsubara frequencies
    iws = 1j * (2 * np.arange(n_iw) + 1) * np.pi / beta

    print("\nCalculating impurity Green's function...", flush=True)
    timer = Timer(prefix="Time: ")
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
                print(f"\n particle excitation: N + 1 = {N_ex}")
            else:
                N_ex = N - 1
                print(f"\n hole excitation: N - 1 = {N_ex}")

            if 0 <= N_ex <= 2*n_sites:
                # Compute
                #   <n| c_i (iw - H + E_n) c_j^+ |n> for particle excitation
                #   <n| c_j^+ (iw + H - E_n) c_i |n> for hole excitation
                if full_diagonalization[N_ex]:
                    print("  Use the Lehmann representation", flush=True)
                    gf_1 = calc_gf_Lehmann(iws, _Cdag, flag_spin_conserve, eigvec, E_n, eigvals[N_ex], eigvecs[N_ex], pm)

                else:
                    print("  Solve linear equations", flush=True)
                    _timer = Timer(prefix="  Time: ")
                    gf_1 = calc_gf_iterative(iws, _Cdag, flag_spin_conserve, eigvec, E_n, hamils[N_ex], pm)
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

