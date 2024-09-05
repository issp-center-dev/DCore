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
    ops['c^+'] = np.array([[0, 0], [1, 0]], dtype=int)
    ops['c'] = np.array([[0, 1], [0, 0]], dtype=int)
    ops['I'] = np.identity(2, dtype=int)
    ops['F'] = np.array([[1, 0], [0, -1]], dtype=int)  # to represent fermionic anticommutation
    return ops


# kronecker product for an arbitrary number of operators
def kron(*ops):
    r = 1
    for op in ops:
        r = sp.kron(r, op)
    return r


def slice_spmatrix(matrix, indices, N1, N2):
    if 0 <= N1 < indices.size and 0 <= N2 < indices.size:
        return sp.lil_matrix(matrix[np.ix_(indices[N1], indices[N2])])
    else:
        return None


def sort_eigenvalues(vals, vecs):
    # Sort eigenvalues and eigenvectors in ascending order
    idx = np.argsort(vals)
    vals = vals[idx]
    vecs = vecs[:, idx]
    return vals, vecs


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


def main():
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

    # Load H0
    print("\nLoading H0 and U_ijkl")
    h0 = np.load(params['file_h0'])
    print("\nH0 info:")
    print_ndarray_info(h0, prefix=" | ")
    # assert h0.shape == (n_flavors, n_flavors)
    assert h0.shape == (2*n_sites, 2*n_sites)

    # Load U_ijkl
    umat = np.load(params['file_umat'])
    print("\nU_ijkl info:")
    print_ndarray_info(umat, prefix=" | ")
    assert umat.shape == (n_flavors, n_flavors, n_flavors, n_flavors)

    # Check if H0 and U_ijkl are real. If so, convert to real dtype.
    print("\nConvert H0 and U_ijkl to real dtype if possible")
    h0 = convert_to_real_dtype(h0, name="H0")
    umat = convert_to_real_dtype(umat, name="U_ijkl")

    dim = 2 ** (2*n_sites)
    print("\nDimension of Hilbert space:", dim)

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
    # print(Cdag)

    # TODO: define density matrix first

    # Make many-body Hamiltonian
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

    print("\nHamiltonian matrix info:")
    print_sparse_matrix_info(hamil, prefix=" | ")

    # Particle number conservation

    # particle numbers to be considered
    if particle_numbers is None:
        particle_numbers = np.arange(2*n_sites+1)

    # 0 <= n <= 2*n_sites
    assert np.all((particle_numbers >= 0) & (particle_numbers <= 2*n_sites))

    print("\nParticle number conservation:")
    print(" N dim[N]")
    dims = np.zeros(2*n_sites+1, dtype=int)
    indices = np.empty(2*n_sites+1, dtype=object)
    # for n in range(dims.size):
    for N in particle_numbers:
        indices[N] = [i for i, state in enumerate(product([0, 1], repeat=2*n_sites)) if sum(state)==N]
        # print(indices)
        dims[N] = len(indices[N])
        print(f" {N} {dims[N]}")
    assert np.sum(dims) == dim

    # Split the Hamiltonian, Cdag, and C matrices into blocks
    hamil = sp.lil_matrix(hamil)  # for slicing
    hamils = np.empty(2*n_sites+1, dtype=object)
    for N in particle_numbers:
        # hamils[N] = sp.lil_matrix(hamil[np.ix_(indices[N], indices[N])])
        hamils[N] = slice_spmatrix(hamil, indices, N, N)
    del hamil

    # retain only Cdag and C at impurity sites
    Cdags = np.empty((n_flavors, 2*n_sites+1), dtype=object)
    Cs = np.empty((n_flavors, 2*n_sites+1), dtype=object)
    for i in range(n_flavors):
        cdag_i = sp.lil_matrix(Cdag[i])
        c_i = sp.lil_matrix(C[i])
        for N in particle_numbers:
            Cdags[i, N] = slice_spmatrix(cdag_i, indices, N+1, N)
            Cs[i, N] = slice_spmatrix(c_i, indices, N-1, N)
        del cdag_i, c_i
    del Cdag, C

    # Use full diagonalization if dim is small
    # full_diagonalization = (n_eigen >= dim - 1)

    # Solve the eigenvalue problem
    print("\nSolving the eigenvalue problem...")
    eigvals = np.empty(2*n_sites+1, dtype=object)
    eigvecs = np.empty(2*n_sites+1, dtype=object)
    full_diagonalization = np.full(2*n_sites+1, False)
    for N in particle_numbers:
        print(f"\nN = {N}  (dim[N] = {dims[N]})")
        # print_sparse_matrix_info(hamils[N], prefix=" | ")
        full_diagonalization[N] = (n_eigen >= dims[N] - 1)

        timer = Timer(prefix=" Time: ")
        # if dims[N] == 1:
        #     eigvals[N] = np.array([hamils[N][0, 0]])
        #     eigvecs[N] = np.array([[1]])
        # elif n_eigen >= dims[N] - 1:
        if full_diagonalization[N]:
            print(" n_eigen >= dim[N]\n  -> Use full diagonalization")
            # n_eigen = dim
            eigvals[N], eigvecs[N] = scipy.linalg.eigh(hamils[N].toarray())
        else:
            print(" n_eigen < dim[N]\n  -> Use Lanczos method")
            eigvals[N], eigvecs[N] = sp.linalg.eigsh(hamils[N], k=n_eigen, which='SA')
            # ‘SA’ : Smallest (algebraic) eigenvalues.
        timer.print()

        check_orthonormality(eigvecs[N], ignore_orthonormality=ignore_orthonormality)

    # assert E.shape == (n_eigen,)
    # assert eigvecs.shape == (dim, n_eigen)

    # One-dimensionalize
    eigs = []
    EIG = namedtuple('EIG', ['val', 'N', 'i'])
    for N in particle_numbers:
        # print(eigvals[N].shape)
        # print(eigvals[N])
        eigs += [EIG(eigval, N, i) for i, eigval in enumerate(eigvals[N])]
    print("\nTotal eigenvalues computed: ", len(eigs))
    # print(eigs)
    # print(eigvals.shape)
    # eigvals = np.concatenate(eigvals)
    # print(eigvals.shape)

    # Sort eigenvalues and eigenvectors
    #   The eigenvalues are not ordered in scipy.sparse.linalg.eigsh()
    #   for complex Hermitian matrices since eigs() is called internally.
    # E, eigvecs = sort_eigenvalues(eigvals, eigvecs)
    eigs = sort_eigs(eigs)

    print("\nEigenvalues:")
    E = np.array([eig.val for eig in eigs])
    print(E)

    # Boltzmann factors
    weights = np.exp(-beta * (E - E[0]))
    weights /= np.sum(weights)
    print("\nWeights (Boltzmann factors / Z):")
    print(weights)

    n_initial_states = np.count_nonzero(weights > weight_threshold)
    print("\nNumber of initial states:", n_initial_states)

    # if weights[-1] > weight_threshold and n_eigen < dim:
    if n_eigen < len(weights) and weights[n_eigen] > weight_threshold:
        sys.stdout.flush()
        print(f"\nWarning: n_eigen={n_eigen} may be too small: The weight for the highest-energy state computed is {weights[-1]}.", file=sys.stderr)

    # Save eigenvalues
    # indexed_E = np.column_stack((np.arange(len(E)), E, weights))
    # header = f"dim = {dim}\nn_eigen = {n_eigen}\ni  E_i  Boltzmann_weight"
    # np.savetxt("eigenvalues.dat", indexed_E, fmt='%d %.8e %.5e', header=header)

    with open("eigenvalues.dat", "w") as f:
        f.write(f"# dim = {dim}\n")
        f.write(f"# n_eigen = {n_eigen} (for each n)\n")
        f.write(f"# N  E_i  Boltzmann_weight\n")
        for i, eig in enumerate(eigs):
            f.write(f"{eig.N}  {eig.val:.8e}  {weights[i]:.5e}\n")


    # Check if eigenvalues are ordered
    # assert np.all(np.diff(E) >= 0)

    # TODO: Save eigenvectors
    # print("Eigenvectors:")
    # print(eigvecs)

    # print("\nWeights:")
    # print(eigvecs * eigvecs.conj())

    # print("\nOrthonormality of eigenvectors:")
    # overlap = eigvecs.conj().T @ eigvecs  # <m|n>
    # print("m n <m|n>")
    # for m, n in np.ndindex(n_eigen, n_eigen):
    #     print(m, n, overlap[m, n])

    # Check orthonormality of eigenvectors
    # print("\nCheck orthonormality of eigenvectors")
    # if not np.allclose(overlap, np.identity(n_eigen)):
    #     print("  not orthonormal")
    #     sys.stdout.flush()
    #     print("Warning: Eigenvectors are not orthonormal!", file=sys.stderr)
    #     if ignore_orthonormality:
    #         print("  -> Continue. Set ignore_orthonormality=False to abort calculation.", file=sys.stderr)
    #     else:
    #         print("  -> Exit. Set ignore_orthonormality=True to ignore this warning.", file=sys.stderr)
    #         sys.exit(1)
    # else:
    #     print("  orthogonal")

    # TODO: orthogonalize eigenvectors

    # Matsubara frequencies
    iws = 1j * (2 * np.arange(n_iw) + 1) * np.pi / beta

    # Calculate impurity Green's function
    print("\nCalculating impurity Green's function...")
    gf = np.zeros((n_flavors, n_flavors, n_iw), dtype=complex)
    for n in range(n_initial_states):
        N = eigs[n].N
        E_n = eigs[n].val
        eigvec = eigvecs[N][:, eigs[n].i]

        print(f"\nInitial state {n+1}/{n_initial_states}  (N = {N})")

        # loop for [particle excitation, hole excitation]
        for _Cdag, particle_excitation in [[Cdags, True], [Cs, False]]:

            if particle_excitation:
                Np = N + 1
                print(f"\n particle excitation: Np = {Np}")
            else:
                Np = N - 1
                print(f"\n hole excitation: Np = {Np}")
            # print(f"Np = {Np}")

            if Np < 0 or Np > 2*n_sites:
                continue

            if full_diagonalization[Np]:
                print("  Use the Lehmann representation")
                # cdag_im[i, m] = <m|c_i^+|n>  for a given n
                cdag_im = np.empty((n_flavors, dims[Np]), dtype=complex)
                for i in range(n_flavors):
                    cdag_im[i] = eigvecs[Np].conj().T @ _Cdag[i, N] @ eigvec

                for l, iw in enumerate(iws):
                    # ene_denom[m] = 1 / (iw - E_m + E_n)
                    # ene_denom_1 = 1 / (iw - E + E[n])
                    if particle_excitation:
                        ene_denom = 1 / (iw - eigvals[Np] + E_n)
                    else:
                        ene_denom = 1 / (iw + eigvals[Np] - E_n)

                    # ene_denom[m] = 1 / (iw + E_m - E_n)
                    # ene_denom_2 = - ene_denom_1.conj()

                    for i, j in np.ndindex(n_flavors, n_flavors):
                        if flag_spin_conserve:
                            if i // 2 != j // 2:  # skip different spins
                                continue

                        # gf[i, j, l] += np.einsum("m, m, m", cdag_im[i].conj().T, ene_denom, cdag_im[i]) * weights[n]
                        gf_1 = np.einsum("m, m, m", cdag_im[i].conj().T, ene_denom, cdag_im[j])
                        # gf_2 = np.einsum("m, m, m", c_im[j].conj().T, ene_denom_2, c_im[i])

                        if particle_excitation:
                            gf[i, j, l] += gf_1 * weights[n]
                        else:
                            gf[j, i, l] += gf_1 * weights[n]  # i <-> j for hole excitation

            else:
                print("  Solve linear equations")
                timer = Timer(prefix="  Time: ")

                # < c_i c_j^+ > for particle excitation
                # < c_i^+ c_j > for hole excitation  (i <-> j later)
                for j in range(n_flavors):
                    # cdag_j = c_j^+ |n>  for a given (j, n)
                    cdag_j = _Cdag[j, N] @ eigvec

                    x0 = None

                    # outer_v = None  # for LGMRES algorithm: Krylov subspace must be prepared

                    for l, iw in enumerate(iws):
                        # A = iw - H + E_n
                        # A = (iw + pm * E[n]) * sp.identity(dim) - pm * hamil
                        if particle_excitation:
                            # A = iw - H + E_n
                            A = (iw + E_n) * sp.identity(dims[Np]) - hamils[Np]
                        else:
                            # A = iw + H - E_n
                            A = (iw - E_n) * sp.identity(dims[Np]) + hamils[Np]
                        assert isinstance(A, sp.spmatrix)
                        assert A.shape == (dims[Np], dims[Np])

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

                        assert x.shape == (dims[Np],)

                        for i in range(n_flavors):
                            if flag_spin_conserve:
                                if i // 2 != j // 2:  # skip different spins
                                    continue

                            if i == j:
                                gf_1 = cdag_j.conj().T @ x
                            else:
                                cdag_i = _Cdag[i, N] @ eigvec
                                gf_1 = cdag_i.conj().T @ x
                                del cdag_i

                            # gf[i, j, l] += gf_1 * weights[n]
                            if particle_excitation:
                                gf[i, j, l] += gf_1 * weights[n]
                            else:
                                gf[j, i, l] += gf_1 * weights[n]  # i <-> j for hole excitation
                    del cdag_j
                timer.print()

    # Save Green's function
    np.save("gf", gf)

if __name__ == '__main__':
    main()

