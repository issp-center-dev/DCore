import numpy as np
import scipy.sparse as sp
import argparse
import json
import sys


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


def sort_eigenvalues(vals, vecs):
    # Sort eigenvalues and eigenvectors in ascending order
    idx = np.argsort(vals)
    vals = vals[idx]
    vecs = vecs[:, idx]
    return vals, vecs


def convert_to_real(array, name="array"):
    # Check if H0 and U_ijkl are complex. If not, convert to real.
    if np.any(np.iscomplex(array)):
        print(f"{name} is complex.")
        return array
    else:
        print(f"{name} is real. -> Convert to real.")
        return array.real


def main():
    parser = argparse.ArgumentParser(description='Process a string argument.')
    parser.add_argument('filename', type=str, help='Input file in json format')
    args = parser.parse_args()

    # Load input parameters
    with open(args.filename, 'r') as f:
        params = json.load(f)

    assert isinstance(params, dict)
    print(params)

    n_flavors = params['n_flavors']
    n_sites = params['n_sites']
    beta = params['beta']
    n_eigen = params['n_eigen']
    n_iw = params['n_iw']
    flag_spin_conserve = params['flag_spin_conserve']

    # Load H0
    print("\nLoading H0 and U_ijkl")
    h0 = np.load(params['file_h0'])
    print(h0.shape)
    # assert h0.shape == (n_flavors, n_flavors)
    assert h0.shape == (2*n_sites, 2*n_sites)

    # Load U_ijkl
    umat = np.load(params['file_umat'])
    print(umat.shape)
    assert umat.shape == (2*n_sites, 2*n_sites, 2*n_sites, 2*n_sites)

    # Check if H0 and U_ijkl are complex. If not, convert to real.
    h0 = convert_to_real(h0, name="H0")
    umat = convert_to_real(umat, name="U_ijkl")

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

    print("\nHamiltonian:")
    print(" | type =", type(hamil))
    print(" | dtype =", hamil.dtype)

    # TODO: Use particle number conservation

    # Use full diagonalization if dim is small
    full_diagonalization = (n_eigen >= dim - 1)

    # Solve the eigenvalue problem
    if full_diagonalization:
        print("\nn_eigen >= dim\n  -> Use full diagonalization")
        n_eigen = dim
        E, eigvecs = np.linalg.eigh(hamil.toarray())
    else:
        print("\nn_eigen < dim\n  -> Use Lanczos method")
        E, eigvecs = sp.linalg.eigsh(hamil, k=n_eigen, which='SA')
        # E, eigvecs = sp.linalg.eigsh(hamil, k=n_eigen, which='LM')
        # E, eigvecs = sp.linalg.eigsh(hamil, k=n_eigen, which='SA', v0=np.ones(dim))
        # ‘SA’ : Smallest (algebraic) eigenvalues.

    assert E.shape == (n_eigen,)
    assert eigvecs.shape == (dim, n_eigen)

    # TODO: Sort eigenvalues and eigenvectors
    # The eigenvalues are not ordered in scipy.sparse.linalg.eigsh()
    # for complex Hermitian matrices since eigs() is called internally.

    # Sort eigenvalues and eigenvectors
    E, eigvecs = sort_eigenvalues(E, eigvecs)

    print("\nEigenvalues:")
    print(E)

    # Boltzmann factors
    weights = np.exp(-beta * (E - E[0]))
    weights /= np.sum(weights)
    print("\nWeights (Boltzmann factors / Z):")
    print(weights)

    weight_threshold = 1e-6
    n_initial_states = np.count_nonzero(weights > weight_threshold)
    print("\nNumber of initial states:", n_initial_states)

    if weights[-1] > weight_threshold and n_eigen < dim:
        sys.stdout.flush()
        print(f"\nWarning: n_eigen={n_eigen} may be too small: The weight for the highest-energy state computed is {weights[-1]}.", file=sys.stderr)

    # Save eigenvalues
    indexed_E = np.column_stack((np.arange(len(E)), E, weights))
    header = f"dim = {dim}\nn_eigen = {n_eigen}\ni  E_i  Boltzmann_weight"
    np.savetxt("eigenvalues.dat", indexed_E, fmt='%d %.8e %.5e', header=header)

    # Check if eigenvalues are ordered
    # assert np.all(np.diff(E) >= 0)

    # TODO: Save eigenvectors
    # print("Eigenvectors:")
    # print(eigvecs)

    # print("\nweights")
    # print(eigvecs * eigvecs.conj())

    # print("\nOrthonormality of eigenvectors:")
    overlap = eigvecs.conj().T @ eigvecs  # <m|n>
    # print("m n <m|n>")
    # for m, n in np.ndindex(n_eigen, n_eigen):
    #     print(m, n, overlap[m, n])

    ignore_orthonormality = True

    # Check orhotonormality of eigenvectors
    print("\nCheck orthonormality of eigenvectors")
    if not np.allclose(overlap, np.identity(n_eigen)):
        print("  not orthonormal")
        sys.stdout.flush()
        print("Warning: Eigenvectors are not orthonormal!", file=sys.stderr)
        if ignore_orthonormality:
            print("  -> Continue. Set ignore_orthonormality=False to abort calculation.", file=sys.stderr)
        else:
            print("  -> Exit. Set ignore_orthonormality=True to ignore this warning.", file=sys.stderr)
            sys.exit(1)
    else:
        print("  orthogonal")

    # Matsubara frequencies
    iws = 1j * (2 * np.arange(n_iw) + 1) * np.pi / beta

    # Calculate impurity Green's function
    gf = np.zeros((n_flavors, n_flavors, n_iw), dtype=complex)
    if full_diagonalization:
        for n in range(n_initial_states):

            # cdag_im[i, m] = <m|c_i^+|n>  for a given n
            cdag_im = np.empty((n_flavors, dim), dtype=complex)
            for i in range(n_flavors):
                cdag_im[i] = eigvecs.conj().T @ Cdag[i] @ eigvecs[:, n]

            # c_im[i, m] = <m|c_i|n>  for a given n
            c_im = np.empty((n_flavors, dim), dtype=complex)
            for i in range(n_flavors):
                c_im[i] = eigvecs.conj().T @ C[i] @ eigvecs[:, n]

            for l, iw in enumerate(iws):
                # ene_denom[m] = 1 / (iw - E_m + E_n)
                ene_denom_1 = 1 / (iw - E + E[n])

                # ene_denom[m] = 1 / (iw + E_m - E_n)
                ene_denom_2 = - ene_denom_1.conj()

                for i, j in np.ndindex(n_flavors, n_flavors):
                    if flag_spin_conserve:
                        if i // 2 != j // 2:  # skip different spins
                            continue

                    # gf[i, j, l] += np.einsum("m, m, m", cdag_im[i].conj().T, ene_denom, cdag_im[i]) * weights[n]
                    gf_1 = np.einsum("m, m, m", cdag_im[i].conj().T, ene_denom_1, cdag_im[j])
                    gf_2 = np.einsum("m, m, m", c_im[j].conj().T, ene_denom_2, c_im[i])
                    gf[i, j, l] += (gf_1 + gf_2) * weights[n]

    else:
        for n in range(n_initial_states):

            # loop for [particle excitation, hole excitation]
            for _Cdag, particle_excitation in [[Cdag, True], [C, False]]:
            # for _Cdag, particle_excitation in [[Cdag, True],]:
            # for _Cdag, particle_excitation in [[C, False],]:

                # < c_i c_j^+ > for particle excitation
                # < c_i^+ c_j > for hole excitation  (i <-> j later)
                for j in range(n_flavors):
                    # cdag_j = c_j^+ |n>  for a given (j, n)
                    cdag_j = _Cdag[j] @ eigvecs[:, n]

                    for l, iw in enumerate(iws):
                        # A = iw - H + E_n
                        # A = (iw + pm * E[n]) * sp.identity(dim) - pm * hamil
                        if particle_excitation:
                            # A = iw - H + E_n
                            A = (iw + E[n]) * sp.identity(dim) - hamil
                        else:
                            # A = iw + H - E_n
                            A = (iw - E[n]) * sp.identity(dim) + hamil
                        assert isinstance(A, sp.spmatrix)
                        assert A.shape == (dim, dim)

                        # Solve A |x> = c_j^+ |n>
                        x = sp.linalg.spsolve(A, cdag_j)
                        # x = np.linalg.solve(A.toarray(), cdag_j)
                        # x, _ = sp.linalg.lgmres(A, cdag_i)
                        assert x.shape == (dim,)

                        for i in range(n_flavors):
                            if flag_spin_conserve:
                                if i // 2 != j // 2:  # skip different spins
                                    continue

                            if i == j:
                                gf_1 = cdag_j.conj().T @ x
                            else:
                                cdag_i = _Cdag[i] @ eigvecs[:, n]
                                gf_1 = cdag_i.conj().T @ x
                                del cdag_i

                            # gf[i, j, l] += gf_1 * weights[n]
                            if particle_excitation:
                                gf[i, j, l] += gf_1 * weights[n]
                            else:
                                gf[j, i, l] += gf_1 * weights[n]  # i <-> j for hole excitation
                    del cdag_j

    # Save Green's function
    np.save("gf", gf)

if __name__ == '__main__':
    main()

