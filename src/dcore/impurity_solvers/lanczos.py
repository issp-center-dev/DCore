import numpy as np
from scipy.linalg import eigh_tridiagonal


def continued_fraction(a, b):
    """
    Compute the generalized continued fraction using the fundamental recurrence formulas.

      b_0 + a_1 / (b_1 + a_2 / (b_2 + a_3 / (b_3 + ...)))

    Parameters:
    - a (ndarray): a coefficients [0, a_1, a_2, ..., a_n].
    - b (ndarray): b coefficients [b_0, b_1, ..., b_n].

    a[0] is not used.

    Returns:
    - The value of the generalized continued fraction.
    """
    assert isinstance(a, np.ndarray)
    assert isinstance(b, np.ndarray)
    assert a.ndim == 1
    assert a.shape == b.shape

    # Initialize A_0 and B_0
    A_prev2 = b[0]  # A_0
    B_prev2 = 1     # B_0

    # Initialize A_1 and B_1 if there's more than one term
    if len(b) > 1:
        A_prev1 = b[1] * A_prev2 + a[1] * 1  # A_1 (using a[1], ignoring a[0])
        B_prev1 = b[1] * 1  # B_1
    else:
        return A_prev2  # If there's only one term, the result is b_0

    # Iteratively calculate A_n and B_n
    for n in range(2, len(b)):
        A_n = b[n] * A_prev1 + a[n] * A_prev2  # using a[n], ignoring a[0]
        B_n = b[n] * B_prev1 + a[n] * B_prev2

        # Update previous terms for next iteration
        # A_prev2, A_prev1 = A_prev1, A_n
        # B_prev2, B_prev1 = B_prev1, B_n

        # All terms are devided by B_n to avoid divergence
        A_prev2 = A_prev1 / B_n
        A_prev1 = A_n / B_n
        B_prev2 = B_prev1 / B_n
        B_prev1 = 1

    # Return A_n / B_n
    return A_prev1


class LanczosEigenSolver:
    def __init__(self, A, max_iter=None, tol=1e-6):
        """
        Initializes the Lanczos eigenvalue solver.

        Parameters:
        A (np.ndarray): A Hermitian matrix.
        max_iter (int): Maximum number of iterations. Default is size of the matrix.
        tol (float): Convergence tolerance. Default is 1e-10.
        """
        self.A = A
        self.n = A.shape[0]
        self.max_iter = max_iter if max_iter is not None else self.n
        self.tol = tol
        self.V = None  # Lanczos basis vectors

    def _tridiagonalize(self, dim, v0=None):
        """
        Tridiagonalizes the matrix A using the Lanczos method.

        Parameters:
        dim (int): Dimension of the tridiagonal matrix.
        v0 (np.ndarray or None): Initial vector for the Lanczos iteration. If None, a random vector is used.

        Returns:
        alphas (np.ndarray): Diagonal elements of the tridiagonal matrix.
        betas (np.ndarray): Off-diagonal elements of the tridiagonal matrix.
        """
        if v0 is None:
            v = np.random.randn(self.n)
        else:
            assert v0.shape == (self.n,)
            v = v0

        v /= np.linalg.norm(v)

        beta = 0
        alphas = np.zeros(dim)  # Diagonal elements of the tridiagonal matrix
        betas = np.zeros(dim - 1)  # Off-diagonal elements of the tridiagonal matrix

        self.V = np.zeros((self.n, dim), dtype=complex)  # Krylov subspace basis vectors
        # V = np.zeros((self.n, dim + 1), dtype=complex)  # Krylov subspace basis vectors
        V = self.V  # alias
        V[:, 0] = v

        for j in range(dim):
            w = self.A @ v
            # if j == 0:
            #     w = self.A @ v
            # else:
            #     w = self.A @ v - beta * V[:, j-1]

            alpha = np.vdot(w, v)  # Use np.vdot for complex numbers

            # w = w - alpha * v
            if j == 0:
                w -= alpha * v
            else:
                w -= alpha * v - beta * V[:, j-1]

            # Orthogonalize w with respect to all previous vectors in V using np.vdot for complex numbers
            for i in range(j + 1):
                w = w - np.vdot(V[:, i], w) * V[:, i]

            beta = np.linalg.norm(w)
            print(f"j={j}, alpha={alpha}, beta={beta}")

            alphas[j] = alpha.real  # Since A is Hermitian, alpha should be real

            if j < dim - 1:  # for next iteration
                betas[j] = beta  # Store off-diagonal element

                # Store the orthogonalized vector as the new Krylov subspace vector
                v = w / beta
                V[:, j+1] = v

        return alphas, betas

    def solve(self, k, ncv=None, v0=None, return_eigenvectors=True):
        assert k <= self.n

        if ncv is None:
            ncv = max(2*k + 1, 20)  # TODO
        ncv = min(ncv, self.n)

        alphas, betas = self._tridiagonalize(ncv, v0=v0)
        assert alphas.size == ncv
        assert betas.size == ncv - 1

        if return_eigenvectors:
            # Compute eigenvalues and eigenvectors of the tridiagonal matrix using eigh_tridiagonal
            eigvals, eigvecs_T = eigh_tridiagonal(alphas, betas)

            # Compute the eigenvectors of the original matrix A by multiplying Krylov subspace vectors with T's eigenvectors
            eigvecs_A = self.V[:, :len(eigvals)] @ eigvecs_T

            return eigvals[:k], eigvecs_A[:, :k]
        else:
            eigvals = eigh_tridiagonal(alphas, betas, eigvals_only=True)
            return eigvals[:k]

    def calc_gf(self, omega, v, ncv):
        """
        Computes the Green's function G(w) = <v| (w-H)^{-1} |v> using the Lanczos method.

        Parameters:
        omega (np.ndarray): complex frequencies.

        Returns:
        gf (np.ndarray): Green's function G(iω_k).
        """
        assert v.shape == (self.n,)

        if ncv is None:
            ncv = 20  # TODO
        ncv = min(ncv, self.n)

        v_norm = np.linalg.norm(v)
        v /= v_norm

        alphas, betas = self._tridiagonalize(ncv, v0=v)
        assert alphas.size == ncv
        assert betas.size == ncv - 1

        gf = np.zeros_like(omega, dtype=complex)
        for i, w in enumerate(omega):

            a = np.zeros(ncv + 1, dtype=float)
            b = np.zeros(ncv + 1, dtype=complex)

            # a[0] = 0  # not used
            a[1] = 1
            a[2:] = -betas**2

            # b[0] = 0
            b[1:] = w - alphas

            gf[i] = continued_fraction(a, b)

        return gf * v_norm**2



    # def solve(self, v0=None, return_eigenvectors=True):
    #     """
    #     Uses the Lanczos method to compute the top k eigenvalues and eigenvectors of the matrix A.

    #     Parameters:
    #     k (int): Number of eigenvalues and eigenvectors to compute.
    #     v0 (np.ndarray or None): Initial vector for the Lanczos iteration. If None, a random vector is used.

    #     Returns:
    #     eigvals (np.ndarray): Approximate top k eigenvalues of A.
    #     eigvecs (np.ndarray): Corresponding top k eigenvectors of A.
    #     """
    #     if v0 is None:
    #         v = np.random.randn(self.n)
    #     else:
    #         v = v0

    #     v = v / np.linalg.norm(v)  # Normalize the initial vector

    #     beta = 0
    #     alphas = np.zeros(self.max_iter)  # Diagonal elements of the tridiagonal matrix
    #     betas = np.zeros(self.max_iter - 1)  # Off-diagonal elements of the tridiagonal matrix

    #     # V = np.zeros((self.n, self.max_iter), dtype=complex)  # Krylov subspace basis vectors
    #     V = np.zeros((self.n, self.max_iter + 1), dtype=complex)  # Krylov subspace basis vectors
    #     V[:, 0] = v

    #     for j in range(self.max_iter):
    #         if j == 0:
    #             w = self.A @ v
    #         else:
    #             w = self.A @ v - beta * V[:, j-1]

    #         # Orthogonalize w with respect to all previous vectors in V using np.vdot for complex numbers
    #         for i in range(j + 1):
    #             w = w - np.vdot(V[:, i], w) * V[:, i]

    #         alpha = np.vdot(w, v)  # Use np.vdot for complex numbers
    #         w = w - alpha * v
    #         beta = np.linalg.norm(w)

    #         alphas[j] = alpha.real  # Since A is Hermitian, alpha should be real

    #         if j < self.max_iter - 1:
    #             betas[j] = beta  # Store off-diagonal element

    #         if beta < self.tol:
    #             alphas = alphas[:j+1]  # Truncate alphas to the correct size
    #             betas = betas[:j]  # Truncate betas to the correct size
    #             break

    #         # Store the orthogonalized vector as the new Krylov subspace vector
    #         v = w / beta
    #         V[:, j+1] = v

    #     if return_eigenvectors:
    #         # Compute eigenvalues and eigenvectors of the tridiagonal matrix using eigh_tridiagonal
    #         eigvals, eigvecs_T = eigh_tridiagonal(alphas, betas)

    #         # Compute the eigenvectors of the original matrix A by multiplying Krylov subspace vectors with T's eigenvectors
    #         eigvecs_A = V[:, :len(eigvals)] @ eigvecs_T

    #         return eigvals, eigvecs_A
    #     else:
    #         eigvals = eigh_tridiagonal(alphas, betas, eigvals_only=True)
    #         return eigvals
