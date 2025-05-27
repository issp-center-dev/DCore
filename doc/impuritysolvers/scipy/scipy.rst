Exact diagonalization solver: ``scipy/sparse``
==============================================

``scipy/sparse`` computes the Green's function by the exact diagonalization (ED) method. Small matrices are diagonalized by full diagonalization method implemented in ``scipy.linalg``, while large mtrices are treated by sparse solvers implemented in ``scipy.sparse.linalg``.

Features
--------

- Arbitrary temperature

- All interactions available in ``DCore`` are supported.

Install
-------

The following Python library needs to be installed:

- SciPy

How to use
----------

Parameters for this solver are as follows:

::

    [impurity_solver]
    name = scipy/sparse
    n_bath{int} = 0
    # fit_gtol{float} = 1e-5
    # dim_full_diag{int} = 10000
    # particle_numbers{str} = None  # e.g., [1,2,3]
    # n_eigen{int} = 100
    # eigen_solver{str} = eigsh
    # gf_solver{str} = bicgstab
    # check_n_eigen{bool} = True
    # check_orthonormality{bool} = True

The first two parameters are mandatory. The remaining parameters, which are commented out, are optional. The values shown represent their default settings.

Description of the parameters are given below:

- ``n_bath``: Number of bath sites. See :doc:`Pomerol solver<../pomerol/pomerol>` for details.

- ``fit_gtol{float}``: Tolerance for the fitting of the hybridization function. See :doc:`Pomerol solver<../pomerol/pomerol>` for details.

- ``dim_full_diag{int}``: Maximum dimension of the matrix to be diagonalized by full diagonalization method. If the matrix is larger than this value, it will be treated by sparse solver.

- ``particle_numbers{str}``: Particle numbers to be considered in a list format, e.g., `[1, 2, 3]`. If None (default), all particle numbers are considered. This parameter is useful to reduce the calculation time. Note that all states necessary to the Green's function calculation, namely, thermally occupied :math:`N`-particle states and :math:`N\pm 1` states, should be included.

- ``n_eigen{int}``: Number of eigenvalues to be computed by the sparse solver.

- ``eigen_solver{str}``: Name of the eigenvalue solver to be used. Available option is `'eigsh' <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh>`_ only.

- ``gf_solver{str}``: Name of the linear equation solver to be used for Green's function calculation. Available options are 'spsolve', 'bicg', 'bicgstab', 'cg', 'cgs', 'gmres', 'lgmres', 'minres', 'qmr', 'gcrotmk', 'tfqmr'. See `official document <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems>`_ for details.

- ``check_n_eigen{bool}``: If True (default), calculation stops if ``n_eigen`` is not sufficient. Note that, even if this check is passed, ``n_eigen`` may be insufficient. One should check the convergence of the Green's function by increasing ``n_eigen``.

- ``check_orthonormality{bool}``: If True (default), orthonormality of the eigenvectors is checked. If the check fails, calculation stops.

Example
-------

to be updated.

Benchmark
---------

to be updated.
