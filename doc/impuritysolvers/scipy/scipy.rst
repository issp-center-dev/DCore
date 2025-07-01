Exact diagonalization solver: ``scipy/sparse``
==============================================

``scipy/sparse`` computes the Green's function by the exact diagonalization (ED) method. Small matrices are diagonalized by full diagonalization method implemented in ``scipy.linalg``, while large mtrices are treated by sparse solvers implemented in ``scipy.sparse.linalg``.

Features
--------

- Arbitrary temperature

- All interactions available in ``DCore`` are supported.

- Support thread parallelization by OpenMP (via NumPy and SciPy).

Install
-------

The following Python library needs to be installed:

- SciPy

Input parameters
----------------

Parameters for this solver are as follows:

::

    [impurity_solver]
    name = scipy/sparse
    n_bath{int} = 0
    # fit_gtol{float} = 1e-5
    # dim_full_diag{int} = 10000
    # particle_numbers{str} = [1,2,3]
    # weight_threshold{float} = 1e-6
    # n_eigen{int} = 100
    # eigen_solver{str} = eigsh
    # gf_solver{str} = bicgstab
    # check_n_eigen{bool} = False
    # check_orthonormality{bool} = False

The first two parameters are mandatory. The remaining parameters, which are commented out, are optional.

The table below shows the detailed description of the parameters.

.. csv-table::
    :header: "Name", "Type", "Default", "Description"

    "n_bath", "int", "0", "Number of bath sites. See :doc:`Pomerol solver<../pomerol/pomerol>` for details."
    "fit_gtol", "float", "1e-5", "Tolerance for the fitting of the hybridization function. See :doc:`Pomerol solver<../pomerol/pomerol>` for details."
    "dim_full_diag", "int", "10000", "Maximum dimension of the matrix to be diagonalized by full diagonalization method. If the matrix is larger than this value, it will be treated by sparse solver."
    "particle_numbers", "str", "auto", "Particle numbers to be considered. Allowed inputs are 'auto', 'all', or '[1,2,3]' etc. By default ('auto'), a minimal set of particle numbers is evaluated by prescreening calculations. If the set is known in advance, one can specify it in a list format, e.g., [1,2,3]. All particle numbers are considered if 'all' is chosen. Note that all states necessary to the Green's function calculation, namely, thermally occupied N-particle states and N±1 states, should be included."
    "weight_threshold", "float", "1e-6", "Threshold for the Boltzmann factor. States with Boltzmann factor smaller than this value are ignored in the Green's function calculation."
    "n_eigen", "int", "100", "Number of eigenvalues to be computed in each :math:`N`-particle subspace by the sparse solver."
    "eigen_solver", "str", "eigsh", "Name of the eigenvalue solver to be used. Available option is `'eigsh' <https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.eigsh.html#scipy.sparse.linalg.eigsh>`_ only."
    "gf_solver", "str", "bicgstab", "Name of the linear equation solver to be used for Green's function calculation. Available options are 'spsolve', 'bicg', 'bicgstab', 'cg', 'cgs', 'gmres', 'lgmres', 'minres', 'qmr', 'gcrotmk', 'tfqmr'. See `SciPy official document <https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#solving-linear-problems>`_ for details."
    "check_n_eigen", "bool", "True", "If True (default), calculation stops if ``n_eigen`` is not sufficient. Note that, even if this check is passed, ``n_eigen`` may be insufficient. One should check the convergence of the Green's function by increasing ``n_eigen``."
    "check_orthonormality", "bool", "True", "If True (default), orthonormality of the eigenvectors is checked. If the check fails, calculation stops."

Standard output
----------------

The standard output of the solver is saved in a file `stdout.log` in directory `work/imp_shell0_ite1`, where `shell0` and `ite1` stand for the `0`-th shell and the `1`-st iteration, respectively.

Some relevant output are explained below.

.. code-block:: text

    Particle-number conservation:
      N dim[N]
      0 1
      1 10
      2 45
      3 120
      4 210
      5 252
      6 210
      7 120
      8 45
      9 10
     10 1

This shows the dimension of the Hamiltonian matrix for each particle number.

When ``particle_numbers='auto'`` (default), a prescreening calculation is done to find a minimal set of particle numbers.

.. code-block:: text

    Prescreening -- Solving the eigenvalue problem...

    N = 0  (dim = 1)
     Time: 0m0.000s

    N = 1  (dim = 10)
     Time: 0m0.001s

    ...

    Finish the eigenvalue problem
    Time: 0m0.006s

    Summary of lowest energy state in each N block:
      N       dim    E_min     weight_rel
     ------------------------------------
      0         1    0.00e+00  1.4e-87
      1        10   -4.00e+00  3.3e-70
      2        45   -8.00e+00  7.7e-53
      3       120   -1.20e+01  1.8e-35
      4       210   -1.60e+01  4.2e-18
      5       252   -2.00e+01  1.0e+00
      6       210   -1.60e+01  4.2e-18
      7       120   -1.20e+01  1.8e-35
      8        45   -8.00e+00  7.7e-53
      9        10   -4.00e+00  3.3e-70
     10         1    1.82e-09  1.4e-87

    Particle numbers of thermally occupied states:
     [5]

    Particle numbers to be considered:
     [4 5 6]

The result of the prescreening is shown above. In this example, only N=5 includes thermally occupied states. For Green's function calculation, particle numbers N plus/minus 1 should be included as shown in the final line.

Next, the eigenvalue problems are solved for the specific particle numbers.

.. code-block:: text

    Solving the eigenvalue problem...

    N = 4  (dim = 210)
     full diagonalization
     Time: 0m0.009s

    N = 5  (dim = 252)
     full diagonalization
     Time: 0m0.010s

    N = 6  (dim = 210)
     full diagonalization
     Time: 0m0.006s

    Finish the eigenvalue problem
    Time: 0m0.216s

The time taken for solving the eigenvalue problem is shown for each particle number. The solver uses full diagonalization method for matrices smaller than ``dim_full_diag``, while it switches to sparse solver specified by ``eigen_solver`` for larger matrices. In this example, only full diagonalization is used.

.. code-block:: text

    Total eigenvalues computed:  672

    Save eigenvalues and eigenvectors in
     'eigenvalues.dat'
     'eigenvectors.dat'

    Number of initial states: 32
      N=4: 0 / 210
      N=5: 32 / 252
      N=6: 0 / 210

The first line shows the total number of eigenvalues computed. The last block shows the number of initial states that have the Boltzmann weight larger than `weight_threshold`. As expected from the prescreening, only N=5 contains initial states.

Then proceed to the calculation of the impurity Green's function.

.. code-block:: text

    Calculating impurity Green's function...

    Initial state 1/32  (N = 5)

     particle excitation: N + 1 = 6
      Use the Lehmann representation

     hole excitation: N - 1 = 4
      Use the Lehmann representation

    Initial state 2/32  (N = 5)

     particle excitation: N + 1 = 6
      Use the Lehmann representation

     hole excitation: N - 1 = 4
      Use the Lehmann representation

Lehmann representation is used when the dimension of :math:`N \pm 1`-particle states is smaller than ``dim_full_diag``. Otherwise, linear equations are solved by the sparse solver specified by ``gf_solver``. The output in this case is as follows:

.. code-block:: text

    Calculating impurity Green's function...

    Initial state 1/32  (N = 5)

     particle excitation: N + 1 = 6
      Solve linear equations
      Time: 0m3.816s

     hole excitation: N - 1 = 4
      Solve linear equations
      Time: 0m3.774s

    Initial state 2/32  (N = 5)

     particle excitation: N + 1 = 6
      Solve linear equations
      Time: 0m3.794s

     hole excitation: N - 1 = 4
      Solve linear equations
      Time: 0m3.782s

Output file
-----------

- **eigenvalues.dat**

  .. code-block:: text

    # dim = 1024
    # n_eigen = 100 (for each n)
    # N  E_i  Boltzmann_weight
    5  -2.00000000e+01  3.12500e-02
    5  -2.00000000e+01  3.12500e-02
    5  -2.00000000e+01  3.12500e-02
    5  -2.00000000e+01  3.12500e-02
    5  -2.00000000e+01  3.12500e-02
    5  -2.00000000e+01  3.12500e-02
    5  -2.00000000e+01  3.12500e-02
    ...

  This file contains the eigenvalues computed and the corresponding Boltzmann weights in ascending order. The numbers from left to right show the particle number, the eigen-energy, and the Boltzmann weight.

- **eigenvectors.dat**

  .. code-block:: text

    # The i-th digit from the left of the 10-dimensional Fock state
    # indicates whether (spin, site) state is occupied (1) or empty (0).
    # 'site' includes bath sites (b) after Wannier orbitals (w).
    #
    #   w1+ w2+ ... b1+ b2+ ... w1- w2- ... b1- b2- ...
    #
    # E=-2.00000000e+01, N=5
     |0000011111> 1.00000000e+00
    # E=-2.00000000e+01, N=5
     |0000111110> 1.00000000e+00
    # E=-2.00000000e+01, N=5
     |0001011101> 1.00000000e+00
    # E=-2.00000000e+01, N=5
     |0001111100> 1.00000000e+00
    # E=-2.00000000e+01, N=5
     |0110010011> 1.00000000e+00
    # E=-2.00000000e+01, N=5
     |0110110010> 1.00000000e+00
    ...

  This file contains the eigenvectors. The non-zero coefficients and the corresponding Fock states are shown for each eigenvector.


Benchmark
---------

to be updated.
