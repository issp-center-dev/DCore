Exact diagonalization solver: ``HPhi``
======================================

``HPhi`` is a numerical exact-diagonalization / Lanczos program for quantum
lattice models.
``DCore`` provides an interface to ``HPhi`` to solve the effective impurity
problem of DMFT with a discretized hybridization function (a finite number of
bath sites).

Features
--------

- Arbitrary temperature (the one-body Green's function is evaluated on the
  Matsubara axis).

- All interactions available in ``DCore`` are supported (the full
  :math:`U_{ijkl}` tensor is passed to ``HPhi``).

- With ``n_bath = 0`` the solver reduces to the Hubbard-I approximation.

Install
-------

The following program needs to be installed:

- `HPhi <https://github.com/issp-center-dev/HPhi>`_

Make sure that the ``HPhi`` executable is available, and pass its path through
the ``exec_path`` parameter below.

How to use
----------

Mandatory parameters:

::

    [impurity_solver]
    name = HPhi
    exec_path{str} = /install_directory/bin/HPhi

Optional parameters:

::

    [impurity_solver]
    n_bath{int} = 0        # number of bath sites per spin-orbital (0 = Hubbard-I)
    exct{int} = 1          # number of low-lying eigenstates used for the
                           # finite-temperature Green's function
    fit_gtol{float} = 1e-5 # tolerance of the bath-hybridization fitting

Notes
-----

- ``n_bath`` controls the number of discretized bath sites that approximate the
  hybridization function. Larger ``n_bath`` improves the approximation but
  enlarges the many-body Hilbert space, whose dimension grows exponentially
  with the total number of sites ``n_orb + n_bath``.

- ``exct`` sets how many low-lying eigenstates are included in the
  finite-temperature average. It should be large enough that the Boltzmann
  weight of the highest computed state is negligible at the target temperature.

- ``HPhi`` requires the number of MPI processes to be a power of four for the
  eigenenergy calculation. If the requested number of processes is not a power
  of four, ``DCore`` automatically falls back to the largest power of four for
  that step and prints a warning.
