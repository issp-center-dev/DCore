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
    exct_weight_threshold{float} = 1e-3   # threshold for the "exct too small" warning
    n_procs_per_hphi{int} = 1             # MPI ranks per HPhi run in the Gf step

Notes
-----

- ``n_bath`` controls the number of discretized bath sites that approximate the
  hybridization function. Larger ``n_bath`` improves the approximation but
  enlarges the many-body Hilbert space, whose dimension grows exponentially
  with the total number of sites ``n_orb + n_bath``.

- ``exct`` sets how many low-lying eigenstates are included in the
  finite-temperature average. It should be large enough that the Boltzmann
  weight of the highest computed state is negligible at the target temperature.
  See `Choosing exct (degenerate ground states)`_ below.

- ``HPhi`` requires the number of MPI processes to be a power of four for the
  eigenenergy calculation. If the requested number of processes is not a power
  of four, ``DCore`` automatically falls back to the largest power of four for
  that step and prints a warning.

- The one-body Green's function is computed from many independent HPhi runs
  (one per excitation). These runs are parallelised in two levels: ``n_outer``
  runs execute concurrently, each using ``n_procs_per_hphi`` MPI ranks, with
  ``n_outer = np / n_procs_per_hphi``. The default ``n_procs_per_hphi = 1`` runs
  each HPhi serially and executes ``np`` runs concurrently (the previous
  behaviour). Increase ``n_procs_per_hphi`` (a power of four) when a single HPhi
  run is heavy, e.g. with many bath sites; it trades concurrency for per-run MPI
  speed. The total ``np`` is shared between the two levels.


Choosing exct (degenerate ground states)
----------------------------------------

``HPhi`` builds the finite-temperature Green's function from the lowest
``exct`` eigenstates, weighting eigenstate :math:`n` by the Boltzmann factor
:math:`e^{-\beta (E_n - E_0)}`.  If ``exct`` is too small, states that are still
thermally relevant are left out of the trace and the resulting self-energy is
**wrong**.

.. warning::

   The default ``exct = 1`` is unsafe for multi-orbital models.  Multi-orbital
   impurities very often have a **degenerate ground state** (for example, one
   electron in two degenerate orbitals is four-fold degenerate).  With
   ``exct = 1`` only one member of the multiplet enters the thermal trace, which
   **breaks the orbital symmetry** and yields a spurious, asymmetric self-energy
   with non-zero off-diagonal components — even though the exact answer is
   symmetric and diagonal.

   Set ``exct`` large enough to cover the whole degenerate ground multiplet
   *and* every excited state with a non-negligible Boltzmann weight at the
   target temperature.  When in doubt, increase ``exct`` until the result stops
   changing (for a tiny problem you may simply use the full Hilbert-space
   dimension :math:`4^{\,\mathrm{n\_orb}+\mathrm{n\_bath}}`).

To help catch this, ``DCore`` inspects the computed eigenenergies after the
eigenvalue step and prints a warning when the highest computed state still
carries a Boltzmann weight above ``exct_weight_threshold`` (default ``1e-3``) —
i.e. when states above the ``exct`` cutoff are likely missing from the thermal
trace.  The warning reports the ground-state degeneracy and the remaining
weight; increase ``exct`` until it disappears.  This is a thermal
**convergence** check, so it also fires for non-degenerate ground states with
low-lying thermally-populated excitations.

.. tip::

   Because both ``HPhi`` and the ``scipy/sparse`` solver are exact-diagonalization
   solvers, they must give the same self-energy for the same impurity model.
   Cross-checking ``HPhi`` against ``scipy/sparse`` on a small model is a good
   way to confirm that ``exct`` (and other settings) are adequate.


Running HPhi in parallel
------------------------

``DCore`` launches ``HPhi`` through the MPI command defined in the ``[mpi]``
section, where ``#`` is replaced by the number of processes passed to
``dcore`` (``dcore --np N``)::

    [mpi]
    command = mpirun -np #

The HPhi solver runs in two MPI phases:

1. **Eigenvalue step** -- a single ``HPhi`` invocation that uses *all* ``np``
   processes. ``HPhi`` requires this count to be a power of four; ``DCore``
   automatically rounds down (with a warning) if it is not.

2. **Green's-function step** -- many independent ``HPhi`` runs (one per
   excitation). These are distributed by the two-level layout described above:
   ``n_outer = np / n_procs_per_hphi`` runs execute concurrently, each launched
   with ``n_procs_per_hphi`` MPI ranks. ``n_procs_per_hphi`` must be a power of
   four.

.. note::

   With an **MPI-only** build of ``HPhi`` (one that must be started through the
   MPI launcher), keep ``n_procs_per_hphi`` :math:`\geq` the smallest value your
   launcher accepts (typically 1 via ``mpirun -np 1`` / ``srun -n 1``), and
   prefer launching each Green's-function run through the launcher rather than
   bare. Running many bare MPI executables concurrently can fail due to MPI
   runtime resource clashes; going through the launcher (``n_procs_per_hphi``
   selects ``n_inner`` ranks per run) avoids this.

On a SLURM cluster (e.g. ISSP System B "ohtaka", 128 cores/node), set the
launcher to ``srun`` and choose the split with ``n_procs_per_hphi``. A minimal
job script that puts each Green's-function run on 4 ranks and runs
``128 / 4 = 32`` of them concurrently per node looks like::

    #!/bin/sh
    #SBATCH -p i8cpu          # interactive queue for testing (adjust for production)
    #SBATCH -N 1
    #SBATCH -n 128
    #SBATCH -t 0:30:00

    module load <your HPhi / python environment>

    # in the input file:
    #   [mpi]
    #   command = srun -n #
    #   [impurity_solver]
    #   name = HPhi
    #   exec_path{str} = /path/to/mpi/HPhi
    #   n_procs_per_hphi{int} = 4
    dcore --np 128 input.ini

Here ``dcore --np 128`` gives ``np = 128`` (a power of four) for the eigenvalue
step, and the Green's-function step runs ``32`` HPhi instances of ``4`` ranks
each. Start small (e.g. ``--np 4`` with ``n_procs_per_hphi = 4``) on the
interactive queue to validate the setup before scaling up.
