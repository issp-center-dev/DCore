How to integrate your own solver
================================

You can use your own impurity solver from DCore by wrapping the executable as appropriate.
The null solver serves as a template of wrapper code.

    :download:`null_solver.py <../../src/dcore/impurity_solvers/null_solver.py>`

Follow the instructions below.

..
    .. literalinclude:: ../../src/dcore/impurity_solvers/null_solver.py
       :language: python

Copy the template
-----------------

In directory **python/impurity_solvers**, copy the template

::

    $ cp null_solver.py your_solver.py

Edit the template
-----------------

-   Change the class name

    .. code-block:: python

        class NullSolver(SolverBase):

    Let's say you changed the class name into **YourSolver**.

-   Set the solver name

    .. code-block:: python

        def name(self):
            return "null"

-   Set input to the solver

    .. code-block:: python

        # (1) Set configuration for the impurity solver
        # input:
        #   self.beta
        #   self.set_G0_iw
        #   self.u_mat
        #
        # Additionally, the following variables may be used:
        #   self.n_orb
        #   self.n_flavor
        #   self.gf_struct
        #   self.n_tau
        #   self.use_spin_orbit

        # (1a) If H0 is necessary:
        # Non-interacting part of the local Hamiltonian including chemical potential
        # Make sure H0 is hermite.
        # Ordering of index in H0 is spin1, spin1, ..., spin2, spin2, ...
        H0 = extract_H0(self._G0_iw, self.block_names)

        # (1b) If Delta(iw) and/or Delta(tau) are necessary:
        # Compute the hybridization function from G0:
        #     Delta(iwn_n) = iw_n - H0 - G0^{-1}(iw_n)
        # H0 is extracted from the tail of the Green's function.
        self._Delta_iw = delta(self._G0_iw)
        Delta_tau = make_block_gf(GfImTime, self.gf_struct, self.beta, self.n_tau)
        for name, block in self._Delta_iw:
            Delta_tau[name] << Fourier(self._Delta_iw[name])

        # (1c) Set U_{ijkl} for the solver
        # for i, j, k, l in product(range(self.n_flavors), repeat=4):
        #     self.u_mat[i, j, k, l]

    Here, you generate **all** necessary input files to run your program.

-   Run the solver

    .. code-block:: python

        # (2) Run a working horse
        with open('./output', 'w') as output_f:
            launch_mpi_subprocesses(mpirun_command, [exec_path, 'input.ini'], output_f)

    The second argument to the function ``launch_mpi_subprocesses`` is the actual command and arguments used for invoking MPI processes from **DCore**.
    These options and arguments to the solver must be given as a list (refer to ``subprocess`` module in Python).

-   Convert output of the solver

    .. code-block:: python

        # (3) Copy results into
        #   self._Sigma_iw
        #   self._Gimp_iw

    Read output files generated by the solver, and set them into **DCore** variables.
    The self-energy and the Green's function are stored as ``BlockGf`` class of **TRIQS** library.
    See other wrappers or `the TRIQS documentation <https://triqs.github.io/triqs/1.4/reference/gfs/py/contents.html>`_.

Register your solver
--------------------

Finally, you register your own solver to **DCore**.
Edit **python/impurity_solvers/__init__.py**
to import your class

.. code-block:: python

    from .your_solver import YourSolver

and add it to the dictionary variable ``solver_classes`` as

.. code-block:: python

    solver_classes = {
        ...
        'your_solver': YourSolver,
    }

Then, you can invoke your solver from **DCore** by

.. code-block:: python

    [impurity_solver]
    name = your_solver
