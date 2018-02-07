.. _howtocthyb:

How to choose parameters for TRIQS/cthyb
========================================


Basic parameters for QMC
------------------------

The accuracy of the TRIQS/cthyb calcuration depends on mainly the following
three parameters.

 * ``n_cycles{int}``

   **This parameter affects the statistical error.**
 
   The number of QMC cycles for comupting any quantities.
   The numerical noise can be reduced by increasing this parameter.
 
   When we use the MPI parallelism, the total number of QMC cycle
   ( ``number-of-processes * n_cycles``) affects the accuracy.
   Therefore,
   **n_cycles in inverse propotion to the number of MPI processes keeps the accuracy**.
 
 * ``n_warmup_cycles{int}``

   **This parameter affects the systematical error.**
 
   The number of QMC cycles for the thermalization before the above main calculation.
   If it is insufficient, the state to be computed may be different from the
   equiribrium state.

   **This parameter is independent from the MPI parallelism.**
 
 * ``length_cycle{int}``

   **This parameter affects the systematical error.**
 
   Each QMC cycles have subcycles.
   The length of this subcycle should be long enough to escape from the auto-correlation.

   **This parameter is independent from the MPI parallelism.**

The computational time is propotional to ``length_cycle*(n_cycles+n_warmup_cycles)``.
We will show the procedure to find appropriate values for these parameter.   
The sample case is the single-band Bethe lattice defined with the following input file.

::

   [model]
   seedname = bethe
   lattice = bethe
   norb = 1
   nelec = 1.0
   t = 1.0
   kanamori = [(4.0, 0.0, 0.0)]

   [system]
   beta = 40.0
   nk = 100

   [impurity_solver]
   name = TRIQS/cthyb
   n_warmup_cycles{int} = 5000
   n_cycles{int} = 500000
   length_cycle{int} = 50

   [control]
   max_step = 1

   [tool]
   omega_check = 30.0

Only the first DMFT cycle is required for checking the effect of the QMC parameter.
Therefore, we specify ``max_step = 1``.

Step 1 : ``n_cycles``
~~~~~~~~~~~~~~~~~~~~~

