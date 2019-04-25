Multiorbital model by a QMC solver
==================================

:download:`dmft_bethe.ini <dmft_bethe.ini>`

An interesting phenomena called spin-freezing transition occurs in multi-orbital models
[`P. Werner, E. Gull, M. Troyer and A. J. Millis, PRL 101, 166405 (2008) <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.101.166405>`_].
The spin-freezing phenomena is signaled by a peculiar frequency dependence of the self-energy: :math:`\mathrm{Im}\Sigma(i\omega_n) \propto \omega_n^{0.5}`.
In this tutorial, we solve the three-orbital model on a Bethe lattice to reproduce the results in the above reference.

.. literalinclude:: dmft_bethe.ini
   :language: ini
                              
Pre-process : ``dcore_pre``
---------------------------

.. First, we have to generate the lattice model as
.. The h5 file stores information of the model including the lattice structure, hopping parameters, interaction parameters.

We first generate a h5 file that is necessary for DMFT calculations.
The script ``dcore_pre`` is invoked for this purpose:

.. code-block:: bash

   $ dcore_pre dmft_bethe.ini

.. Then it outputs model HDF5 file (``bethe.h5``).
.. Parameters in [model] and [system] blocks are reads in the input file.

If succeeded, a h5 file named *seedname*.h5 (``bethe.h5`` in the present case) is generated.

DMFT loop : ``dcore``
---------------------

The DMFT loop is performed by ``dcore`` program.
In this tutorial, we use continuous-time QMC implementation of ALPS/CT-HYB.
The runtime of the impurity solver is set to 300 sec.
You should not use the Hubbard-I solver for a metallic system.
One can run the program with 24 MPI processes as follows.

.. code-block:: bash

   $ export MPIRUN="mpirun"
   $ dcore dmft_bethe.ini --np 24

.. Then it generates the result HDF5 file.

Any environment variables in ``command`` of the mpi section and ``exec_path`` of the impurity_solver section are expanded at runtime.
In addition, ``#`` in ``command`` of the mpi section is replaced by the number of processes specified at runtime (24 in this case).
In the above example, we defined the environment varible ``MPIRUN``.

Note that ``dcore`` must be lauched without using the mpirun command as
it launches MPI processes internally for heavy tasks.
The QMC solver is executed with the number of MPI processes as well.

.. We run this sample with 24 MPI processes.
Each self-consistent step takes around 5 min,
most of which is spent for solving an effective impurity problem by QMC.
40 iterations take around 200 min.
Results for the self-energy and Green's function in each iteration are accumulated into a h5 file named *seedname*.out.h5 (``bethe.out.h5`` in the present case).

One can check the convergence of DMFT iterations by using ``dcore_check`` program as follows.

.. code-block:: bash

   $ dcore_check dmft_bethe.ini

The extension can be pdf, eps, jpg, png, etc.

.. We can find the following standard output.

``dcore_check`` program prints the value of the chemical potential at each iteration on the standard output:

::

    @ Reading dmft_bethe.ini ...
  Loading Sigma_iw...
  Loading dc_imp and dc_energ...
    Total number of Iteration: 40

    Iter  Chemical-potential
    1 -0.751277456949
    2 0.434702172833
    3 1.24052366776
    4 1.87430591459
    5 2.3902083604
    6 2.80560950907
    7 3.15120273079
    8 3.39906953401
    9 3.59109395044
    10 3.77841000694
    11 3.94933085606
    12 4.04565623891
    13 4.14089310218
    14 4.21761388235
    15 4.19672207134
    16 4.22954996706
    17 4.21218913905
    18 4.23609175782
    19 4.29360816707
    20 4.30206162285
    21 4.30581583599
    22 4.31236778925
    23 4.3507968933
    24 4.34734048538
    25 4.34881158874
    26 4.35435696735
    27 4.32563489267
    28 4.3284996731
    29 4.32645567868
    30 4.31936571106
    31 4.30926807535
    32 4.34097777828
    33 4.3454259465
    34 4.29481990383
    35 4.32541999195
    36 4.33632971069
    37 4.38081434746
    38 4.38713422713
    39 4.34154680207
    40 4.38279281529
   Output check/sigma.dat
   Output check/sigma_ave.png
   Output check/iter_mu.dat
   Output check/iter_mu.png
   Output check/iter_sigma-ish0.png
   Output check/iter_sigma.dat

    Done

.. We also can see the imaginary-time self-energy at last seven iterations.

``dcore_check`` generates several figures as well as data files in text format.
For instance, ``check/iter_sigma-ish0.png`` shows how the renormalization factor converges for each orbital.

.. image:: check/iter_sigma-ish0.png
   :width: 800
   :align: center

If those results are not converged, one can restart DMFT iterations.

You can plot the data for positive Matsubara frequencies as follows (like Fig. 3 of PRL 101, 166405 (2008)) using
:download:`a gnuplot command file <plot.plt>`.


The reference data extracted from PRL 101, 166405 (2008) are available 
:download:`here <sigma-PRL101-166405.txt>`.
The plot should look like the following.

.. image:: sigma.svg
   :width: 600
   :align: center
