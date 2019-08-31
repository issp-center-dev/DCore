The first example: 2D Hubbard model
===================================

The first example is the two-dimensional Hubbard model.
We use the Hubbard-I approximation for solving the effective impurity problem and see the emergence of the Mott gap.
The input file is given below.

:download:`dmft_square.ini <dmft_square.ini>`

.. literalinclude:: dmft_square.ini
   :language: ini

Pre-process : ``dcore_pre``
---------------------------

.. First, we have to generate the lattice model as
.. The h5 file stores information of the model including the lattice structure, hopping parameters, interaction parameters.

We first generate an HDF5 file that is necessary for DMFT calculations.
The script ``dcore_pre`` is invoked for this purpose:

.. code-block:: bash

   $ dcore_pre dmft_square.ini

.. Then it outputs model HDF5 file (``square.h5``).
.. Parameters in [model] and [system] blocks are reads in the input file.

Then, an HDF5 file named *seedname*.h5 (``square.h5`` in the present case) will be generated.

DMFT loop : ``dcore``
---------------------

One can perform a DMFT self-consistent calculation with ``dcore`` program.
In this tutorial, we use the Hubbard-I solver just for simplicity.
One can run the program by

.. code-block:: bash

   $ dcore dmft_square.ini --np 1

with a single process. 

It takes several minutes. You may run it with MPI to reduce the computational time.
Results for the self-energy and Green's function in each iteration are accumulated into an HDF5 file named *seedname*.out.h5 (``square.out.h5`` in the present case).


One can check the convergence of DMFT iterations by using ``dcore_check`` program as follows.

.. code-block:: bash

   $ dcore_check dmft_square.ini

``dcore_check`` program prints the value of the chemical potential at each iteration on the standard output:

::

    @ Reading dmft_square.ini ...
  Loading Sigma_iw...
  Loading dc_imp and dc_energ...
    Total number of Iteration: 5

    Iter  Chemical-potential
    1 0.0
    2 0.141978800943
    3 0.597913733347
    4 0.700078346042
    5 0.742275654406
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

Spectral function : ``dcore_post``
----------------------------------
We can calculate the density of states and the momentum-dependent single-particle excitations using ``dcore_post`` program.
In the Hubbard-I solver, the self-energy on the real-frequency axis can be directly computed (no analytical continuation is required).
Hence, the impurity problem is solved once more in ``dcore_post``.

The calculation is done by the following command:

.. code-block:: bash

   $ dcore_post dmft_square.ini --np 1

After finishing the calculation,
``square_akw.dat``, ``square_akw.gp`` and ``square_dos.dat`` are generated.
The data of momentum-resolved spectral functions are output into ``square_akw.dat``.
We can easily plot the result by using the script file ``square_akw.gp`` for gnuplot:

.. code-block:: bash

   $ gnuplot square_akw.gp

In the plot shown below, the left and right panels correspond to up and down spins, respectively.

.. image:: akw.png
   :width: 700
   :align: center

The result for the density of states is output into ``square_dos.dat``.
We can plot it using gnuplot as follows:

.. code-block:: gnuplot

   set xlabel "Energy"
   set ylabel "DOS"
   plot "square_dos.dat" w l

.. image:: dos.png
   :width: 700
   :align: center

Another impurity solver: CTHYB-SEG
---------------------------

The input file for ALPS/cthyb-seg is given as below.

:download:`dmft_square_ctseg.ini <dmft_square_ctseg.ini>`

.. literalinclude:: dmft_square_ctseg.ini
   :language: ini

The momentum-resolved spectral functions obtained by executing the similar procedure given in this tutorial are given as below.

.. image:: akw_ctseg.png
   :width: 700
   :align: center
