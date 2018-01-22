Tutorial with single-band 2D Hubbard model
==========================================

:download:`dmft_square.ini <square/dmft_square.ini>`

.. literalinclude:: square/dmft_square.ini

Pre-process : ``dcore_pre``
---------------------------

.. First, we have to generate the lattice model as
.. The h5 file stores information of the model including the lattice structure, hopping parameters, interaction parameters.

We first generate a h5 file that is necessary for DMFT calculations.
The script ``dcore_pre`` is invoked for this purpose:

.. code-block:: bash

   $ dcore_pre dmft_square.ini

.. Then it outputs model HDF5 file (``square.h5``).
.. Parameters in [model] and [system] blocks are reads in the input file.
      
If succeeded, a h5 file named *seedname*.h5 (``square.h5`` in the present case) is generated.

DMFT loop : ``dcore``
---------------------

The DMFT loop is performed by ``dcore`` program.
In this tutorial, we use the Hubbard-I solver just for simplicity.
One can run the program by

.. code-block:: bash

   $ dcore dmft_square.ini

.. Then it generates the result HDF5 file.
      
It takes several minutes. You may run it with MPI to save time.
Results for the self-energy and Green's function in each iteration are accumulated into a h5 file named *seedname*.out.h5 (``square.out.h5`` in the present case).

One can check convergence of DMFT loop by using ``dcore_check`` program.
You can run it with the following command, if X window system is available:

.. code-block:: bash

   $ dcore_check dmft_square.ini

If X window is not available or you prefer plotting in a file, use ``--output`` option to specify output file name

.. code-block:: bash

   $ dcore_check dmft_square.ini --output=check.pdf

The extension can be eps, jpg, etc.

.. We can find the following standard output.

``dcore_check`` program prints the value of the chemical potential at each iteration on the standard output:

::

   Total number of Iteration: 7

   Iter  Chemical-potential
   1 0.0
   2 0.141978800943
   3 0.46478279315
   4 0.637322531819
   5 0.646372779247
   6 0.680315738712
   7 0.708829559685

.. We also can see the imaginary-time self-energy at last seven iterations.
      
``dcore_check`` also plots the self-energy for the last seven iterations in Matsubara-frequency domain.

.. image:: square/convergence.png
   :width: 800
   :align: center

If those results are not converged, one can repeat the DMFT iteration using the same ini file. ``dcore`` program automatically finds results in the previous run and continue iterations.

Spectral function : ``dcore_post``
----------------------------------
We can calculate the density of states and the momentum-dependent single-particle excitations using ``dcore_post`` program.
For Hubbard-I solver, the self-energy is first calculated in ``dcore_post``.
The calculation is done by the following command:

.. code-block:: bash

   $ dcore_post dmft_square.ini

After finishing the calculation,
``square_akw.dat``, ``square_akw.gp`` and ``square_dos.dat`` are generated.
The data of spectral function are output into ``square_akw.dat``.
By using ``square_akw.gp``, we can easily plot the result:

.. code-block:: bash

   $ gnuplot square_akw.gp

.. image:: square/akw.png
   :width: 700
   :align: center

The data of density of states are output into ``square_dos.dat``.
We can plot the results by using gnuplot as follows:

.. code-block:: gnuplot

   gnuplot> set xlabel "Energy"
   gnuplot> set ylabel "DOS"
   gnuplot> plot "square_dos.dat" w l

.. image:: square/dos.png
   :width: 700
   :align: center
