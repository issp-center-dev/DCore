Tutorial with single-band 2D Hubbard model
==========================================

:download:`dmft_square.ini <square/dmft_square.ini>`

.. literalinclude:: square/dmft_square.ini
                              
Pre-process : ``dcore_pre``
---------------------------

First, we have to generate the lattice model as

.. code-block :: bash

   $ pydmf_pre dmft_square.ini

Then it outputs model HDF5 file (``square.h5``).

SCF cycle : ``dcore``
---------------------

The DMFT loop is performed by ``dcore`` program.
In this example, we use the Hubbard-I approximation for the impurity solver.

.. code-block :: bash

   $ pydmf dmft_square.ini

Then it generates the result HDF5 file.
We can check the convergence of the SCF cycle by using ``dcore_check`` program.

.. code-block :: bash

   $ pydmf_check dmft_square.ini

We can find the following standard output.

::

   Total number of Iteration: 10

   Iter  Chemical-potential
   1 0.0
   2 0.350471842936
   3 0.499576874374
   4 0.590504345604
   5 0.64888137395
   6 0.686284855018
   7 0.71112006693
   8 0.727786131759
   9 0.739550620867
   10 0.748282590054

We also can see the imaginary-time self-energy at last seven iterations.

.. image:: square/convergence.png
   :width: 500
   :align: center
   
Spectrul function : ``dcore_post``
----------------------------------

.. code-block :: bash

   $ pydmf_post dmft_square.ini
   $ gnuplot square_akw.gp

.. image:: square/akw.png
   :width: 700
   :align: center

.. code-block :: gnuplot

   gnuplot> set xlabel "Energy"
   gnuplot> set ylabel "DOS"
   gnuplot> plot "square_dos.dat" w l

.. image:: square/dos.png
   :width: 700
   :align: center
