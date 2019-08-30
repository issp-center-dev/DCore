CT-QMC with segment implementation: ``ALPS/CT-HYB-SEGMENT``
===========================================================

`ALPS/CT-HYB-SEGMENT solver <https://github.com/ALPSCore/CT-HYB-SEGMENT>`_ implements the CT-HYB QMC algorithm and is developed by H. Hafermann, P. Werner, E. Gull.
Only density-density interactions are taken into account.

A minimum working example is the following. Please be sure to turn on density_density option in model block (otherwise the simulation fails).
We assumed that the executable is installed at "$HOME/opt/CT-HYB/bin/alps_cthyb".
The solver terminates after MAX_TIME (300 seconds in the following example) or the total number of SWEEPS has been reached.


.. code-block:: ini

  [model]
  density_density = True
  [impurity_solver]
  name = ALPS/cthyb-seg
  exec_path{str} = $HOME/opt/CT-HYB/bin/alps_cthyb
  MAX_TIME{int} = 300
  cthyb.N_MEAS{int} = 50
  cthyb.THERMALIZATION{int}=1000
  cthyb.SWEEPS{int}=100000000
  cthyb.TEXT_OUTPUT{bool} = False

ALPS/CT-HYB-SEGMENT has many input parameters.
A complete list of the parameters can be obtained by the help command

.. code-block:: bash

  $HOME/opt/CT-HYB/bin/alps_cthyb --help

The following optional parameters may be useful:

.. code-block:: ini

  [impurity_solver]
  cthyb.MEASURE_nn{bool} = True   # static density-density correlation functions
  cthyb.MEASURE_nnw{bool} = True  # density-density correlation functions in frequency domain
  cthyb.MEASURE_nnt{bool} = True  # density-density correlation functions <n(0) n(t)>
  cthyb.MEASURE_g2w{bool} = True  # measure two-particle Green's function in frequency space

Please also refer to `the wiki page <https://github.com/ALPSCore/CT-HYB-SEGMENT/wiki/Changes-of-Parameters>`_ for some descriptions on parameters.

The DCore interface generates input files for ALPS/CT-HYB-SEGMENT into a working directory at work/imp_shell<ish>_ite<ite> (ish is the index of the shell and ite is the iteration).
Then, ALPS/CT-HYB-SEGMENT is executed in the working directory, and numerical results are stored there.
For example, the occupation number and the double occupancy are saved in the file 'observables.dat'.

Simple example: 2D Hubbard model
---------------------------

As an example, we consider the two-dimensional Hubbard model. The input file is given as below.
:download:`dmft_square_ctseg.ini <dmft_square_ctseg.ini>`

.. literalinclude:: dmft_square_ctseg.ini
   :language: ini

The momentum-resolved spectral functions obtained by executing the similar procedure given in the tutorial (See :doc:`here<../../../tutorial/square/square>`) are given as below.

.. image:: akw_ctseg.png
   :width: 700
   :align: center
