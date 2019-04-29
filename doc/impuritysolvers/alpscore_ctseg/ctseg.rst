CT-QMC with segment implementation: ``ALPS/CT-HYB-SEGMENT``
===========================================================

This impurity solver implements the CT-HYB-SEGMENT algorithm.
Only density-density interactions are taken into account.
The off-diagonal components of the Green's functions are not computed correctly.

A minimum working example is the following.
We assumed that the executable is installed at "$HOME/opt/CT-HYB/bin/alps_cthyb".
The solver runs for 300 seconds.


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

ALPS/CT-HYB-SEGMENT have many input parameters.  
Please refer to `the reference manual <https://github.com/ALPSCore/CT-HYB-SEGMENT/wiki/Changes-of-Parameters>`_ for a list of available input parameters.

The DCore interface generates input files for ALPS/CT-HYB-SEGMENT into a working directory at work/imp_shell<ish>_ite<ite> (ish is the index of the shell and ite is the iteration).
Then, ALPS/CT-HYB-SEGMENT is excecuted in the working directory.
