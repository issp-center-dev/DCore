.. _howtoalpscthyb:

CT-QMC: ``ALPS/CT-HYB``
=======================

`ALPS/CT-HYB solver <https://github.com/ALPSCore/CT-HYB>`_ implements the CT-HYB QMC algorithm and is developed by H. Shinaoka, E. Gull and P. Werner.
The spin-orbit coupling is supported by default.
One can call ALPS/CT-HYB from DCore through the "ALPS/cthyb" interface.
For the moment, this interface is the first choice for QMC solvers in DCore
as it is used by developers of DCore on a regular basis.

The old triqs_interface for DCore v1.x has retired.
The new interface involved in DCore 2.x requires only the installation of ALPS/CT-HYB 1.x.
Please follow `the instruction on the official site <https://github.com/ALPSCore/CT-HYB/wiki>`_.

A minimum working example is the following.
We assumed that the executable is installed at "$HOME/opt/CT-HYB/bin/hybmat".
The solver runs for 300 seconds.

.. code-block:: ini

  [impurity_solver]
  name = ALPS/cthyb
  timelimit{int} = 300
  exec_path{str} = $HOME/opt/CT-HYB/bin/hybmat

.. note::

  When ``timelimit{int}`` is too small to generate a sufficient number of the Monte Carlo samples, the solver crashes with a segmentation fault.
  If you treat with a large system (e.g., many orbitals), you may need large ``timelimit{int}``.

ALPS/CT-HYB have many input parameters.
Please refer to `the reference manual <https://github.com/ALPSCore/CT-HYB/wiki/Input-parameters>`_ for a list of available input parameters.
The DCore interface supports support all of them.
For instance, one can enable verbose model as follows.
Note that one must specify type for input parameters of ALPS/CT-HYB.

.. code-block:: ini

  [impurity_solver]
  verbose{int} = 1

The DCore interface generates input files for ALPS/CT-HYB into a working directory at work/imp_shell<ish>_ite<ite> (ish is the index of the shell and ite is the iteration).
Then, ALPS/CT-HYB is excecuted in the working directory.
