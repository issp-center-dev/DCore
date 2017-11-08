pydmft
======

Perform DMFT.

Usage
-----

::

   $ pydmft dmft.ini seedname

Input-file format
-----------------

Example
~~~~~~~

.. literalinclude:: ../tutorial/bethe-t2g/dmft.ini

Details
~~~~~~~

* [system] block

======= ======= ======= =================================================
Name    Type    Default Description
======= ======= ======= =================================================
beta    Float   1.0     Inverse temperature.
n_iw    Integer 2048    Number of Matsubara frequencies.
n_tau   Integer 10000   Number of imaginary-time points.
dc_type Integer -1      Type of double-counting correction.
fix_mu  Bool    False   Whether or not to use a fixed chemical potential.
======= ======= ======= =================================================
  
* [impurity_solver] block

==== ======= =========== ===================================================
Name Type    Default     Description
==== ======= =========== ===================================================
name String  TRIQS/cthyb Name of impurity solver. Choosen from "TRIQS/cthyb"
                         "TRIQS/hubbrad-I", and "ALPS/cthyb".
==== ======= =========== ===================================================

**... and other parameters (Solver dependent)**


* [control] block

========= ======= ======= ============================================
Name      Type    Default Description
========= ======= ======= ============================================
max_step  Integer 100     Max number of SCF steps.
sigma_mix Float   0.5     Mixing parameter for self-energy.
delta_mix Float   0.5     Mixing parameter for hybridization function.
========= ======= ======= ============================================

