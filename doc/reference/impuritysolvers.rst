.. _impuritysolvers:

Impurity solvers
==================

``null``
~~~~~~~~~~~~~~~~~~~~

This impurity solver returns zero self-energy, which allows to check the non-interacting limit.


``TRIQS/hubbard-I``
~~~~~~~~~~~~~~~~~~~~
This impurity solver implements the Hubbard-I approximation.


``pomerol``
~~~~~~~~~~~~~~~~~~~~
This impurity solver implements the Hubbard-I approximation.


``ALPS/CT-HYB``
~~~~~~~~~~~~~~~~~~~~
This impurity solver implements the CT-HYB QMC algorithm and is used by DCore developers on a regular basis.


``ALPS/CT-HYB-SEGMENT``
~~~~~~~~~~~~~~~~~~~~~~~~~~
This impurity solver implements the CT-HYB segment algorithm.
Only density-density interactions are taken into account.
The off-diagonal components of the Green's functions are not computed correctly.


``TRIQS/cthyb``
~~~~~~~~~~~~~~~~~~~~
This impurity solver implements the CT-HYB QMC algorithm.



