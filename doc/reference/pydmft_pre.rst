.. _pydmftpre:

pydmft_pre
==========

Prepare model for pyDMFT.
Wannier90 as well as the following preset models:

* Chain lattice

  .. math::

     {\hat H} = \sum_{i} \sum_{\alpha} \sum_{\sigma}
     t c_{i+1 \alpha \sigma}^\dagger c_{i \alpha \sigma}
     +h.c. + {\hat H}_{\rm int}
  
* Square lattice

  .. math::

     {\hat H} = \sum_{i,j} \sum_{\alpha} \sum_{\sigma}
     t c_{i \alpha \sigma}^\dagger c_{i \alpha \sigma}
     +h.c. + {\hat H}_{\rm int}
  
* Cubic lattice

  .. math::

     {\hat H} = \sum_{i} \sum_{\alpha} \sum_{\sigma}
     t c_{i \alpha \sigma}^\dagger c_{i \alpha \sigma}
     +h.c. + {\hat H}_{\rm int}
  
* Bethe lattice

Usage
-----

::

   $ pydmft_pre input-file

Input-file format
-----------------

Example
~~~~~~~

.. literalinclude:: ../tutorial/bethe-t2g/model.in

Details
~~~~~~~

* [model] block

============= ============= =========== ================================================================
Name          Type          Default     Description
============= ============= =========== ================================================================
t             Float         1.0         Transfer integral (Nearest neighbor) 
tp            Float         0.0         Transfer integral (Second nearest)
U             Float         0.0         On-site Coulomb potential
J             Float         0.0         On-site Hund potential
orbital_model String        single      Chosen from "single", "eg", "t2g", "full-d"
nk            Integer       8           Number of *k* along each line
nk0           Integer       0           Number of *k* (Only wannier90)
nk1           Integer       0           Number of *k* (Only wannier90)
nk2           Integer       0           Number of *k* (Only wannier90)
ncor          Integer       1           Number of correlation shell.
lattice       String        chain       Chosen from "chain", "square", "cubic", "bethe", and "wannier90"
nelec         Float         1.0         Number of electrons per unit cell.
seedname      String        pydmft      Name of the system.
                                        It should be the same as the seedname of wannier90.
cshell        Integer array [(0,1),...] Anguler momentum, and the number of orbitals of each
                                        correlation shell. 
============= ============= =========== ================================================================

.. * equiv-0, equiv-1, equiv-2, differ at each shell, Equivalence of shell.

