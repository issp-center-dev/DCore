.. _inputformat:

Input-file format and usage
===========================

Usage
-----

::

   $ dcore_pre input-file
   $ dcore input-file
   $ dcore_post input-file

They use the same input file.
   
Input-file format
-----------------

.. literalinclude:: ../tutorial/bethe-t2g/dmft_bethe.ini

[model] block
~~~~~~~~~~~~~

dcore_pre and dcore_post read this block.

============= ============= =========== ================================================================
Name          Type          Default     Description
============= ============= =========== ================================================================
t             Float         1.0         Transfer integral (Nearest neighbor) 
t'            Float         0.0         Transfer integral (Second nearest)
U             Float         0.0         On-site Coulomb potential
J             Float         0.0         On-site Hund potential
orbital_model String        single      Chosen from "single", "eg", "t2g", "full-d"
ncor          Integer       1           Number of correlation shell (Only wannier90).
lattice       String        chain       Chosen from "chain", "square", "cubic", "bethe", and "wannier90"
nelec         Float         1.0         Number of electrons per unit cell.
seedname      String        dcore      Name of the system.
                                        It should be the same as the seedname of wannier90.
cshell        Integer array [(0,1),...] Anguler momentum, and the number of orbitals of each
                                        correlation shell (Only wannier90). 
============= ============= =========== ================================================================

Prepare model for DCore.
Wannier90 as well as the following preset models:

* chain

* square

* cubic

* bethe
  Semicircular DOS with energy ranges [-2t:2t].  

* wannier90
  Read hopping parameters from the Wannier90 output.  

.. math::

   {\hat H} = \sum_{i j} \sum_{\alpha \beta} \sum_{\sigma}
   t_{i \alpha j \beta} c_{i \alpha \sigma}^\dagger c_{j \beta \sigma}
   +h.c. + {\hat H}_{\rm int}
  
.. image:: model.png
   :width: 700
   :align: center
        
[system] block
~~~~~~~~~~~~~~

dcore_pre and dcore read this block.

======= ======= ======= =================================================
Name    Type    Default Description
======= ======= ======= =================================================
beta    Float   1.0     Inverse temperature.
n_iw    Integer 2048    Number of Matsubara frequencies.
n_tau   Integer 10000   Number of imaginary-time points.
dc_type Integer -1      Type of double-counting correction.
fix_mu  Bool    False   Whether or not to use a fixed chemical potential.
nk      Integer 8       Number of *k* along each line
nk0     Integer 0       Number of *k* (Only wannier90)
nk1     Integer 0       Number of *k* (Only wannier90)
nk2     Integer 0       Number of *k* (Only wannier90)
======= ======= ======= =================================================
  
[impurity_solver] block
~~~~~~~~~~~~~~~~~~~~~~~

dcore and dcore_post read this block.
  
==== ======= =========== ===================================================
Name Type    Default     Description
==== ======= =========== ===================================================
name String  TRIQS/cthyb Name of impurity solver. Choosen from "TRIQS/cthyb"
                         "TRIQS/hubbrad-I", and "ALPS/cthyb".
==== ======= =========== ===================================================

**... and other parameters (Solver dependent)**

[control] block
~~~~~~~~~~~~~~~

dcore reads this block.

========= ======= ======= ============================================
Name      Type    Default Description
========= ======= ======= ============================================
max_step  Integer 100     Max number of SCF steps.
sigma_mix Float   0.5     Mixing parameter for self-energy.
delta_mix Float   0.5     Mixing parameter for hybridization function.
========= ======= ======= ============================================

[tool] block
~~~~~~~~~~~~

dcore_post reads this block.

============= ================= ===================== ======================================================
Name          Type              Default               Description
============= ================= ===================== ======================================================
bvec          Float array       [(1.0,0.0,0.0), (0.0, Reciplocal lattice vectors
                                1.0,0.0),(0.0,0.0,1.0
                                )]
nnode         Integer           2                     Number of node for the *k* path
nk_line       Integer           8                     Number of *k* along each line
knode         Sting Float array [(G,0.0,0.0,0.0),     The name and the fractional coordinate of each k-node.
                                 (X,1.0,0.0,0.0)]
omega_min     Float             -1                    Minimum value of real frequency
omega_max     Float             1                     Max value of real frequency
Nomega        Integer           100                   Number of real frequencies
broadening    Float             0.1                   An additional Lorentzian broadening
eta           Float             0.01                  Imaginary frequency shift
do_dos        Bool              False                 Whether or not calculate DOS
do_band       Bool              False                 Whether or not calculate band
============= ================= ===================== ======================================================
