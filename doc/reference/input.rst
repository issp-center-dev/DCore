.. _inputformat:

Input-file format and usage
===========================

Usage
-----

The following programs can read the same input file.

Pre-processing : ``dcore_pre``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program generates model HDF5 file (*seedname*.h5).
It must be executed before the main program, ``dcore``, runs.
This program reads ``[model]`` and ``[system]`` block.

::

   $ dcore_pre input-file

Main program : ``dcore``
~~~~~~~~~~~~~~~~~~~~~~~~

This program performs DMFT cycle and output the self energy etc. into a HDF
file (*seedname*.out.h5).
This program reads ``[model]``, ``[system]``, ``[impurity-solver]`` and ``[control]`` block.

::

   $ dcore input-file

Convergence-check : ``dcore_check``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program can be used for checking the convergence of the SCF-cycle.
This program reads ``[model]`` and ``[tool]`` block.

::

   $ dcore_check input-file

``dcore_check`` shows the history of the chemical potential and the
first component of the self energy at imaginary frequency, :math:`\Sigma_{0 0}(i \omega_n)`
at the last seven iterations.

.. image:: ../tutorial/square/convergence.png
   :width: 500
   :align: center

Post-processing : ``dcore_post``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program compute the total and *k*\ -resolved spectrul function from the outputted
HDF5 file (*seedname*.out.h5).
This program reads ``[model]``, ``[system]``, ``[impurity-solver]`` and ``[tool]`` block.

::

   $ dcore_post input-file

Input-file format
-----------------

.. literalinclude:: ../tutorial/bethe-t2g/dmft_bethe.ini

[model] block
~~~~~~~~~~~~~

dcore_pre, dcore_check and dcore_post read this block.

============= ============= ===================== ================================================================
Name          Type          Default               Description
============= ============= ===================== ================================================================
t             Float         1.0                   Transfer integral (Nearest neighbor)
t'            Float         0.0                   Transfer integral (Second nearest)
U             Float         0.0                   On-site Coulomb potential
J             Float         0.0                   On-site Hund potential
orbital_model String        single                Chosen from "single", "eg", "t2g", "full-d"
ncor          Integer       1                     Number of correlation shell (Only wannier90).
lattice       String        chain                 Chosen from "chain", "square", "cubic", "bethe", and "wannier90"
nelec         Float         1.0                   Number of electrons per unit cell.
seedname      String        dcore                 Name of the system.
                                                  It should be the same as the seedname of wannier90.
cshell        Integer array [(0,1),...]           Anguler momentum, and the number of orbitals of each
                                                  correlation shell (Only wannier90).
bvec          Float array   [(1.0,0.0,0.0), (0.0, Reciplocal lattice vectors
                            1.0,0.0),(0.0,0.0,1.0
                            )]
============= ============= ===================== ================================================================

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

   {\hat H} = \sum_{i j} \sum_{\alpha \beta}^{N_{\rm band}} \sum_{\sigma=\uparrow, \downarrow}
   t_{i \alpha j \beta} c_{i \alpha \sigma}^\dagger c_{j \beta \sigma}
   +h.c. + {\hat H}_{\rm int},

where :math:`{\hat H}_{\rm int}` is the Kanamori interaction term
given by

.. math::

    {\hat H}_{\rm int} = \sum_i \Big[
    \sum_{\alpha} U n_{i\alpha\uparrow} n_{i\alpha\downarrow}
    + \frac{1}{2} \sum_{\alpha \neq \beta, \sigma \sigma'}
    (U' - J\delta_{\sigma\sigma'}) n_{i\alpha\sigma} n_{i\beta\sigma'}
    - \sum_{\alpha \neq \beta} J(c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow} c_{i\beta\downarrow}^{\dagger} c_{i\beta\uparrow}
    + c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow}^{\dagger} c_{i\beta\uparrow} c_{i\beta\downarrow}
    \Big]

The parameter :math:`U'` is fixed at :math:`U'=U-2J`.

.. image:: model.png
   :width: 700
   :align: center

[system] block
~~~~~~~~~~~~~~

dcore_pre and dcore read this block.

======= =========== ===================== =================================================
Name    Type        Default               Description
======= =========== ===================== =================================================
beta    Float       1.0                   Inverse temperature.
n_iw    Integer     2048                  Number of Matsubara frequencies.
n_tau   Integer     10000                 Number of imaginary-time points.
dc_type Integer     -1                    Type of double-counting correction.
fix_mu  Bool        False                 Whether or not to use a fixed chemical potential.
nk      Integer     8                     Number of *k* along each line
nk0     Integer     0                     Number of *k* (Only wannier90)
nk1     Integer     0                     Number of *k* (Only wannier90)
nk2     Integer     0                     Number of *k* (Only wannier90)
prec_mu Float       0.0001                Threshold for calculating chemical potential
                                          with the bisection method.
======= =========== ===================== =================================================

[impurity_solver] block
~~~~~~~~~~~~~~~~~~~~~~~

dcore and dcore_post read this block.

==== ======= =========== ===================================================
Name Type    Default     Description
==== ======= =========== ===================================================
name String  TRIQS/cthyb Name of impurity solver. Choosen from "TRIQS/cthyb"
                         "TRIQS/hubbrad-I", and "ALPS/cthyb".
==== ======= =========== ===================================================

**... and other parameters (Solver dependent)**.
We have to specify additional parameters with types (e.g. ``n_cycles{int} = 500000``).
For more details, please see the reference page of
`TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/reference/solve_parameters.html#solve-parameters>`_,
`ALPS/cthyb <https://github.com/shinaoka/triqs_interface#program-parameters>`_, etc..

[control] block
~~~~~~~~~~~~~~~

dcore reads this block.

========= ======= ======= ==================================================
Name      Type    Default Description
========= ======= ======= ==================================================
max_step  Integer 100     Max number of SCF steps.
sigma_mix Float   0.5     Mixing parameter for self-energy.
delta_mix Float   0.5     Mixing parameter for hybridization function.
restart   Bool    False   Whether or not restart from a previous calculation
========= ======= ======= ==================================================

[tool] block
~~~~~~~~~~~~

dcore_check and dcore_post reads this block.

============= ================= ===================== ======================================================
Name          Type              Default               Description
============= ================= ===================== ======================================================
nnode         Integer           2                     Number of node for the *k* path
nk_line       Integer           8                     Number of *k* along each line
knode         Sting Float array [(G,0.0,0.0,0.0),     The name and the fractional coordinate of each k-node.
                                 (X,1.0,0.0,0.0)]
omega_min     Float             -1                    Minimum value of real frequency
omega_max     Float             1                     Max value of real frequency
Nomega        Integer           100                   Number of real frequencies
broadening    Float             0.1                   An additional Lorentzian broadening
eta           Float             0.0                   Imaginary frequency shift for the Pade approximation
n_pade        Integer           100                   Number of frequencies for the Pade approximation
============= ================= ===================== ======================================================
