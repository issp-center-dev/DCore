pydmft
======

Perform DMFT.

Usage
-----

::

   $ pydmft input.ini seedname

Input-file format
-----------------

* ``system``
* ``beta``, float, 1.0, Inverse temperature.
* ``n_iw``, int, 2048, Number of Matsubara frequencies.
* ``n_tau``, int, 10000, Number of imaginary-time points.
* ``dc_type``, int, -1, Type of double-counting correction.
* ``fix_mu``, bool, False, Whether or not to use a fixed chemical potential.
* ``impurity_solver``
* ``N_l``, int, 50, Number of Legendre polynomials.
* ``name``, str, 'TRIQS/cthyb', Name of impurity solver.
* ``control``
* ``max_step``, int, 100, Max number of SCF steps.
* ``sigma_mix``, float, 0.5, Mixing parameter for self-energy.
* ``delta_mix``, float, 0.5, Mixing parameter for hybridization function.

