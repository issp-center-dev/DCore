.. Frequently-Asked Questions
.. ==========================


Impurity solvers
----------------

What is the difference between ALPSCore/CT-HYB and TRIQS/cthyb?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ALPS/CT-HYB`` is developed by one of the main developers of DCore, H. Shinaoka.
Both of ``ALPS/CT-HYB`` and ``TRIQS/cthyb`` implement the hybridization-expansion continuous-time quantum Monte Carlo method.
The main difference is the reliability of measurement of the single-particle Green's function.
ALPSCore/CT-HYB uses a more elaborate algorithm (worm sampling).
The non-worm conventional sampling, which is implemented in ``TRIQS/cthyb``,
may give wrong results in some situations (e.g. SOI coupling with orbital-diagonal bath).




..
   ``dcore`` crashes abnormally when using cthyb
   ---------------------------------------------

   Please retry.
