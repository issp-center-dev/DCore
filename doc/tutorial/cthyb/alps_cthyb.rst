.. _howtoalpscthyb:

How to choose parameters for ALPS/cthyb
========================================


Basic parameters for QMC
------------------------

DCore invokes the `ALPSCore/CT-HYB solver <https://github.com/ALPSCore/CT-HYB>`_ developed by H. Shinaoka, E. Gull and P. Werner
via a `triqs-compatible Python interface <https://github.com/shinaoka/triqs_interface>`_.
In the following, we give a short description on some of the parameters (see `the official document <https://github.com/shinaoka/triqs_interface/blob/master/README.md>`_ for more details).

In DCore calculations, you may specify the following two parameters.

 * ``max_time{int}``

   **This parameter affects the statistical error and run time.**
   **In most cases, this is the only one parameter you have to specify.**
 
   Total simulation time (in units of second) including thermalization process.
   Please run longer to reduce statistical noise.

 
 * ``thermalization_time{int}``

   **This parameter affects the systematical error.**
 
   Length of the thermalization steps in units of second.
   The default value is 10 % of max_time, which is a reasonable choice in most cases.
   Thus, we recommend to use the default value unless the solver complains about thermalization time.
