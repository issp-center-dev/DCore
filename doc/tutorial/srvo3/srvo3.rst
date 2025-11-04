.. _srvo3:

SrVO\ :sub:`3`
==============

.. note::

   This tutorial requires parallel computing.

Crystal structure of SrVO\ :sub:`3` (drawn by `VESTA <http://jp-minerals.org/vesta/en/>`_).

.. image:: struct_srvo3.png
   :width: 200
   :align: center

Construction of Wannier functions
-----------------------------------

Maximally localized Wannier functions for the *t*\ :sub:`2g` manifold can be constructed by using DFT code.
Please download precomputed data for Wannier functions

:download:`srvo3_hr.dat <qe/srvo3_hr.dat>`

and save it to your working directly. This data was computed by using Quantum ESPRESSO and Wannier90.
The procedure of Wanniernization is detailed in :doc:`qe/qe` and :doc:`openmx/openmx`.

DMFT calculation
----------------

:download:`srvo3.ini <srvo3.ini>`

.. literalinclude:: srvo3.ini
   :language: ini

To use TRIQS/cthyb, please use 

.. literalinclude:: srvo3_triqs_cthyb.ini
   :language: ini

instead.
Please see :ref:`howtocthyb` for the details of the parameter setting for TRIQS/cthyb. 

To generate reference data, we used 48 MPI processes (ISSP system B).
The total computational time for 12 iterations is around 30 mins.

The accuracy of QMC results can be improved by setting a longer simulation time (ALPS/CT-HYB) or a larger n_cycles (TRIQS/cthyb).
For ALPS/CT-HYB, the parameter time_limit is given in seconds.

DMFT setup: dcore_pre
~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_pre srvo3.ini

Running self-consistent DMFT calculation: dcore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore srvo3.ini --np 48

Post-processing and data analysis: dcore_anacont and dcore_spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_anacont srvo3.ini
   $ dcore_spectrum srvo3.ini --np 48
   $ cd post
   $ gnuplot akw.gp

.. image:: akw_srvo3.png
   :width: 500
   :align: center

The left and right panels show results for up and down spins, respectively.
