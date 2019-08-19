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

Please see :ref:`howtocthyb` for the details of the parameter setting.

.. note::

   The parameter ``n_cycles{int}`` should be tuned in inverse proportion to the number of MPI processes.
   The following result is obtained with 432 MPI processes at ``n_cycles{int} = 10000``
   (70 seconds per DMFT cycle on ISSP system B).
   If we want to compute by using 32 MPI processes at the same accuracy,
   ``n_cycles{int}`` should be 10000\*432/32=135000.

DMFT setup: dcore_pre
~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_pre srvo3.ini

Running self-consistent DMFT calculation: dcore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore srvo3.ini --np 32

Post-processing and data analysis: dcore_post
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_post srvo3.ini --np 32
   $ gnuplot nis_akw.gp

.. image:: akw_srvo3.png
   :width: 500
   :align: center

"+" indicates the original band structure.
