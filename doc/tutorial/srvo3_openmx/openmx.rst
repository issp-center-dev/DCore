OpenMX
======

.. include:: ../../warning_compatibility.rst

OpenMX is a first-principles program package based on
the numerical localized basis set and pseudopotentials.
OpenMX itself can generate hopping parameter in the wannier90 format.
In this tutorial, we demonstrate the calculation of SrVO\ :sub:`3`.

.. note::

   This tutorial requires large computational resources or long simulation time.

SCF computation and Wannier with OpenMX
---------------------------------------

:download:`scf.dat <scf.dat>`

.. literalinclude:: scf.dat

.. code-block:: bash

   $ openmx scf.dat

Then, convert the OpenMX output to the wannier90 format.
It can be performed with ``openmx2dcore`` utility as:

.. code-block:: bash

   $ openmx2dcore.py SrVO3 srvo3

DMFT calculation
----------------

:download:`srvo3.ini <srvo3_openmx.ini>`

.. literalinclude:: srvo3_openmx.ini
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

   $ dcore srvo3.ini

Post-processing and data analysis: dcore_anacont and dcore_spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_anacont srvo3.ini
   $ dcore_spectrum srvo3.ini
   $ cd post
   $ sed -e "s/every 10/every 1/g" akw.gp
   $ gnuplot akw.gp

.. image:: akw_srvo3.png
   :width: 500
   :align: center

"+" indicates the original band structure.
