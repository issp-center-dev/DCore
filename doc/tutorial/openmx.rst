Tutorial with OpenMX
====================

OpenMX is a first-principles program package based on
the numerical localized basis set and pseudopotentials.
OpenMX itself can generate hopping parameter in the wannier90 format.
In this tutorial, we demonstrate the calculation of SrVO\ :sub:`3`.

SCF computation and Wannier with OpenMX
---------------------------------------

:download:`scf.dat <openmx/scf.dat>`

.. literalinclude:: openmx/scf.dat

DMFT
----

First, convert the OpenMX output to the wannier90 format.
It can be performed with ``openmx2dcore`` utility as:

.. code-block:: bash
                
   $ openmx2dcore SrVO3.HWR srvo3

Then ``srvo3_hr.dat`` is generated.
The rest is the same as :ref:`example for QE <qesrvo3>`.   
