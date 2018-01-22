.. _respacksrvo3:

RESPACK (SrVO\ :sub:`3`)
========================

.. todo::

   The QE-SrVO\ :sub:`3` tutorial is not complete.

Crystal structure of SrVO\ :sub:`3` (drawn by `VESTA <http://jp-minerals.org/vesta/en/>`_).

.. image:: qe/struct_srvo3.png
   :width: 200
   :align: center

SCF calculation of Quantum ESPRESSO
-----------------------------------

:download:`scf_srvo3.in <qe/scf_srvo3.in>`

.. literalinclude:: qe/scf_srvo3.in

The pseudopotentials are downloaded from
`Sr.pbe-spn-kjpaw_psl.0.2.3.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/Sr.pbe-spn-kjpaw_psl.0.2.3.upf>`_,
`V.pbe-spn-kjpaw_psl.0.2.3.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/V.pbe-spn-kjpaw_psl.0.2.3.upf>`_, and
`O.pbe-n-kjpaw_psl.0.1.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/O.pbe-n-kjpaw_psl.0.1.upf>`_.

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in scf_srvo3.in 

Wannierization
--------------   
   
Perform non-scf calculation for generating Bloch orbitals that are used
in the wannierization.

:download:`nscf_respack.in <respack/nscf_respack.in>`

.. literalinclude:: respack/nscf_respack.in

Then, run ``pw.x`` as

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in nscf_respack.in 

and convert outputs so that RESPACK can read them.
   
.. code-block:: bash
                
   $ qe2respack.sh srvo3.save/

This program locates in the ``util/qe2respack/`` directory of RESPACK.
   
Execute ``calc_wannier`` (in ``src/wannier/`` of RESPACK) for the actual wannierization.
with the following input file:

:download:`nscf_respack.in <respack/respack.in>`

.. literalinclude:: respack/respack.in

.. code-block:: bash
                
   $ calc_wannier < respack.in

(Optional) Check wannierization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

If we want to check the quality of the wannierization,
we plot the original and the wannier-interpolated band structure simultaneously.

First, we compute the band structure with the following input file:

:download:`band_srvo3.in <qe/band_srvo3.in>`

.. literalinclude:: qe/band_srvo3.in

.. code-block:: bash

   $ mpiexec -np 4 pw.x -in band_srvo3.in

:download:`bands_srvo3.in <qe/bands_srvo3.in>`

.. literalinclude:: qe/bands_srvo3.in

.. code-block:: bash

   $ mpiexec -np 4 bands.x -in bands_srvo3.in

.. code-block:: gnuplot

   plot [][11:18] "bands.out.gnu" u 1:2 w p tit "Orig", 12.3116 tit "E_F", "srvo3_band.dat" u ($1*0.6146):2 tit "Wannier" w l

.. image:: qe/band_srvo3.png
   :width: 500
   :align: center

DMFT calculation
----------------   
   
:download:`srvo3.ini <qe/srvo3.ini>`

.. literalinclude:: qe/srvo3.ini
                              
DMFT setup: dcore_pre
~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_pre srvo3.ini

Running self-consistent DFT+DMFT : dcore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore srvo3.ini

Post-processing and data analysis: dcore_post
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_post srvo3.ini
