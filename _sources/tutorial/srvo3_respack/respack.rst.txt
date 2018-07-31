.. _respacksrvo3:

Downfolding with RESPACK (SrVO\ :sub:`3`)
=========================================

.. note::

   This tutorial requires large computational resources or the long simulation time.

Crystal structure of SrVO\ :sub:`3` (drawn by `VESTA <http://jp-minerals.org/vesta/en/>`_).

.. image:: ../srvo3_qe/struct_srvo3.png
   :width: 200
   :align: center

SCF calculation of Quantum ESPRESSO
-----------------------------------

:download:`scf_srvo3_r.in <scf_srvo3_r.in>`

.. literalinclude:: scf_srvo3_r.in

The pseudopotentials are downloaded from

http://www.quantum-simulation.org/potentials/sg15_oncv/upf/O_ONCV_PBE-1.0.upf

http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Sr_ONCV_PBE-1.0.upf

http://www.quantum-simulation.org/potentials/sg15_oncv/upf/V_ONCV_PBE-1.0.upf

They are part of
`The SG15 Optimized Norm-Conserving Vanderbilt (ONCV) pseudopotentials <http://www.quantum-simulation.org/potentials/sg15_oncv/>`_.
For the downfolding with RESPACK, we should use the norm-conserving pseudopotentials rather than
the ultrasoft pseudopotentials or PAW.

The SCF calculation of the electronic charge is performed as follows:

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in scf_srvo3.in 

Wannierization
--------------   
   
Next, we perform non-scf calculation for generating Bloch orbitals that are used
in the wannierization.

:download:`nscf_respack.in <nscf_respack.in>`

.. literalinclude:: nscf_respack.in

Then, run ``pw.x`` as

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in nscf_respack.in 

and convert outputs so that RESPACK can read them.
   
.. code-block:: bash
                
   $ qe2respack.sh srvo3.save/

This program locates in the ``util/qe2respack/`` directory of RESPACK.
   
Execute ``calc_wannier`` (in ``src/wannier/`` of RESPACK) for the actual wannierization.
with the following input file:

:download:`nscf_respack.in <respack.in>`

.. literalinclude:: respack.in

.. code-block:: bash
                
   $ calc_wannier < respack.in

(Optional) Check wannierization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

If you want to check the quality of the wannierization,
you can plot the original and the wannier-interpolated band structure simultaneously.

First, compute the band structure with the following input file:

:download:`band_srvo3_r.in <band_srvo3_r.in>`

.. literalinclude:: band_srvo3_r.in

.. code-block:: bash

   $ mpiexec -np 4 pw.x -in band_srvo3_r.in

:download:`bands_srvo3.in <../srvo3_qe/bands_srvo3.in>`

.. literalinclude:: ../srvo3_qe/bands_srvo3.in

.. code-block:: bash

   $ mpiexec -np 4 bands.x -in bands_srvo3.in

.. code-block:: gnuplot

   plot [][11:18] "bands.out.gnu" u 1:2 w p tit "Orig", 12.3116 tit "E_F", "dir-wan/dat.iband" u ($1*2.5731):2 tit "Wannier" w l

.. image:: ../srvo3_qe/band_srvo3.png
   :width: 500
   :align: center

Dielectric matrix and Effective interaction
-------------------------------------------           

Next, we move on the calculation of the dielectric matrix with cRPA.
We use the program ``calc_chiqw`` (in ``src/calc_chiqw`` in RESPACK) as

.. code-block:: bash
                
   $ mpiexec calc_chiqw < respack.in

where the input file is the same as above.

After we compute the dielectric matrix, we calculate the effective interaction :math:`U` and :math:`J` as   

.. code-block:: bash
                
   $ calc_w3d < respack.in
   $ calc_j3d < respack.in

The output of these program should be transformed into the wannier90 format by using the utility program
``respack2wan90.py`` in ``bin/`` directory of triqs.   

.. code-block:: bash
                
   $ respack2wan90.py srvo3

The command-line argument (``"srvo3"`` in this case) must be the same as ``seedname`` in the DCore input.
   
DMFT calculation
----------------   
   
:download:`srvo3.ini <srvo3_respack.ini>`

.. literalinclude:: srvo3_respack.ini
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

   $ dcore_pre srvo3_respack.ini

Running the DMFT calculation: dcore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore srvo3_respack.ini

Post-processing and data analysis: dcore_post
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_post srvo3_respack.ini
   $ sed -e "s/every 10/every 3/g" srvo3_akw.gp
   $ gnuplot nis_akw.gp

.. image:: akw_srvo3.png
   :width: 500
   :align: center

"x" indicates the original band structure.
