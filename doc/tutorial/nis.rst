NiS
===

Crystal structure of NiS (drawn by `VESTA <http://jp-minerals.org/vesta/en/>`_).

.. image:: qe/struct_nis.png
   :width: 200
   :align: center

SCF calculation of Quantum ESPRESSO
-----------------------------------

:download:`scf_nis.in <qe/scf_nis.in>`

.. literalinclude:: qe/scf_nis.in

The pseudopotentials are downloaded from
`Ni.pbe-n-kjpaw_psl.0.1.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/Ni.pbe-n-kjpaw_psl.0.1.upf>`_ and
`S.pbe-n-kjpaw_psl.0.1.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/S.pbe-n-kjpaw_psl.0.1.upf>`_.

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in scf_nis.in 

Wannierization
--------------   
   
Generate Bloch orbitals for the Wannier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Perform non-scf calculation for generating Bloch orbitals that are used
in the wannierization.

:download:`nscf_nis.in <qe/nscf_nis.in>`

.. literalinclude:: qe/nscf_nis.in

This *k*\ -grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4

where ``kmesh.pl`` is located in the `utility/` directory of Wannier90.

Then, run ``pw.x`` as

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in nscf_nis.in 

Preprocess for Wannier90
~~~~~~~~~~~~~~~~~~~~~~~~

Pre-process with wannier90 program.
It is always required before pw2wannier.x runs.

:download:`srvo3.win <qe/srvo3.win>`

.. literalinclude:: qe/srvo3.win

This *k* grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4 wan

.. code-block:: bash
                
   $ wannier90.x -pp srvo3 

QE to wannier90 interface
~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`pw2wan_nis.in <qe/pw2wan_nis.in>`

.. literalinclude:: qe/pw2wan_nis.in

.. code-block:: bash
                
   $ mpirun -np 4 pw2wan.x -in pw2wan_nis.in 

Wannier90
~~~~~~~~~   

Execute ``wannier90.x`` for the actual wannierization.
The input file is the same as that for the preprocessing run.

.. code-block:: bash
                
   $ wannier90 srvo3 

(Optional) Check wannierization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

If we want to check the quarity of the wannierization,
we plot the original and the wannier-interpolated band structure simalteneously.

First, we compute the band structure with the following input file:

:download:`band_nis.in <qe/band_nis.in>`

.. literalinclude:: qe/band_nis.in

.. code-block:: bash

   $ mpiexec -np 4 pw.x -in band_nis.in

:download:`bands_nis.in <qe/bands_nis.in>`

.. literalinclude:: qe/bands_nis.in

.. code-block:: bash

   $ mpiexec -np 4 bands.x -in bands_nis.in

.. code-block:: gnuplot

   plot [][3:13] "bands.out.gnu" u 1:2 w p tit "Orig", 10.913 tit "E_F", "nis_band.dat" u ($1*0.549):2 tit "Wannier" w l

.. image:: qe/band_nis.png
   :width: 500
   :align: center

DMFT calculation
----------------   
   
:download:`nis.ini <qe/nis.ini>`

.. literalinclude:: qe/nis.ini
                              
DMFT setup: pydmft_pre
~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ pydmf_pre nis.ini

Running self-consistent DFT+DMFT : pydmft
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ pydmf nis.ini

Post-processing and data analysis: pydmft_post
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ pydmf_post nis.ini
