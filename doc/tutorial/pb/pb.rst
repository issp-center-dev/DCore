.. _pbtutorial:

Pb with spin-orbit interaction
==============================

In this tutorial, we perform the DMFT calculation of Pb with the spin-orbit interaction.

SCF calculation of Quantum ESPRESSO
-----------------------------------

:download:`scf_pb.in <scf_pb.in>`

.. literalinclude:: scf_pb.in

The pseudopotential is included in 
`rel-pbe.0.3.1.tgz <http://theossrv1.epfl.ch/uploads/Main/NoBackup/rel-pbe.0.3.1.tgz>`_.

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in scf_pb.in 

Wannierization
--------------   
   
Generate Bloch orbitals for the Wannier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Perform non-scf calculation for generating Bloch orbitals that are used
in the wannierization.

:download:`nscf_pb.in <nscf_pb.in>`

.. literalinclude:: nscf_pb.in

This *k*\ -grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4

where ``kmesh.pl`` is located in the `utility/` directory of Wannier90.

Then, run ``pw.x`` as

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in nscf_pb.in 

Pre-process for Wannier90
~~~~~~~~~~~~~~~~~~~~~~~~~

Pre-process with wannier90 program.
It is always required before pw2wannier.x runs.
We will wannierize Pb 6p orbitals.

:download:`pb.win <pb.win>`

.. literalinclude:: pb.win

.. note::

   For the following DMFT calculation,
   the projection of the Wannier function should be the correct order, namely,

   ::

      begin projections
      First_shell(u)
      First_shell(d)
      Second_shell(u)
      Second_shell(d)
      :
      end projections

The *k* grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4 wan

.. code-block:: bash
                
   $ wannier90.x -pp pb 

QE to wannier90 interface
~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`pw2wan_pb.in <pw2wan_pb.in>`

.. literalinclude:: pw2wan_pb.in

.. code-block:: bash
                
   $ mpirun -np 4 pw2wan.x -in pw2wan_pb.in 

Wannier90
~~~~~~~~~   

Execute ``wannier90.x`` for the actual wannierization.
The input file is the same as that for the pre-processing run.

.. code-block:: bash
                
   $ wannier90 pb

DMFT calculation
----------------   
   
:download:`pb.ini <pb.ini>`

.. literalinclude:: pb.ini
   :language: ini
                              
DMFT setup: dcore_pre
~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_pre pb.ini

Running self-consistent DFT+DMFT : dcore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore pb.ini

Post-processing and data analysis: dcore_post
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block :: bash

   $ dcore_post pb.ini
   $ gnuplot pb_akw.gp

.. image:: akw_pb.png
   :width: 500
   :align: center

"+" indicates the original band structure.
