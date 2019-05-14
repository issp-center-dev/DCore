.. _pbtutorial:

Spin-orbit interaction: The case with Pb
========================================

.. include:: ../../warning_compatibility.rst

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
   the projections of the Wannier functions should be given in the correct order, namely,

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
                
   $ wannier90.x pb

DMFT calculation
----------------   
   
:download:`pb.ini <pb.ini>`

.. literalinclude:: pb.ini
   :language: ini


.. note::
   
   In this tutorial, ``Hubbard-I`` solver is selected to calculate fast.
   However, in actual calculation, ``QMC`` method is more suitable for metallic state.
   Thus, we recommend that you run the DMFT calculation and check the results by using ``QMC`` solver (``TRIQS/cthyb`` or ``ALPS/cthyb``).
   In the following, modified parts of ``pb.ini`` to use the ``QMC`` solver are shown.

   - TRIQS/cthyb
     
      .. code-block:: ini
      
         [impurity_solver]
         name = TRIQS/cthyb
         n_cycles{int} = 100000
         n_warmup_cycles{int} = 10000
         length_cycle{int} = 500
         move_double{bool} = True
         verbosity{int} = 10

   - ALPS/cthyb
     
      .. code-block:: ini
      
         [impurity_solver]
         name = ALPS/cthyb
         max_time{int} = 1200
         thermalization_time{int} = 120
         basis_rotation{int} = 1
         verbosity{int} = 1

   - Common change in ``[tool]`` box
   
      .. code-block:: ini
      
         [tool]
         omega_pade = 10.0

   You also notice that simple interaction is selected in this tutorial.
   When you want to calculate by using more realistic interactions,
   please try to use **RESPACK** to obtain the interactions.
      
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
