.. _DFTDMFTtutorial:

Tutorial with QuantumESPRESSO and wannier90
===========================================

Quantum ESPRESSO is a first-principles program package based on
the plane wave and pseudopotentials.

wannier90 is a library and post-processing tool for generating wannier orbitals
from results by many DFT programs.

In this tutorial we will perform DFT+DMFT calculations of SrVO\ :sub:`3`.

SCF calculation of Quantum ESPRESSO
-----------------------------------

:download:`scf.in <qe/scf.in>`

.. literalinclude:: qe/scf.in

The pseudopotentials are downloaded from
`Sr.pbe-spn-kjpaw_psl.0.2.3.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/Sr.pbe-spn-kjpaw_psl.0.2.3.upf>`_,
`V.pbe-spn-kjpaw_psl.0.2.3.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/V.pbe-spn-kjpaw_psl.0.2.3.upf>`_, and
`O.pbe-n-kjpaw_psl.0.1.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/O.pbe-n-kjpaw_psl.0.1.upf>`_.

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in scf.in 

Wannierization
--------------   
   
Generate Bloch orbitals for the Wannier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`nscf.in <qe/nscf.in>`

.. literalinclude:: qe/nscf.in

This *k* grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4

where ``kmesh.pl`` is located in the `utility/` directory of Wannier90.

Run ``pw.x`` as

.. code-block:: bash
                
   $ mpirun -np 4 pw.x -in nscf.in 

Preprocess for Wannier90
~~~~~~~~~~~~~~~~~~~~~~~~

:download:`srvo3.win <qe/srvo3.win>`

.. literalinclude:: qe/srvo3.win

This *k* grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4 wan

.. code-block:: bash
                
   $ wannier90 -pp srvo3 

QE to wannier90 interface
~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`pw2wan.in <qe/pw2wan.in>`

.. literalinclude:: qe/pw2wan.in

.. code-block:: bash
                
   $ mpirun -np 4 pw2wan.x -in pw2wan.in 

Wannier90
~~~~~~~~~   

.. code-block:: bash
                
   $ wannier90 srvo3 

DMFT calculation
----------------   
   
DMFT setup: pydmft_pre
~~~~~~~~~~~~~~~~~~~~~~

:download:`model.in <qe/model.in>`

.. literalinclude:: qe/model.in
                              
.. code-block :: bash

   $ pydmf_pre model.in

Running self-consistent DFT+DMFT : pydmft
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`input.ini <qe/input.ini>`

.. literalinclude:: qe/input.ini
                              
.. code-block :: bash

   $ pydmf input.ini srvo3

Post-processing and data analysis: pydmft_post
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`input.ini <qe/input.ini>`

.. literalinclude:: qe/input.ini
                              
.. code-block :: bash

   $ pydmf_post input.ini srvo3
