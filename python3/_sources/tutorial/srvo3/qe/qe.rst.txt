.. _qesrvo3:

Wanniernization using Quantum ESPRESSO
=======================================

SCF calculation of Quantum ESPRESSO
-----------------------------------

:download:`scf_srvo3.in <scf_srvo3.in>`

.. literalinclude:: scf_srvo3.in

The pseudopotentials are downloaded from
`Sr.pbe-spn-kjpaw_psl.0.2.3.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/Sr.pbe-spn-kjpaw_psl.0.2.3.upf>`_,
`V.pbe-spn-kjpaw_psl.0.2.3.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/V.pbe-spn-kjpaw_psl.0.2.3.upf>`_, and
`O.pbe-n-kjpaw_psl.0.1.upf <http://theossrv1.epfl.ch/uploads/Main/NoBackup/O.pbe-n-kjpaw_psl.0.1.upf>`_.

.. code-block:: bash

   $ mpirun -np 4 pw.x -in scf_srvo3.in

Wannierization
--------------

Generate Bloch orbitals for the Wannier
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Perform non-scf calculation for generating Bloch orbitals that are used
in the wannierization.

:download:`nscf_srvo3.in <nscf_srvo3.in>`

.. literalinclude:: nscf_srvo3.in

This *k*\ -grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4

where ``kmesh.pl`` is located in the `utility/` directory of Wannier90.

Then, run ``pw.x`` as

.. code-block:: bash

   $ mpirun -np 4 pw.x -in nscf_srvo3.in

Pre-process for Wannier90
~~~~~~~~~~~~~~~~~~~~~~~~~

Pre-process with wannier90 program.
It is always required before pw2wannier.x runs.

:download:`srvo3.win <srvo3.win>`

.. literalinclude:: srvo3.win

This *k* grid is generated as follows:

.. code-block:: bash

   $ kmesh.pl 4 4 4 wan

.. code-block:: bash

   $ wannier90.x -pp srvo3

QE to wannier90 interface
~~~~~~~~~~~~~~~~~~~~~~~~~

:download:`pw2wan_srvo3.in <pw2wan_srvo3.in>`

.. literalinclude:: pw2wan_srvo3.in

.. code-block:: bash

   $ mpirun -np 4 pw2wan.x -in pw2wan_srvo3.in

Wannier90
~~~~~~~~~

Execute ``wannier90.x`` for the actual wannierization.
The input file is the same as that for the pre-processing run.

.. code-block:: bash

   $ wannier90.x srvo3

(Optional) Check wannierization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to check the quality of the wannierization,
you can plot the original and the wannier-interpolated band structure simultaneously.

First, compute the band structure with the following input file:

:download:`band_srvo3.in <band_srvo3.in>`

.. literalinclude:: band_srvo3.in

.. code-block:: bash

   $ mpiexec -np 4 pw.x -in band_srvo3.in

:download:`bands_srvo3.in <bands_srvo3.in>`

.. literalinclude:: bands_srvo3.in

.. code-block:: bash

   $ mpiexec -np 4 bands.x -in bands_srvo3.in

.. .. code-block:: gnuplot
.. code-block:: guess

   plot [][11:18] "bands.out.gnu" u 1:2 w p tit "Orig", 12.3116 tit "E_F", "srvo3_band.dat" u ($1*0.6146):2 tit "Wannier" w l

.. image:: band_srvo3.png
   :width: 500
   :align: center
