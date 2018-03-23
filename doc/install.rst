
.. highlight:: bash

.. _installation:
               
Installation
============

Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` library and
   :ref:`DFTTools <triqsdfttools:dft>` must be installed in advance.
   In the following, we will suppose that they are installed in the ``path_to_triqs`` directory.

#. You will also need at least one impurity solver.
   At present, ``DCore`` supports the following programs:

   - `Hubbard-I solver <https://triqs.ipht.cnrs.fr/1.x/applications/hubbardI/>`_

   - `TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/index.html>`_

   - `ALPSCore/CT-HYB <https://github.com/ALPSCore/CT-HYB>`_ + `triqs_interface <https://github.com/shinaoka/triqs_interface>`_

   .. note::

      If you want to treat the spin-orbit coupling in `TRIQS/cthyb` solver,
      it must be built with the following CMake options:

      ::

         -DHYBRIDISATION_IS_COMPLEX=ON -DLOCAL_HAMILTONIAN_IS_COMPLEX=ON

   .. note::

      One must build ALPSCore, TRIQS, TRIQS/DFTTools, and other solvers with the same C++ compiler and the same C++ standard (C++14).

      ::

   .. note::

      ALPSCore/CT-HYB is a standalone program with an interface based on text and HDF5 files.
      "triqs_interface" is a TRIQS-compatible Python interface of ALPSCore/CT-HYB.
      This allows to use this impurity solver transparently from DCore (built with TRIQS python architecture).
      The main developer of ALPSCore/CT-HYB, H. Shinaoka, is one of the developers of DCore.

      Both of ALPSCore/CT-HYB and TRIQS/cthyb implement the hybridization-expansion continuous-time quantum Monte Carlo method.
      The main difference is the reliability of measurement of the single-particle Green's function.
      ALPSCore/CT-HYB uses a more elaborate algorithm (worm sampling).
      The non-worm conventional sampling, which is implemented in TRIQS/cthyb,
      may give wrong results in some situations (e.g. SOI coupling with orbital-diagonal bath).

      ::

Installation steps 
------------------

#. Download the sources from github :: 
 
     $ git clone https://github.com/issp-center-dev/DCore.git dcore.src
 
#. Create an empty build directory where you will compile the code:: 
 
     $ mkdir dcore.build && cd dcore.build 
 
#. In the build directory call cmake specifying where the TRIQS library is installed:: 
 
     $ cmake -DTRIQS_PATH=path_to_triqs ../dcore.src 
 
#. Compile the code, run the tests and install the application:: 
 
     $ make 
     $ make test 
     $ make install 
 
Version compatibility 
--------------------- 
 
The current version of DCore supports TRIQS 1.4, and ALPSCore 2.1 or later.
Be careful that the version of the TRIQS library and of the dft tools must be
compatible (more information on the :ref:`TRIQS website <triqslibs:welcome>`).
If you want to use a version of the dft tools that is not the latest one, go
into the directory with the sources and look at all available versions:: 
 
     $ cd src && git tag 
 
Checkout the version of the code that you want, for instance:: 
 
     $ git co 1.4
 
Then follow the steps 2 to 5 described above to compile the code. 
