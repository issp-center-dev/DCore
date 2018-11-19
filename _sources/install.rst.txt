
.. highlight:: bash

.. _installation:

Installation
============

Prerequisites
-------------

#. **[mandatory]** ``DCore`` is implemented using :ref:`TRIQS <triqslibs:welcome>` library.
   You first need to install ``TRIQS`` library prior to installing all other programs.
   In the following, we suppose that ``TRIQS`` is installed in directory *path_to_triqs*.

   .. note::

     ``DCore`` supports ``TRIQS`` version **1.4.2** (the current version 2.0 is not supported).
     Since the installation manual in :ref:`TRIQS <triqslibs:welcome>` is not compatible with **1.4.2**, we briefly describe below how to install it.

     #. You can download the zip file or tar file from https://github.com/TRIQS/triqs/releases/tag/1.4.2.
        The version **1.4.2** is not compatible with the latest version of h5py (>=2.8.0).
        If you encounter any problem, you can clone the repository and checkout **1.4.x** branch, which contains bug fixes.
        We suppose that source files are located in *path_to_triqs_src* directory.

     #. In an empty directory, type the following command:

        .. code-block:: bash

          $ cmake -DCMAKE_INSTALL_PREFIX=path_to_triqs path_to_triqs_src
          $ make
          $ make test
          $ make install

#. **[mandatory]**
   You also need :ref:`DFTTools <triqsdfttools:dft>`, which runs on the ``TRIQS`` library.

   .. note::

     The current version in the GitHub repository is not compatible with ``TRIQS`` version 1.4.2.
     You need to get old code that is compatible with 1.4.2.

     #. Clone the repository, and checkout commit **d005756**.

        .. code-block:: bash

          $ git clone git@github.com:TRIQS/dft_tools.git
          $ cd dft_tools
          $ git checkout d005756

     #. To build the source files, make an empty directory, move into it, and then type the following command:

        .. code-block:: bash

          $ cmake -DCMAKE_INSTALL_PREFIX=path_to_triqs\
            -DTRIQS_PATH=path_to_triqs path_to_dft_tools
          $ make
          $ make test
          $ make install

#. **[optional]** You will also need at least one impurity solver.
   At present, ``DCore`` supports the following programs:

   - `Hubbard-I solver <https://triqs.ipht.cnrs.fr/1.x/applications/hubbardI/>`_

   - `TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/index.html>`_
     [Note: checkout tag **1.4.2** on master branch, or download the released package from `here <https://github.com/TRIQS/cthyb/releases/tag/1.4.2>`_]

   - `ALPSCore/CT-HYB <https://github.com/ALPSCore/CT-HYB>`_ + `triqs_interface <https://github.com/shinaoka/triqs_interface>`_

   .. note::

      If you want to treat the spin-orbit coupling in `TRIQS/cthyb` solver,
      it must be built with the following CMake options:

      ::

         -DHYBRIDISATION_IS_COMPLEX=ON -DLOCAL_HAMILTONIAN_IS_COMPLEX=ON

   .. note::

      One must build ALPSCore, TRIQS, TRIQS/DFTTools, and other solvers with the same C++ compiler and the same C++ standard (C++14).

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

.. Version compatibility
.. ---------------------
..
.. The current version of DCore supports TRIQS 1.4, and ALPSCore 2.1 or later.
.. Be careful that the version of the TRIQS library and of the dft tools must be
.. compatible (more information on the :ref:`TRIQS website <triqslibs:welcome>`).
.. If you want to use a version of the dft tools that is not the latest one, go
.. into the directory with the sources and look at all available versions::
..
..      $ cd src && git tag
..
.. Checkout the version of the code that you want, for instance::
..
..      $ git co 1.4
..
.. Then follow the steps 2 to 5 described above to compile the code.
