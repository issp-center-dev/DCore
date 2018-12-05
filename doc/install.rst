
.. highlight:: bash

.. _installation:

Installation
============

Download
--------

You can download the source files in two ways.

- **[Release version]**

  Go to the `release page <https://github.com/issp-center-dev/DCore/releases>`_ and download the latest tar or zip file.
  You can unpack, for example, the tar file by

  .. code-block:: bash

      $ tar zxvf v1.x.x.tar.gz

- **[Develop version]**

  Newly implemented features are available in the source code in GitHub repository. You can clone the repository by

  .. code-block:: bash

      $ git clone https://github.com/issp-center-dev/DCore.git dcore.src

  The master branch basically corresponds to the latest released package, and the develop branch includes new features. Note that the develop branch may not be well tested and contains some bugs.

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

      If you want to treat the spin-orbit coupling in ``TRIQS/cthyb`` solver,
      it must be built with the following CMake options:

      ::

         -DHYBRIDISATION_IS_COMPLEX=ON -DLOCAL_HAMILTONIAN_IS_COMPLEX=ON

   .. note::

      One must build ALPSCore, TRIQS, TRIQS/DFTTools, and other solvers with the same C++ compiler and the same C++ standard (C++14).

   .. note::

      ALPSCore/CT-HYB is a standalone program with an interface based on text and HDF5 files.
      "triqs_interface" is a TRIQS-compatible Python interface of ``ALPSCore/CT-HYB``.
      This allows to use this impurity solver transparently from DCore (built with TRIQS python architecture).
      The main developer of ``ALPSCore/CT-HYB``, H. Shinaoka, is one of the developers of DCore.

      Both of ``ALPSCore/CT-HYB`` and ``TRIQS/cthyb`` implement the hybridization-expansion continuous-time quantum Monte Carlo method.
      The main difference is the reliability of measurement of the single-particle Green's function.
      ALPSCore/CT-HYB uses a more elaborate algorithm (worm sampling).
      The non-worm conventional sampling, which is implemented in ``TRIQS/cthyb``,
      may give wrong results in some situations (e.g. SOI coupling with orbital-diagonal bath).

Installation steps
------------------

#. Create an empty directory where you will compile the code

   ::

     $ mkdir dcore.build && cd dcore.build

#. In the build directory, call cmake command with an option to specify the path to TRIQS library

   ::

     $ cmake -DTRIQS_PATH=path_to_triqs path_to_dcore_src

   Here, *path_to_triqs* refers to your ``TRIQS`` install directory, and *path_to_dcore_src* refers to the ``DCore`` source directory.
   If you want to install DCore to somewhere other than the triqs directory, please set CMAKE_INSTALL_PREFIX to your install directory.
   If the cmake command succeeded, you will see the following message.

   ::

     -- Build files have been written to: /.../dcore.build

#. Compile the code by

   ::

     $ make

#. We recommend that you run the tests to check if the compiling is properly finished. Type

   ::

     $ make test

   The test results look like

   ::

     Running tests...
     Test project /Users/k-yoshimi/CLionProjects/DCore/build
           Start  1: typed_parser
      1/21 Test  #1: typed_parser .....................   Passed    0.62 sec
           Start  2: openmx
      2/21 Test  #2: openmx ...........................   Passed    0.66 sec
           Start  3: respack
      3/21 Test  #3: respack ..........................   Passed    0.64 sec
           Start  4: pre_preset
      4/21 Test  #4: pre_preset .......................   Passed    0.91 sec
           Start  5: pre_wannier
      5/21 Test  #5: pre_wannier ......................   Passed    2.42 sec
           Start  6: pre_wannier_so
      6/21 Test  #6: pre_wannier_so ...................   Passed    0.77 sec
           Start  7: pre_respack
      7/21 Test  #7: pre_respack ......................   Passed    0.86 sec
           Start  8: pre_respack_so
      8/21 Test  #8: pre_respack_so ...................   Passed    0.75 sec
           Start  9: main_chain_hubbardI
      9/21 Test  #9: main_chain_hubbardI ..............   Passed    2.92 sec
           Start 10: main_chain_hubbardI_so
     10/21 Test #10: main_chain_hubbardI_so ...........   Passed    3.51 sec
           Start 11: main_chain_triqsqmc
     11/21 Test #11: main_chain_triqsqmc ..............***Failed    0.70 sec
           Start 12: main_chain_triqsqmc_so
     12/21 Test #12: main_chain_triqsqmc_so ...........***Failed    0.68 sec
           Start 13: check_chain_hubbardI
     13/21 Test #13: check_chain_hubbardI .............   Passed    1.49 sec
           Start 14: check_chain_hubbardI_so
     14/21 Test #14: check_chain_hubbardI_so ..........   Passed    1.48 sec
           Start 15: check_chain_triqsqmc
     15/21 Test #15: check_chain_triqsqmc .............***Failed    0.74 sec
           Start 16: check_chain_triqsqmc_so
     16/21 Test #16: check_chain_triqsqmc_so ..........***Failed    0.76 sec
           Start 17: check_chain_alpsqmc
     17/21 Test #17: check_chain_alpsqmc ..............***Failed    0.73 sec
           Start 18: post_chain_hubbardI
     18/21 Test #18: post_chain_hubbardI ..............   Passed    2.28 sec
           Start 19: post_chain_hubbardI_so
     19/21 Test #19: post_chain_hubbardI_so ...........   Passed    2.18 sec
           Start 20: post_chain_triqsqmc
     20/21 Test #20: post_chain_triqsqmc ..............***Failed    0.71 sec
           Start 21: post_chain_triqsqmc_so
     21/21 Test #21: post_chain_triqsqmc_so ...........***Failed    0.72 sec

     67% tests passed, 7 tests failed out of 21

     Total Test time (real) =  26.56 sec

     The following tests FAILED:
     	 11 - main_chain_triqsqmc (Failed)
     	 12 - main_chain_triqsqmc_so (Failed)
     	 15 - check_chain_triqsqmc (Failed)
     	 16 - check_chain_triqsqmc_so (Failed)
     	 17 - check_chain_alpsqmc (Failed)
     	 20 - post_chain_triqsqmc (Failed)
     	 21 - post_chain_triqsqmc_so (Failed)

   In the above example, all tests related to QMC are failed, because only the Hubbard-I solver is installed. The tests from 1 to 8 should be passed regardless of impurity solvers you installed.

#. Finally, install by

   ::

     $ make install

   ``DCore`` is installed in the directory *path_to_triqs*.

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
