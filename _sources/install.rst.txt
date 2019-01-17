
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
   If you want to install ``DCore`` to somewhere other than the triqs directory, please add the option ``-DCMAKE_INSTALL_PREFIX=install_path``, where *install_path* is your install directory.
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
     Test project /home/otsuki/gitclones/dcore.build
           Start  1: typed_parser
      1/20 Test  #1: typed_parser .....................   Passed    0.27 sec
           Start  2: openmx
      2/20 Test  #2: openmx ...........................   Passed    0.25 sec
           Start  3: respack
      3/20 Test  #3: respack ..........................   Passed    0.20 sec
           Start  4: pre_preset
      4/20 Test  #4: pre_preset .......................   Passed    0.56 sec
           Start  5: pre_wannier
      5/20 Test  #5: pre_wannier ......................   Passed    1.63 sec
           Start  6: pre_wannier_so
      6/20 Test  #6: pre_wannier_so ...................   Passed    0.27 sec
           Start  7: pre_respack
      7/20 Test  #7: pre_respack ......................   Passed    0.25 sec
           Start  8: pre_respack_so
      8/20 Test  #8: pre_respack_so ...................   Passed    0.26 sec
           Start  9: main_chain_hubbardI
      9/20 Test  #9: main_chain_hubbardI ..............   Passed    1.66 sec
           Start 10: main_chain_hubbardI_so
     10/20 Test #10: main_chain_hubbardI_so ...........   Passed    1.93 sec
           Start 11: main_chain_triqsqmc
     11/20 Test #11: main_chain_triqsqmc ..............   Passed   19.04 sec
           Start 12: main_chain_triqsqmc_so
     12/20 Test #12: main_chain_triqsqmc_so ...........   Passed   31.01 sec
           Start 13: check_chain_hubbardI
     13/20 Test #13: check_chain_hubbardI .............   Passed    0.73 sec
           Start 14: check_chain_hubbardI_so
     14/20 Test #14: check_chain_hubbardI_so ..........   Passed    0.77 sec
           Start 15: check_chain_triqsqmc
     15/20 Test #15: check_chain_triqsqmc .............   Passed    0.75 sec
           Start 16: check_chain_triqsqmc_so
     16/20 Test #16: check_chain_triqsqmc_so ..........   Passed    0.79 sec
           Start 17: post_chain_hubbardI
     17/20 Test #17: post_chain_hubbardI ..............   Passed    2.70 sec
           Start 18: post_chain_hubbardI_so
     18/20 Test #18: post_chain_hubbardI_so ...........   Passed    2.16 sec
           Start 19: post_chain_triqsqmc
     19/20 Test #19: post_chain_triqsqmc ..............   Passed    2.40 sec
           Start 20: post_chain_triqsqmc_so
     20/20 Test #20: post_chain_triqsqmc_so ...........   Passed    2.14 sec

     100% tests passed, 0 tests failed out of 20

     Total Test time (real) =  69.77 sec


   Some tests (from 9 to 20) may be failed, if some impurity solvers are not installed. The tests from 1 to 8 should be passed regardless of impurity solvers you installed.

#. Finally, install by

   ::

     $ make install

   ``DCore`` is installed in the directory *path_to_triqs* (or *install_path* if specified).

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
