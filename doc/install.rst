
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

      $ tar zxvf v2.x.x.tar.gz

- **[Develop version]**

  Newly implemented features are available in the source code in GitHub repository. You can clone the repository by

  .. code-block:: bash

      $ git clone https://github.com/issp-center-dev/DCore.git dcore.src

  The master branch basically corresponds to the latest released package, and the develop branch includes new features.
  Note that develop branches may not be well tested and contains some bugs.

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

   - `ALPS/CT-HYB <https://github.com/ALPSCore/CT-HYB>`_

   - `TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/index.html>`_
     [Note: checkout tag **1.4.2** on master branch, or download the released package from `here <https://github.com/TRIQS/cthyb/releases/tag/1.4.2>`_]
   and others. Please see a complete list of the supported impurity solvers in :doc:`reference/impuritysolvers`.

   .. note::

      If you want to treat the spin-orbit coupling in ``TRIQS/cthyb`` solver,
      it must be built with the following CMake options:

      ::

         -DHYBRIDISATION_IS_COMPLEX=ON -DLOCAL_HAMILTONIAN_IS_COMPLEX=ON

   .. note::

      One must build TRIQS, TRIQS/DFTTools, and TRIQS solvers using the same C++ compiler with the same C++ standard (C++14).
      One does not necessarily have to build ALPS/CT-HYB with the same C++ compiler as that used for TRIQS.

   .. note::

      ``ALPS/CT-HYB`` is developed by one of the main developers of DCore, H. Shinaoka.
      Both of ``ALPS/CT-HYB`` and ``TRIQS/cthyb`` implement the hybridization-expansion continuous-time quantum Monte Carlo method.
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

     $ cmake -DTRIQS_PATH=path_to_triqs path_to_dcore_src -DCMAKE_INSTALL_PREFIX=path_to_dcore_install_directory

   Here, *path_to_triqs* refers to your ``TRIQS`` install directory, and *path_to_dcore_src* refers to the ``DCore`` source directory.
   Please set CMAKE_INSTALL_PREFIX to the directory where DCore will be installed.
   If the cmake command succeeded, you will see the following message.

   ::

     -- Build files have been written to: /.../dcore.build

#. Build DCore by

   ::

     $ make

#. We recommend that you run the tests to check if the compiling is properly finished. Type

   ::

     $ make test

   The test results look like

   ::

     Running tests...
     /usr/local/Cellar/cmake/3.13.4/bin/ctest --force-new-ctest-process
     Test project /Users/hiroshi/build/dcore
           Start  1: typed_parser
      1/12 Test  #1: typed_parser .....................   Passed    1.38 sec
           Start  2: tools
      2/12 Test  #2: tools ............................   Passed    0.75 sec
           Start  3: openmx
      3/12 Test  #3: openmx ...........................   Passed    0.70 sec
           Start  4: respack
      4/12 Test  #4: respack ..........................   Passed    0.81 sec
           Start  5: pre_preset
      5/12 Test  #5: pre_preset .......................   Passed    1.69 sec
           Start  6: pre_wannier
      6/12 Test  #6: pre_wannier ......................   Passed    3.02 sec
           Start  7: pre_wannier_so
      7/12 Test  #7: pre_wannier_so ...................   Passed    0.95 sec
           Start  8: pre_respack
      8/12 Test  #8: pre_respack ......................   Passed    0.94 sec
           Start  9: pre_respack_so
      9/12 Test  #9: pre_respack_so ...................   Passed    1.03 sec
           Start 10: alps_cthyb
     10/12 Test #10: alps_cthyb .......................   Passed    3.02 sec
           Start 11: chain_hubbardI_so
     11/12 Test #11: chain_hubbardI_so ................   Passed   28.28 sec
           Start 12: chain_hubbardI
     12/12 Test #12: chain_hubbardI ...................   Passed   23.49 sec
     
     100% tests passed, 0 tests failed out of 12
   
     Total Test time (real) =  66.11 sec
  
   In the above example, all tests have passed in 66 sec.

#. Finally, install by

   ::

     $ make install

   ``DCore`` is installed in the directory *path_to_dcore_install_directory*.

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
