
.. highlight:: bash

.. _installation:

Installation
============

Prerequisites
-------------

#. Python3

#. You will also need at least one impurity solver.

   For example, the following programs are supported:

   * :doc:`Hubbard-I solver<impuritysolvers/triqs_hubbard_one/hubbard_one>`
   * :doc:`ALPS/CT-HYB<impuritysolvers/alpscore_cthyb/cthyb>`
   * :doc:`ALPS/CT-HYB-SEGMENT<impuritysolvers/alpscore_ctseg/ctseg>`
   * :doc:`TRIQS/cthyb<impuritysolvers/triqs_cthyb/cthyb>`
   * :doc:`Pomerol solver<impuritysolvers/pomerol/pomerol>`

   We recommend to use the Hubbard-I solver for tests because this is fast.
   See :doc:`impuritysolvers` for a complete list of supported impurity solvers and their user manuals.

#. [OPTIONAL] `TRIQS 3.x <https://triqs.github.io/triqs/>`_ and `TRIQS/DFTTools 3.x <https://triqs.github.io/dft_tools/>`_.

   TRIQS can improve the performance of DCore.
   The current version of DCore supports TRIQS 3.x.
   Please make sure that the triqs and triqs_dft_tools modules are loadable in your Python environment.
   You may use `MateriAppsInstaller <https://github.com/wistaria/MateriAppsInstaller>`_, a collection of install scripts, to install prerequisites (TRIQS).
   If you want to use TRIQS, please set an environment variable ``DCORE_TRIQS_COMPAT`` to 0.

   .. code-block:: bash
         
      $ export DCORE_TRIQS_COMPAT=0


Installation
------------------

(NOTE: Using a virtual environment such as `venv <https://docs.python.org/3/library/venv.html>`_ is recommended to avoid conflicts with other Python packages.)

You can install the latest version of DCore using ``pip`` command as follows.

   ::

      $ pip3 install dcore -U

Here, ``-U`` stands for ''Upgrade''. The installed packages are upgraded to the latest version, if you already have packages installed.


Installation (only for developers)
--------------------------------------

You can download the source files in two ways.

- **[Release version]**

  Go to the `release page <https://github.com/issp-center-dev/DCore/releases>`_ and download the latest tar or zip file.
  You can unpack, for example, the tar file by

  .. code-block:: bash

      $ tar zxvf v3.x.x.tar.gz

- **[Develop version]**

  Newly implemented features are available in the source code in GitHub repository. You can clone the repository by

  .. code-block:: bash

      $ git clone https://github.com/issp-center-dev/DCore.git dcore.src

  The master branch basically corresponds to the latest released package, and the develop branch includes new features.
  Note that develop branches may not be well tested and contains some bugs.


  Please execute the following command in the source directory to install DCore.

  .. code-block:: bash

      $ pip3 install .

  If both of them did not work, you could build a binary package and install it as follows

  .. code-block:: bash

      $ rm -rf dist
      $ python3 setup.py bdist_wheel
      $ pip3 install dist/dcore-*.whl

  One can run unit tests using the installed DCore by executing the following commands.

  Non-MPI tests can be run as follows.

  .. code-block:: bash

     $ pip3 install pytest
     $ pytest tests/non-mpi/*/*.py

  MPI tests can be run as follows.

  .. code-block:: bash

     $ pytest tests/mpi/*/*.py

  MPI tests invoke MPI parallelized impourity solvers.
  If your system MPI command is not "mpirun", please provide the name of the correct one to DCore at runtime in an input file.
  The default value is "mpirun -np #" (# is replaced by the number of processors).
  For instance, if the MPI command of your system is "mpijob", please set the environment variable "DCORE_MPIRUN_COMMAND" as folows.

  .. code-block:: bash

     $ export DCORE_MPIRUN_COMMAND="mpijob -np #"

  Note that it is not allowed to run MPI programs interactively on some system.
  In this case, please run MPI tests as a parallel job with one process.

  You can build documentations as follows.

  .. code-block:: bash

     $ pip3 install sphinx wild_sphinx_theme matplotlib
     $ python3 -m dcore.option_tables doc/reference
     $ sphinx-build -b html doc html
