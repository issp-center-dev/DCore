
.. highlight:: bash

.. _installation_triqs14x:

Installation
============

#. You first need to install `TRIQS <https://triqs.github.io/triqs/>`_ library.
   ``TRIQS`` will be installed in directory *path_to_triqs*.

   .. note::

     ``DCore`` supports ``TRIQS`` version **1.4.2** (the support for the current version 2.1.x is still experimental).
     Since the installation manual in the TRIQS official website is not compatible with **1.4.2**, we briefly describe below how to install it.
     Please study `the official document for TRIQS 1.4.x <https://triqs.github.io/triqs/1.4.x/>`_.

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

#. You also need `DFTTools <https://triqs.github.io/dft_tools>`_, which runs on the ``TRIQS`` library.

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
