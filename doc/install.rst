
.. highlight:: bash

Installation
============


Prerequisites
-------------

#. The :ref:`TRIQS <triqslibs:welcome>` and
   :ref:`DFT-Tools <triqsdfttools:dft>` toolbox.
   In the following, we will suppose that it is installed in the ``path_to_triqs`` directory.

#. Likely, you will also need at least one impurity solver,
   e.g. the `Hubbard-I solver <https://triqs.ipht.cnrs.fr/1.x/applications/hubbardI/>`_.

Installation steps 
------------------

#. Download the sources from github:: 
 
     $ git clone https://github.com/issp-center-dev/DCore.git src
 
#. Create an empty build directory where you will compile the code:: 
 
     $ mkdir build && cd build 
 
#. In the build directory call cmake specifying where the TRIQS library is installed:: 
 
     $ cmake -DTRIQS_PATH=path_to_triqs ../src 
 
#. Compile the code, run the tests and install the application:: 
 
     $ make 
     $ make test 
     $ make install 
 
Version compatibility 
--------------------- 
 
Be careful that the version of the TRIQS library and of the dft tools must be 
compatible (more information on the :ref:`TRIQS website <triqslibs:welcome>`). 
If you want to use a version of the dft tools that is not the latest one, go
into the directory with the sources and look at all available versions:: 
 
     $ cd src && git tag 
 
Checkout the version of the code that you want, for instance:: 
 
     $ git co 1.4
 
Then follow the steps 2 to 5 described above to compile the code. 
