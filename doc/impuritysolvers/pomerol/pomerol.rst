Hubbard-I approximation: ``pomerol``
====================================

``pomerol`` is an exact diagonalization library implemented in c++.
``DCore`` provides an interface to ``pomerol`` library to use it as a solver for the Hubbard-I approximation.

Features
--------

- All interactions available in ``DCore`` are supported.

- [todo] tail evaluation of Gf

- [todo] dynamical susceptibility

- [experimental] two-particle Green's function

Install
-------

The following library/program needs to be installed:

- `pomerol library <https://github.com/aeantipov/pomerol>`_

- `pomerol2dcore <https://github.com/j-otsuki/pomerol2dcore>`_

How to use
----------

::

    [impurity_solver]
    name = pomerol
    exec_path{str} = /install_directory/bin/pomerol2dcore
