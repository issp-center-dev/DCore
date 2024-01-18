Exact diagonalization solver: ``HPhi``
=========================================

``HPhi`` is a numerical solver package that contains the Lanczos method implemented in c++.
``DCore`` provides an interface to ``HPhi`` to compute an approximate solution of DMFT with discritized hybridization function.

Features
--------

- All interactions available in ``DCore`` are supported.

- Low temperature

- More orbitals and bath sites can be treated compared to the full diagonalization solver ``pomerol``.

Install
-------

The following library/program needs to be installed:

- `HPhi <https://github.com/issp-center-dev/HPhi>`_

How to use
----------

Mandatory parameters:

::

    [impurity_solver]
    name = HPhi
    exec_path{str} = /install_directory/bin/pomerol2dcore

Optional parameters that are common to ``pomerol`` solver (see `/impuritysolvers/pomerol/pomerol`_ for explanations):

::

    n_bath{int} = 3  # 0 for default
    fit_gtol{float} = 1e-6  # 1e-5 for default

Optional parameters for ``HPhi`` solver:

::

    np{int} = 1
    Lanczos_max{int} = 2000
    exct{int} = 10
    LanczosTarget{int} = 2
    LanczosEps{int} = 10

Example
-------

:doc:`The square-lattice model in tutorial <../../tutorial/square/square>` is solved by the HPhi solver using the following input parameter set:

:download:`dmft_square_hphi.ini <dmft_square_hphi.ini>`

.. literalinclude:: dmft_square_hphi.ini
   :language: ini

It is recommended to set ``convergence_tol`` parameter in [control] block to stop the DMFT loop automatically.
The figure below shows the renormalization factor as a function of ``n_bath``.
Convergence to the CTHYB result is obtained around ``n_bath=3``.

.. image:: renorm.png
   :width: 700
   :align: center


Trouble shooting
----------------

Warning

::

    Warning: At T = 0.1, exp[-beta(ene_max-ene_mix)]=1.00e+00 is larger than eta=0.0001.

You need more eigenstates by increasing the value of ``exct``. But, please keep in mind that the Lanczos method should be used at lower temperatures.
