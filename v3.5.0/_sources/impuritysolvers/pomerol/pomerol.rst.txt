Exact diagonalization solver: ``pomerol``
=========================================

``pomerol`` is an exact diagonalization (ED) library implemented in c++.
``DCore`` provides an interface to ``pomerol`` library to compute an approximate solution of DMFT with discritized hybridization function.

Features
--------

- Arbitrary temperature

- All interactions available in ``DCore`` are supported.

- [experimental] two-particle Green's function

- [todo] tail evaluation of Gf

Install
-------

The following library/program needs to be installed:

- `pomerol library <https://github.com/aeantipov/pomerol>`_

- `pomerol2dcore <https://github.com/j-otsuki/pomerol2dcore>`_

How to use
----------

Mandatory parameters:

::

    [impurity_solver]
    name = pomerol
    exec_path{str} = /install_directory/bin/pomerol2dcore

Optional parameters:

::

    n_bath{int} = 3  # 0 for default
    fit_gtol{float} = 1e-6  # 1e-5 for default

The default value of ``n_bath`` is 0, namely, no bath site is taken into account (Hubbard-I approximation).
For ``n_bath>0``, hybridization function Delta(iw) is fitted by

.. math::

    \Delta^{n_\mathrm{bath}}_{o_1 o_2}(i\omega) = \sum_{l=1}^{n_\mathrm{bath}} \frac{V_{o_1 l} V_{l o_2}}{i\omega - \epsilon_l}

Then, the finite-size system consisting of the impurity site and ``n_bath`` bath sites are solve by ED method.
The size of the Hilbert space increases exponentially according to :math:`2^{n_\textrm{spn-orb}}` where :math:`n_\textrm{spn-orb}=2*n_\mathrm{orb} + 2*n_\mathrm{bath}`.
Because of storage limitation, :math:`n_\textrm{spn-orb} \simeq 16` is the limits in this solver.

Example
-------

:doc:`The square-lattice model in tutorial <../../tutorial/square/square>` is solved by the pomerol solver using the following input parameter set:

:download:`dmft_square_pomerol.ini <dmft_square_pomerol.ini>`

.. literalinclude:: dmft_square_pomerol.ini
   :language: ini

It is recommended to set ``convergence_tol`` parameter in [control] block to stop the DMFT loop automatically.
The figure below shows the renormalization factor as a function of ``n_bath``.
Convergence to the CTHYB result is obtained around ``n_bath=3``.

.. image:: renorm.png
   :width: 700
   :align: center
