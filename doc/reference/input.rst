.. _inputformat:

Input-file format and usage
===========================

Usage
-----

The following programs can read the same input file.

Pre-processing : ``dcore_pre``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program generates model HDF5 file (*seedname*.h5).
It must be executed before the main program, ``dcore``, runs.
This program reads ``[model]`` and ``[system]`` block.

::

   $ dcore_pre input-file

Main program : ``dcore``
~~~~~~~~~~~~~~~~~~~~~~~~

This program performs DMFT cycle and output the self energy etc. into a HDF
file (*seedname*.out.h5).
This program reads ``[model]``, ``[system]``, ``[impurity-solver]`` and ``[control]`` block.

::

   $ dcore input-file

Convergence-check : ``dcore_check``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program can be used for checking the convergence of the DMFT loop.
This program reads ``[model]`` and ``[tool]`` block.

::

   $ dcore_check input-file

``dcore_check`` shows the history of the chemical potential and the
average of the self energy at imaginary frequency,

.. math::

   \Sigma_{\rm Ave} (i \omega_n) = 
   \left[\sum_i^{\rm shell} \sum_{\alpha \beta}^{N_{\rm orb}^i} \Sigma_{\alpha \beta}(i\omega)\right]
   /\left[\sum_i^{\rm shell} N_{\rm orb}^{i}\right],

at the last seven iterations.

.. image:: ../tutorial/square/convergence.png
   :width: 500
   :align: center

The maximum frequency of this plot is specified with the parameter ``omega_check``
in the ``[tool]`` block.

Also, this program generates a text file, *seedname*\ ``_sigma.dat``, which contains
the local self energy at the final step as follows:

::

   # Local self energy at imaginary frequency
   # [Column] Data
   # [1] Frequency
   # [2] Re(Sigma_{shell=0, spin=up, 0, 0})
   # [3] Im(Sigma_{shell=0, spin=up, 0, 0})
   # [4] Re(Sigma_{shell=0, spin=down, 0, 0})
   # [5] Im(Sigma_{shell=0, spin=down, 0, 0})
   -157.001093 0.994751 0.006358 0.994751 0.006358
   -156.844013 0.994751 0.006365 0.994751 0.006365
   -156.686934 0.994751 0.006371 0.994751 0.006371
   :
           
Post-processing : ``dcore_post``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program compute the total and *k*\ -resolved spectral function from the outputted
HDF5 file (*seedname*.out.h5).
This program reads ``[model]``, ``[system]``, ``[impurity-solver]`` and ``[tool]`` block.

::

   $ dcore_post input-file

List of input and output
------------------------

================= ================================================== ====================
Program           Input parameters                                   Output HDF files
================= ================================================== ====================
``dcore_pre``     [model], [system]                                  *seedname*.h5
``dcore``         [model], [system], [impurity-solver], [control]    *seedname*.out.h5
``dcore_check``   [model], [tool]                                    ---
``dcore_post``    [model], [system], [impurity-solver], [tool]       ---
================= ================================================== ====================

Input-file format
-----------------

.. literalinclude:: ../tutorial/nis/nis.ini
   :language: ini

[model] block
~~~~~~~~~~~~~

dcore_pre, dcore_check and dcore_post read this block.

.. include:: model_desc.txt

Prepare model for DCore.
Wannier90 as well as the following preset models:

* chain

* square

* cubic

* bethe
  Semicircular DOS with energy ranges [-2t:2t].

* wannier90
  Read hopping parameters from the Wannier90 output.

.. image:: model.png
   :width: 700
   :align: center

.. math::

   {\hat H} = \sum_{i j} \sum_{\alpha \beta}^{N_{\rm band}} \sum_{\sigma=\uparrow, \downarrow}
   t_{i \alpha j \beta} c_{i \alpha \sigma}^\dagger c_{j \beta \sigma}
   +h.c. + {\hat H}_{\rm int},

where

.. math::

   {\hat H}_{\rm int} = \frac{1}{2}
   \sum_{i, \alpha \beta \gamma \delta,\sigma \sigma'}
   U^{i}_{\alpha \beta \gamma \delta}
   c_{i \alpha \sigma}^\dagger c_{i \beta \sigma'}^\dagger c_{i \delta \sigma'} c_{i \gamma \sigma}.

The interaction matrix :math:`U^{i}_{\alpha \beta \gamma \delta}`
is specified by the parameter ``interaction``.

* If ``interaction = kanamori``

  In this case, the Kanamori-type interaction is used, i.e.

  .. math::
     :nowrap:
     
     \begin{align}
     U_{\alpha \alpha \alpha \alpha} &= U,
     \\
     U_{\alpha \beta \alpha \beta} &= U' \qquad (\alpha \neq \beta),
     \\
     U_{\alpha \beta \beta \alpha} &= J \qquad (\alpha \neq \beta),
     \\
     U_{\alpha \alpha \beta \beta} &= J \qquad (\alpha \neq \beta),
     \end{align}
  
  where :math:`U, U', J` at each correlated shell are specified by the parameter ``kanamori`` as

  ::

     interaction = kanamori
     kanamori = [(U_1, U'_1, J_1), (U_2, U'_2, J_2), ... ]

  For example, if there are two correlated shells that have
  :math:`(U, U', J) = (4, 2, 1)` and :math:`(U, U', J) = (6, 3, 1.5)`, respectovely,
  the input parameter becomes

  ::

     interaction = kanamori
     kanamori = [(4.0, 2.0, 1.0), (6.0, 3.0, 1.5)]

* If ``interaction = slater_f``

  In this case, the interaction matrix is constructed by the effective Slater integrals
  :math:`F_0, F_2, F_4, F_6`.
  These Slater integrals amd angular mementum at each correlated shell
  are specified by the parameter ``slater_f`` as follows

  ::

     interaction = slater_f
     slater_f = [(angular_momentum, F_0, F_2, F_4, F6), ... ]

  For example, if there are two correlated shells,
  one has d-orbital with :math:`(F_0, F_2, F_4) = (2, 1, 0.5)` amd
  the other has p-orbital with :math:`(F_0, F_2) = (3, 1.5)`,
  the input parameter becomes

  ::

     interaction = slater_f
     slater_f = [(2, 2.0, 1.0, 0.5, 0.0), (1, 3.0, 1.5 0.0, 0.0)]

  .. note::
   
     Even when we compute the s, p, d orbital, we have to specify all of :math:`F_0, F_2, F_4, F_6`.
     
* If ``interaction = slater_uj``
      
  Also in this case, the Slater-type interaction is used.
  The effective Slater integrals are computed with the following formulae:

  * :math:`l = 1`

    .. math::
   
       F_0 = U, \quad
       F_2 = 5 J

  * :math:`l=2`

    .. math::
   
       F_0 = U, \quad
       F_2 = \frac{14 J}{1.0 + 0.63},\quad
       F_4 = 0.63 F_2
  
  * :math:`l=3`

    .. math::
   
       F_0 = U, \quad
       F_2 = \frac{6435 J}{286 + 195 \times 451 / 675 + 250 \times 1001 / 2025},\quad
       F_4 = \frac{451 F_2}{675},\quad
       F_6 = \frac{1001 F_2}{2025}

  The :math:`U`, :math:`J` and the angular momentum at each correlated shell
  are specified by the parameter ``slater_uj`` as

  ::

     interaction = slater_uj
     slater_uj = [(angular_momentum1, U1, J1), (angular_momentum2, U2, J2), ... ]
  
* If ``interaction = respack``

  Use the output by `RESPACK <https://sites.google.com/view/kazuma7k6r>`_.
  **Under construction.**

If we want to treat only the density-density part

.. math::

   {\hat H}_{\rm int} = \frac{1}{2}
   \sum_{i, \alpha, \sigma \sigma'}
   U^{i}_{\alpha \alpha \alpha \alpha}
   c_{i \alpha \sigma}^\dagger c_{i \beta \sigma'}^\dagger c_{i \beta \sigma'} c_{i \alpha \sigma}
   + \frac{1}{2}
   \sum_{i, \alpha \neq \beta, \sigma \sigma'}
   U^{i}_{\alpha \beta \alpha \beta}
   c_{i \alpha \sigma}^\dagger c_{i \beta \sigma'}^\dagger c_{i \beta \sigma'} c_{i \alpha \sigma}
   + \frac{1}{2}
   \sum_{i, \alpha \neq \beta, \sigma}
   U^{i}_{\alpha \beta \beta \alpha}
   c_{i \alpha \sigma}^\dagger c_{i \beta \sigma}^\dagger c_{i \alpha \sigma} c_{i \beta \sigma},

we specify the parameter ``density-density`` as

::

   density-density = True

.. note::

   It can not be used in conjunction to the Hubbard-I solver or
   the double-counting correction.

[system] block
~~~~~~~~~~~~~~

dcore_pre and dcore read this block.

.. include:: system_desc.txt

If the parameter ``with_dc`` is specified to ``True``,
the following part of the self-energy is subtracted to avoid the double-counting error of
the self energy.

.. math::

   \Sigma_{i, \alpha \sigma \beta \sigma'}^{\rm dc-imp}
   = \delta_{\sigma \sigma'} \sum_{\gamma \delta \sigma_1}
   U_{\alpha \gamma \beta \delta}
   \langle c_{\gamma \sigma_1}^\dagger c_{\delta \sigma_1}\rangle_0
   - \sum_{\gamma \delta}
   U_{\alpha \gamma \delta \beta}
   \langle c_{\gamma \sigma'}^\dagger c_{\delta \sigma}\rangle_0,

where :math:`\langle \cdots \rangle_0` indicates the expectation value at the initial (Kohn-Sham) state.
   
[impurity_solver] block
~~~~~~~~~~~~~~~~~~~~~~~

dcore and dcore_post read this block.

.. include:: impurity_solver_desc.txt

**... and other parameters (Solver dependent)**.
We have to specify additional parameters with types (e.g. ``n_cycles{int} = 500000``).
For more details, please see the reference page of
`TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/reference/solve_parameters.html#solve-parameters>`_,
`ALPS/cthyb <https://github.com/shinaoka/triqs_interface#program-parameters>`_, etc..

[control] block
~~~~~~~~~~~~~~~

dcore reads this block.

.. include:: control_desc.txt

[tool] block
~~~~~~~~~~~~

dcore_check and dcore_post reads this block.

.. include:: tool_desc.txt
