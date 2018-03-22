.. _inputformat:

Input-file format
=================

The input file is constructed of five blocks, ``[model]``, ``[system]``, ``[impurity]``, ``[control]`` and ``[tool]``.
The example of the input file is shown as follows:

.. literalinclude:: ../tutorial/nis/nis.ini
   :language: ini

The details of each block will be described below.
	      
[model] block
~~~~~~~~~~~~~

``dcore_pre``, ``dcore_check`` and ``dcore_post`` read this block.

.. include:: model_desc.txt

A model for DCore is defined in this block.
You can choose the type of lattice by setting ``lattice``.
The following preset models are defined:

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

Model Hamiltonian is defined as
   
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
  :math:`(U, U', J) = (4, 2, 1)` and :math:`(U, U', J) = (6, 3, 1.5)`, respectively,
  you need to set the input parameters as

  ::

     interaction = kanamori
     kanamori = [(4.0, 2.0, 1.0), (6.0, 3.0, 1.5)]

* If ``interaction = slater_f``

  In this case, the interaction matrix is constructed by the effective Slater integrals 
  :math:`F_0, F_2, F_4, F_6`.
  These Slater integrals and the angular momentum at each correlated shell
  are specified by the parameter ``slater_f`` as follows

  ::

     interaction = slater_f
     slater_f = [(angular_momentum, F_0, F_2, F_4, F6), ... ]

  For example, if there are two correlated shells,
  one has d-orbital with :math:`(F_0, F_2, F_4) = (2, 1, 0.5)` and
  the other has p-orbital with :math:`(F_0, F_2) = (3, 1.5)`,
  you need to set the input parameter as

  ::

     interaction = slater_f
     slater_f = [(2, 2.0, 1.0, 0.5, 0.0), (1, 3.0, 1.5 0.0, 0.0)]

  .. note::
     
     You must specify all of :math:`F_0, F_2, F_4, F_6`.
     
* If ``interaction = slater_uj``
      
  In this case, the Slater-type interaction is used.
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

If you want to treat only the density-density part

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

you specify the parameter ``density-density`` as

::

   density-density = True

.. note::

   It can not be used in conjunction to the Hubbard-I solver or
   the double-counting correction.

[system] block
~~~~~~~~~~~~~~

``dcore_pre`` and ``dcore`` read this block.

.. include:: system_desc.txt

If the parameter ``with_dc`` is specified to ``True``,
the following part of the self-energy is subtracted to avoid the double-counting error of
the self-energy.

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

dcore_check and dcore_post read this block.

.. include:: tool_desc.txt
