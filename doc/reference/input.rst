Input-file format
=================

The input file consists of five parameter blocks named ``[model]``, ``[system]``, ``[impurity_solver]``, ``[control]``, and ``[tool]``.

..
    The example of the input file is shown as follows:

    .. literalinclude:: ../tutorial/nis/nis.ini
       :language: ini

The following table shows which blocks are used by each program.

.. csv-table::
    :header: "", ``dcore_pre``, ``dcore``, ``dcore_check``, ``dcore_post``
    :widths: 5, 5, 5, 5, 5

    [model]          ,  Yes, Yes, Yes, Yes
    [system]         ,     , Yes, Yes, Yes
    [impurity_solver],     , Yes,    , Yes
    [control]        ,     , Yes
    [tool]           ,     ,    , Yes, Yes
    [mpi]            ,     , Yes,    , Yes

For example, we can see that ``dcore_pre`` needs to be re-executed only when [model] block is changed.

The parameters included in each block are explained below.

[model] block
-------------

This block includes parameters for defining a model to be solved.

.. ``dcore_pre``, ``dcore_check`` and ``dcore_post`` read this block.

.. include:: model_desc.txt

lattice
^^^^^^^

.. You can choose the type of lattice by setting ``lattice``.

For model calculations, the following preset models are defined:

* **chain**

* **square**

* **cubic**

* **bethe** (semicircular DOS with energy ranges [-2t:2t])

.. image:: model.png
   :width: 700
   :align: center

For DFT+DMFT calculations, hopping parameters in the Wannier90 format can be imported by

* **wannier90**

  Place the Wannier90 file in the current directory with the name *seedname*\_hr.dat.

For experts, the lattice data may be prepared by your own. In this case, use

* **external**

  In this mode, you should make all necessary data in ``dft_input`` group of *seedname*.h5.
  The data structure follows **DFTTools**. For details, see
  `the reference manual of DFTTools <https://triqs.github.io/dft_tools/1.4/reference/h5structure.html>`_.


  The pre-process ``dcore_pre`` does not touch the data in ``dft_input`` group, and write only additional data such as interactions into ``DCore`` group.

interaction
^^^^^^^^^^^

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
     slater_f = [(angular_momentum, F_0, F_2, F_4, F_6), ... ]

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

you specify the parameter ``density_density`` as

::

   density_density = True

.. note::

   It can not be used in conjunction to the Hubbard-I solver or
   the double-counting correction.

local potential
^^^^^^^^^^^^^^^

An arbitrary local potential can be implemented using parameters ``local_potential_*``.
The format looks like

::

    [model]
    local_potential_matrix = {0: 'pot0.txt', 1: 'pot1.txt'}
    local_potential_factor = 0.01

Here, ``local_potential_matrix`` describes, in the python dictionary format, a set of the inequivalent shell index *ish* and the filename which defines the local potential matrix.
The parameter ``local_potential_factor`` defines a prefactor to the potential matrix.

For example, the Zeeman term along z-axis for S=1/2 is represented by

::

    $ cat pot0.txt
    # spin orb1 orb2  Re Im
    0 0 0   0.5 0.
    1 0 0  -0.5 0.

and the magnetic field is specified by ``local_potential_factor``.

[system] block
--------------

This block includes thermodynamic parameters and some technical parameters such as the number of Matsubara frequencies.

.. ``dcore`` read this block.

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
-----------------------

This block specifies an impurity solver to be used and necessary parameters for running the solver program.

.. ``dcore`` and ``dcore_post`` read this block.

.. include:: impurity_solver_desc.txt

Additionally, we have to specify solver-dependent parameters in the way like ``n_cycles{int} = 500000``.
For details, see :doc:`the reference manual for each solver <../impuritysolvers>`.

..
    `TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/reference/solve_parameters.html#solve-parameters>`_,
    `ALPS/cthyb <https://github.com/shinaoka/triqs_interface#program-parameters>`_, etc..

[control] block
---------------

This block includes parameters that control the self-consistency loop of DMFT.

.. ``dcore`` reads this block.

.. include:: control_desc.txt

[tool] block
------------

This block includes parameters that are solely used by ``dcore_post``.

.. ``dcore_check`` and ``dcore_post`` read this block.

.. include:: tool_desc.txt

[mpi] block
------------

This block includes parameters which are read by ``dcore`` and ``dcore_post``.

.. include:: mpi_desc.txt

When an option ``-DMPIEXEC=<MPIRUN>`` is passed to the ``cmake`` command,
The default value of ``command`` will be replaced with ``<MPIRUN>``.
