Input-file format
=================

The input file consists of six parameter blocks named [model], [system], [impurity\_solver], [control], [tool], and [mpi].

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

For example, we can see that ``dcore_pre`` reads only [model] block. Therefore, ``dcore_pre`` needs to be re-executed only when the parameters in [model] block are changed.

The parameters included in each block are explained below.

.. contents::
    :local:

[model] block
-------------

This block includes parameters for defining a model to be solved.

.. ``dcore_pre``, ``dcore_check`` and ``dcore_post`` read this block.

.. include:: model_desc.txt


See separate pages from the link below for detailed descriptions of some parameters.

:doc:`lattice`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:doc:`interaction`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:doc:`slater_basis`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:doc:`local_potential`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


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
