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

This program can be used for checking the convergence of the SCF-cycle.
This program reads ``[model]`` and ``[tool]`` block.

::

   $ dcore_check input-file

``dcore_check`` shows the history of the chemical potential and the
first component of the self energy at imaginary frequency, :math:`\Sigma_{0 0}(i \omega_n)`
at the last seven iterations.

.. image:: ../tutorial/square/convergence.png
   :width: 500
   :align: center

Post-processing : ``dcore_post``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This program compute the total and *k*\ -resolved spectrul function from the outputted
HDF5 file (*seedname*.out.h5).
This program reads ``[model]``, ``[system]``, ``[impurity-solver]`` and ``[tool]`` block.

::

   $ dcore_post input-file

List of input and output
------------------------

================= ================================================== ====================
Program           Input parameters                                   Output files
================= ================================================== ====================
``dcore_pre``     [model], [system]                                  *seedname*.h5
``dcore``         [model], [system], [impurity-solver], [control]    *seedname*.out.h5
``dcore_check``   [model], [tool]                                    ---
``dcore_post``    [model], [system], [impurity-solver], [tool]       ---
================= ================================================== ====================

Input-file format
-----------------

.. literalinclude:: ../tutorial/bethe-t2g/dmft_bethe.ini

[model] block
~~~~~~~~~~~~~

dcore_pre, dcore_check and dcore_post read this block.

.. include:: model_desc.rst

Prepare model for DCore.
Wannier90 as well as the following preset models:

* chain

* square

* cubic

* bethe
  Semicircular DOS with energy ranges [-2t:2t].

* wannier90
  Read hopping parameters from the Wannier90 output.

.. math::

   {\hat H} = \sum_{i j} \sum_{\alpha \beta}^{N_{\rm band}} \sum_{\sigma=\uparrow, \downarrow}
   t_{i \alpha j \beta} c_{i \alpha \sigma}^\dagger c_{j \beta \sigma}
   +h.c. + {\hat H}_{\rm int},

where :math:`{\hat H}_{\rm int}` is the Kanamori interaction term
given by

.. math::

    {\hat H}_{\rm int} = \sum_i \Big[
    \sum_{\alpha} U n_{i\alpha\uparrow} n_{i\alpha\downarrow}
    + \frac{1}{2} \sum_{\alpha \neq \beta, \sigma \sigma'}
    (U' - J\delta_{\sigma\sigma'}) n_{i\alpha\sigma} n_{i\beta\sigma'}
    - \sum_{\alpha \neq \beta} J(c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow} c_{i\beta\downarrow}^{\dagger} c_{i\beta\uparrow}
    + c_{i\alpha\uparrow}^{\dagger} c_{i\alpha\downarrow}^{\dagger} c_{i\beta\uparrow} c_{i\beta\downarrow}
    \Big]

The parameter :math:`U'` is fixed at :math:`U'=U-2J`.

.. image:: model.png
   :width: 700
   :align: center

[system] block
~~~~~~~~~~~~~~

dcore_pre and dcore read this block.

.. include:: system_desc.rst

[impurity_solver] block
~~~~~~~~~~~~~~~~~~~~~~~

dcore and dcore_post read this block.

.. include:: impurity_solver_desc.rst

**... and other parameters (Solver dependent)**.
We have to specify additional parameters with types (e.g. ``n_cycles{int} = 500000``).
For more details, please see the reference page of
`TRIQS/cthyb <https://triqs.ipht.cnrs.fr/applications/cthyb/reference/solve_parameters.html#solve-parameters>`_,
`ALPS/cthyb <https://github.com/shinaoka/triqs_interface#program-parameters>`_, etc..

[control] block
~~~~~~~~~~~~~~~

dcore reads this block.

.. include:: control_desc.rst

[tool] block
~~~~~~~~~~~~

dcore_check and dcore_post reads this block.

.. include:: tool_desc.rst
