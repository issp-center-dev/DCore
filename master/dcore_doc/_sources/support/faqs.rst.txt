.. Frequently-Asked Questions
.. ==========================

FAQ/Troubleshooting
===================

.. contents::
   :local:
   :depth: 2


General
-------

Can **DCore** compute the double occupancy?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It depends on impurity solvers.
For example, ALPSCore/CT-HYB-SEGMENT supports it, but other CT-QMC codes do not at present.

In general, capability of computing **local** physical quantities totally depends on impurity solvers.
See the documentation of each solver for which local quantities can be computed.


Can **DCore** compute the internal energy?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

No, it is not implemented yet.

If you want to compute it, you need to access to G(k,iw) in HDF5 file, and take summation over k and iw in an appropiate way.


..
    Installation
    ------------






``dcore_pre``
-------------

When should I re-execute ``dcore_pre``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you change any parameters in **[model] block**.
``dcore_pre`` will re-construct the model database in *seedname*.h5.

For changes of parameters other than [model] block, on the other hand, you can restart ``dcore``, skipping ``dcore_pre``.


``dcore``
---------

How should I judge convergence of the DMFT loop?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use the ``dcore_check`` and look at the generated figures.
We do not implement automatic convergence check, because results by QMC solvers include statistical errors and simple convergence criteria does not work.

The DMFT loop does not converge
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the chemical potential and/or self-energy look oscillating as a function of the iteration number, decrease ``sigma_mix`` parameter in [control] block.

You might also need to improve the accuracy of QMC sampling by increasing the measurement time.

If the oscillation persists, consider expanding the unitcell to address a symmetry broken solution.


Can I enforce zero magnetic moment?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes. Use ``time_reversal`` flag in [control] block.
See :doc:`../reference/input`.


Can I fix the sign of the magnetic moment in ordered states?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes. Set the initial static self-energy using ``initial_static_self_energy`` parameter in [control] block.
This works as a local potential in the first iteration.
See :doc:`../reference/input` for details.


Can I set an initial self-energy for the DMFT loop?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes. Use ``initial_self_energy`` parameter in [control] block.
See :doc:`../reference/input`.

Let's say you want to start the DMFT loop using results obtained for different parameter set.
To do this, you perform ``dcore_check`` after convergence, and copy **sigma.dat** into a new working directory (with the name sigma_init.dat below).
Then, assign this file to ``initial_self_energy`` parameter as

::

    [control]
    initial_self_energy = sigma_init.dat



``dcore_check``
---------------

Can I change the file format for figures?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes. Use ``--ext`` option to specify the file extension.
If you prefer PDF, for example, the command reads

::

    $ dcore_check input.init --ext pdf

The argument of ``--ext`` can be png, pdf, eps, jpg, ... (the ones supported by ``matplotlib`` module).



Impurity solvers
----------------

What is the difference between ALPSCore/CT-HYB and TRIQS/cthyb?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``ALPS/CT-HYB`` is developed by one of the main developers of DCore, H. Shinaoka.
Both of ``ALPS/CT-HYB`` and ``TRIQS/cthyb`` implement the hybridization-expansion continuous-time quantum Monte Carlo method.
The main difference is the reliability of measurement of the single-particle Green's function.
ALPSCore/CT-HYB uses a more elaborate algorithm (worm sampling).
The non-worm conventional sampling, which is implemented in ``TRIQS/cthyb``,
may give wrong results in some situations (e.g. SOI coupling with orbital-diagonal bath).

Solver-dependent parameter is not recognized
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The variable type should be specified to be recognized as a solver parameter.
For example, integer variable with name *num* is written as

::

    [impurity_solver]
    num{int} = 100

The type can be int, str, and float.

Can I use my impurity solver?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Yes. Follow the instruction in :doc:`../impuritysolvers/how_to_integrate`





..
   ``dcore`` crashes abnormally when using cthyb
   ---------------------------------------------

   Please retry.
