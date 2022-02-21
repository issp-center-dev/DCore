:orphan:

``lattice`` parameter
^^^^^^^^^^^^^^^^^^^^^^

.. You can choose the type of lattice by setting ``lattice``.

Predefined models
-----------------

For model calculations, the following preset models are defined

::

    lattice = chain
    lattice = square
    lattice = cubic
    lattice = bethe

The first three models use tight-binding parameters up to second nearest neighbors. One can specify them by

::

    t = 1
    t' = -0.1

The last one, ``bethe``, use only parameter ``t``, which defines the width (4t) of the semicircular density of states.

.. image:: model.png
   :width: 700
   :align: center


Wannier90 format
----------------

For DFT+DMFT calculations and model calculations of other lattice structures, hopping parameters in the Wannier90 format can be imported by

::

    lattice = wannier90

Place the Wannier90 file in the current directory with the name *seedname*\_hr.dat. For details of the Wannier90 format, see :ref:`SciPost 2021 <dcore_paper>`.


Connection to external programs
-------------------------------

For experts, the lattice data may be prepared by your own. In this case, use

::

    lattice = external

In this mode, you should make all necessary data in ``dft_input`` group in *seedname*.h5.
The data structure follows **DFTTools**. For details, see
`the reference manual of DFTTools <https://triqs.github.io/dft_tools/1.4/reference/h5structure.html>`_.


The pre-process ``dcore_pre`` does not touch the data in ``dft_input`` group, and write only additional data such as interactions into ``DCore`` group.
