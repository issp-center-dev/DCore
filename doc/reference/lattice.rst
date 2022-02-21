:orphan:

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

  Place the Wannier90 file in the current directory with the name *seedname*\_hr.dat. For details of the Wannier90 format, see :ref:`SciPost 2021 <dcore_paper>`.

For experts, the lattice data may be prepared by your own. In this case, use

* **external**

  In this mode, you should make all necessary data in ``dft_input`` group of *seedname*.h5.
  The data structure follows **DFTTools**. For details, see
  `the reference manual of DFTTools <https://triqs.github.io/dft_tools/1.4/reference/h5structure.html>`_.


  The pre-process ``dcore_pre`` does not touch the data in ``dft_input`` group, and write only additional data such as interactions into ``DCore`` group.
