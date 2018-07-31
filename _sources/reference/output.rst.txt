.. _outputformat:

Output-file format
==================

``dcore_pre``
~~~~~~~~~~~~~

This program generates *seedname*.h5 file.
It has two groups, ``dft_input`` and ``Dcore``.
See `DFTTools <https://triqs.ipht.cnrs.fr/applications/dft_tools/reference/h5structure.html>`_ for the details of the data structure in ``dft_input`` group.
The data included in ``Dcore`` group is list below:

.. include:: dcore_pre.txt

``dcore``
~~~~~~~~~

This program generates *seedname*.out.h5 file.
All data are stored in ``dmft_out`` group.
The following list summarizes data structure in the ``dmft_out`` group:

.. include:: dcore_out.txt

See
`TRIQS manual <https://triqs.ipht.cnrs.fr/1.x/reference/gfs/py/full.html#hdf5>`_
for the data structure of the Green's function and self-energy.
