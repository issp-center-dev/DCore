.. _outputformat:

Output-file format
==================

``dcore_pre``
~~~~~~~~~~~~~

This program generates *seedname*.h5 file.
It has two groups ``dft_input`` and ``Dcore`` .
See `DFTTools <https://triqs.ipht.cnrs.fr/applications/dft_tools/reference/h5structure.html>`_ for the details of the data structure about ``dft_input`` group.
``Dcore`` group has

.. include:: dcore_pre.txt

``dcore`` 
~~~~~~~~~

This program generates *seedname*.out.h5 file.
It has ``dmft_out`` group which is constructed by the following dataset and groups:

.. include:: dcore_out.txt
