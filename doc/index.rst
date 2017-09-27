.. index:: DFTTools

.. module:: pytriqs.applications.pydmft

.. _pydmft:

pyDMFT
======

This :ref:`TRIQS-based <triqslibs:welcome>` application is aimed
at ab-initio calculations for 
correlated materials, combining realistic DFT band-structure
calculations with the dynamical mean-field theory. Together with the
necessary tools to perform the DMFT self-consistency loop for
realistic multi-band problems, the package provides a full-fledged
charge self-consistent interface to the `Wien2K package
<http://www.wien2k.at>`_. In addition, if Wien2k is not available, it
provides a generic interface for one-shot DFT+DMFT calculations, where
only the single-particle Hamiltonian in orbital space has to be
provided.

Learn how to use this package in the :ref:`documentation`.


.. toctree::
   :maxdepth: 2
