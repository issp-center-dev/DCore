.. _structure:

Structure of :program:`pyDMFT`
==============================

.. image:: images/structure.png
   :width: 700
   :align: center

The central part of :program:`DFTTools`, which is performing the
steps for the DMFT self-consistency cycle, is written following the
same philosophy as the :ref:`TRIQS <triqslibs:welcome>` toolbox. At
the user level, easy-to-use python modules are provided that allow to
write simple and short scripts performing the actual calculation.
The usage of those modules is presented in the user guide of this
:ref:`documentation`. Before considering the user guide, we suggest
to read the following introduction on the general structure of
the :program:`DFTTools` package.

The interface layer
-------------------

Since the input for this DMFT part has to be provided by DMFT
calculations, there needs to be another layer that connects the
python-based modules with the DFT output. Naturally, this layer
depends on the DFT package at hand. At the moment, there is an
interface to the Wien2k band structure package, and a very light 
interface that can be used in a more general setup. Note that this
light interface layer **does not** allow full charge self-consistent
calculations. 

Standard interface
""""""""""""""""""

In addition to the specialized Wien2k interface, :program:`DFTTools`
provides also a very light-weight general interface. It basically
consists of a very simple :class:`HkConverter`. As input it requires a
Hamiltonian matrix :math:`H_{mn}(\mathbf{k})` written already in
localized-orbital indices :math:`m,n`, on a :math:`\mathbf{k}`-point
grid covering the Brillouin zone, and just a few other informations
like total number of electrons, how many correlated atoms in the unit
cell, and so on. It converts this Hamiltonian into a hdf5 format and
sets some variables to standard values, such that it can be used with
the python modules performing the DMFT calculation. How the
Hamiltonian matrix :math:`H_{mn}(\mathbf{k})` is actually calculated,
is **not** part of this interface.

Wannier90 interface
"""""""""""""""""""

In addition to the specialized Wien2k interface, :program:`DFTTools`
provides also a very light-weight general interface. It basically
consists of a very simple :class:`HkConverter`. As input it requires a
Hamiltonian matrix :math:`H_{mn}(\mathbf{k})` written already in
localized-orbital indices :math:`m,n`, on a :math:`\mathbf{k}`-point
grid covering the Brillouin zone, and just a few other informations
like total number of electrons, how many correlated atoms in the unit
cell, and so on. It converts this Hamiltonian into a hdf5 format and
sets some variables to standard values, such that it can be used with
the python modules performing the DMFT calculation. How the
Hamiltonian matrix :math:`H_{mn}(\mathbf{k})` is actually calculated,
is **not** part of this interface.

Post-processing
---------------

The main result of DMFT calculation is the interacting Greens function
and the Self energy. However, one is normally interested in
quantities like band structure, density of states, or transport
properties. In order to calculate these things, :program:`DFTTools`
provides the post-processing modules :class:`SumkDFTTools`. It
contains routines to calculate

* (projected) density of states
* partial charges
* correlated band structures (*spaghettis*)
* transport properties such as optical conductivity, resistivity,
  or thermopower.

.. warning::
   At the moment neither :ref:`TRIQS<triqslibs:welcome>` nor :program:`DFTTools`
   provides Maximum Entropy routines! You can use the Pade
   approximation implemented in the :ref:`TRIQS <triqslibs:welcome>` library, or you use your own
   home-made Maximum Entropy code to do the analytic continuation from
   Matsubara to the real-frequency axis.
