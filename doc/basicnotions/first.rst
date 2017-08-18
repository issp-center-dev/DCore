What you should know
====================

Probably, you can hardly wait to perform your first DFT+DMFT calculation
with the :program:`DFTTools` package. This documentation and user guide
should make it as easy as possible to get started quickly.
However, it is mutually important to sort out a few prerequisites first.

What is :program:`DFTTools`?
----------------------------

:program:`DFTTools` connects the :ref:`TRIQS <triqslibs:welcome>` library
to realistic materials calculations based
on density functional theory (DFT). It allows an efficient implementation
of DFT plus dynamical mean-field theory (DMFT) calculations and it supplies
tools and methods to construct Wannier functions and to perform the
DMFT self-consistency cycle in this basis set. Post-processing tools,
such as band-structure plotting or the calculation of transport properties
are also implemented. The package comes with a fully charge self-consistent
interface to the Wien2k band structure code, as well as a generic interface.
We assume that you are already know about DFT and the usage of Wien2k.

Have a look at :ref:`DFT+DMFT page <dftplusdmft>` for a brief introduction on
the DFT+DMFT method and on how the theory is reflected in the
:ref:`basic structure <structure>` of the :program:`DFTTools` package.


Understand the philosophy of :program:`DFTTools`
------------------------------------------------

The purpose of :program:`DFTTools` is to provide the necessary tools
required for a DFT+DMFT calculations. Putting those tools together to a working
DFT+DMFT implementation is the task of the user. We do not
supply an universal script which runs with the click of a button, simply because
each material requires a different treatment or different settings.
Building your own script offers a great deal of flexibility and customizability.
Additionally, the :ref:`DFTTools user guide <documentation>` is there to support you
during this process.

It should go without saying, but the verification of outputs and the inspection
of results on their meaningfulness is the responsibility of the user.

The :program:`DFTTools` package is a toolbox and **not** a black box!


Learn how to use :ref:`TRIQS <triqslibs:welcome>` (and the :ref:`CTHYB <triqscthyb:welcome>` solver)
----------------------------------------------------------------------------------------------------

As :program:`DFTTools` is a :ref:`TRIQS <triqslibs:welcome>` based application
it is beneficial to invest a few hours to become familiar with
the :ref:`TRIQS <triqslibs:welcome>` basics first. The
:ref:`TRIQS tutorial <triqslibs:tutorials>` covers
the most important aspects of :ref:`TRIQS <triqslibs:welcome>`. We recommend
downloading our hands-on training in the form of ipython notebooks from
the `tutorials repository on GitHub <https://github.com/TRIQS/tutorials>`_.
Tutorials 1 to 6 are on the :ref:`TRIQS <triqslibs:welcome>` library, whereas tutorials
7 and 8 are more specific to the usage of the :ref:`CTHYB <triqscthyb:welcome>`
hybridization-expansion solver. In general, those tutorials will take at least a full day to finish.

Afterwards you can continue with the :ref:`DFTTools user guide <documentation>`.


Maximum Entropy (MaxEnt)
------------------------

Analytic continuation is needed for many :ref:`post-processing tools <analysis>`, e.g. to
calculate the spectral function, the correlated band structure (:math:`A(k,\omega)`)
and to perform :ref:`transport calculations <Transport>`.
You can use the Pade approximation available in the :ref:`TRIQS <triqslibs:welcome>` library, however,
it turns out to be very unstable for noisy numerical data. Most of the time, the MaxEnt method
is used to obtain data on the real-frequency axis. At the moment neither :ref:`TRIQS <triqslibs:welcome>` nor
:program:`DFTTools` provide such routines.
