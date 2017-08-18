.. _dftplusdmft:

Introduction to DFT+DMFT
========================

When describing the physical and also chemical properties of
crystalline materials, there is a standard model that is used with
great success for a large variety of systems: band theory. In simple
terms it states that electrons in a crystal form bands of allowed
states in momentum space. These states are then filled by the
electrons according to Pauli's principle up the Fermi level. With this
simple picture one can explain the electronic band structure of simple
materials such as elementary copper or aluminum.

Following this principle one can easily classify all existing
materials into metals and insulators, with semiconductors being
special insulators with a small gap in the excitation
spectrum. Following this band theory, a system is a metal if there is
an odd number of electrons in the valence bands, since this leads to a
partially filled band, cutting the Fermi energy and, thus, producing a
Fermi surface, i.e metallic behavior. On the other hand, an even
number of electrons leads to completely filled bands with a finite
excitation gap to the conduction bands, i.e. insulating behavior.

This classification works pretty well for a large class of
materials, where the electronic band structures are reproduced by
methods based on wave function theories. Certain details such as the
precise value of Fermi velocities and electronic masses, or the actual
value of the gap in semi conductors may show difference between theory
and experiment, but theoretical results agree at least qualitatively
with measured data. 

However, there are certain compounds where this
classification into metals and insulators fails dramatically. This
happens in particular in systems with open d- and f-shells. There,
band theory predicts metallic behavior because of the open-shell
setting, but in experiments many-not all-of these materials show
actually insulating behavior. This cannot be explained by band theory
and the Pauli principle alone, and a different mechanism has to be
invoked. The bottom line is that these materials do not conduct
current 
because of the strong Coulomb repulsion between the electrons. With
reference to Sir Nevill Mott, who contributed substantially to the
explanation of this effect in the 1930's, these materials are in
general referred to as Mott insulators.

Density-functional theory in a (very small) nutshell
----------------------------------------------------

Density-functional theory tells that the ground state density
determines uniquely all physical properties of a system, independent
of the degree of correlations. Moreover, the theorems of Hohenberg,
Kohn, and Sham state that the full interacting many-body problem can
be replaced by independent electrons moving in an effective
single-particle potential. These leads to the famous Kohn-Sham
equations to be solved in DFT:

.. math::
   H_{KS}\psi_{\nu\mathbf{k}}(\mathbf{r})=\left[-\frac{1}{2m_e}\nabla^2+V_{KS}[\rho](\mathbf{r})\right]\psi_{\nu\mathbf{k}}(\mathbf{r})
   = \varepsilon_{\nu\mathbf{k}}\psi_{\nu\mathbf{k}}(\mathbf{r}).

Without going into details of the Kohn-Sham potential :math:`V_{KS}=V(\mathbf{r})+V_H(\mathbf{r})+V_{xc}(\mathbf{r})`
that is discussed in the literature on DFT, let us just note that the
main result of DFT calculations are the Kohn-Sham energies
:math:`\varepsilon_{\nu\mathbf{k}}` and the Kohn-Sham orbitals :math:`\psi_{\nu\mathbf{k}}(\mathbf{r})`. 
This set of equations is exact, however, the exchange correlation
potential :math:`V_{xc}(\mathbf{r})` is not known explicitly. In
order to do actual calculations, it needs to be approximated in some
way. The local density approximation is one of the most famous
approximations used in this context. This approximation works well for
itinerant systems and semiconductors, but fails completely for
strongly-correlated systems.

From DFT to DMFT
----------------

In order to extend our calculations to strong correlations, we need to
go from a description by bands to a description in terms of
(localized) orbitals: Wannier functions.

In principle, Wannier functions :math:`\chi_{\mu\sigma}(\mathbf{r})`
are nothing else than a Fourier transform of the Bloch basis set from
momentum space into real space,

.. math::
   \chi_{\mu\sigma}(\mathbf{r})=\frac{1}{V}\sum_\mathbf{k} e^{-i\mathbf{k}\mathbf{r}}\sum_\nu U_{\mu\nu}\psi_{\mathbf{k}\nu}^\sigma

where we introduced also the spin degree of freedom :math:`\sigma`. The
unitary matrix :math:`U_{\mu\nu}` is not uniquely defined, but allows for a
certain amount of freedom in the calculation of Wannier function. A
very popular choice is the constraint that the resulting Wannier
functions should be maximally localized in space. Another route,
computationally much lighter and more stable, are projective Wannier
functions. This scheme is used for the Wien2k interface in this
package.

A central quantity in this scheme is the projection operator
:math:`P_{m\nu}(\mathbf{k})`, where :math:`m` is an orbital index and
:math:`\nu` a Bloch band index. 
Its definition and how it is calculated can be found in the original
literature or in the extensive documentation of the
:program:`dmftproj` program shipped with :program:`DFTTools`.

Using projective Wannier functions for DMFT
-------------------------------------------

In this scheme-that is used for the interface to Wien2k-the operators
:math:`P_{m\nu}(\mathbf{k})` are not unitary, since the two dimensions
:math:`m` and :math:`\nu` are not necessarily the same. They
allow, however, to project the local DFT Green function from Bloch band
space into Wannier space,

.. math::
   G^0_{mn}(i\omega) =
   \sum_{\mathbf{k}}\sum_{\nu\nu'}P_{m\nu}(\mathbf{k})G^{DFT}_{\nu\nu'}(\mathbf{k},i\omega)P^*_{\nu'
   n}(\mathbf{k})

with the DFT Green function

.. math::
   G^{DFT}_{\nu\nu'}(\mathbf{k},i\omega) = \frac{1}{i\omega +\mu-\varepsilon_{\nu\mathbf{k}}}\delta_{\nu\nu'}

This non-interacting Green function :math:`G^0_{mn}(i\omega)` defines,
together with the interaction Hamiltonian, the Anderson impurity
model. The DMFT self-consistency cycle can now be formulated as
follows:

#. Take :math:`G^0_{mn}(i\omega)` and the interaction Hamiltonian and
   solve the impurity problem, to get the interacting Greens function
   :math:`G_{mn}(i\omega)` and the self energy
   :math:`\Sigma_{mn}(i\omega)`. For the details of how to do
   this in practice, we refer to the documentation of one of the
   Solver applications, for instance the :ref:`CTHYB solver <triqscthyb:welcome>`.

#. The self energy, written in orbital space, has to be corrected by
   the double counting correction, and upfolded into Bloch band space:

   .. math::
      \Sigma_{\nu\nu'}(\mathbf{k},i\omega) = \sum_{mn}P^*_{\nu
      m}(\mathbf{k}) (\Sigma_{mn}(i\omega) -\Sigma^{DC})P_{n\nu'}(\mathbf{k})

#. Use this :math:`\Sigma_{\nu\nu'}(\mathbf{k},i\omega)` as the DMFT
   approximation to the true self energy in the lattice Dyson
   equation:

   .. math::
      G^{latt}_{\nu\nu'}(\mathbf{k},i\omega)  = \frac{1}{i\omega+\mu
      -\varepsilon_{\nu\mathbf{k}}-\Sigma_{\nu\nu'}(\mathbf{k},i\omega)}

#. Calculate from that the local downfolded Greens function in orbital space:

   .. math::
      G^{loc}_{mn}(i\omega) = \sum_{\mathbf{k}}\sum_{\nu\nu'}P_{m\nu}(\mathbf{k})G^{latt}_{\nu\nu'}(\mathbf{k},i\omega)P^*_{\nu'
      n}(\mathbf{k})

#. Get a new :math:`G^0_{mn}(i\omega)` for the next DMFT iteration
   from

   .. math::
      G^0_{mn}(i\omega) = \left[
      \left(G^{loc}_{mn}(i\omega)\right)^{-1} + \Sigma_{mn}(i\omega)
      \right]^{-1}

   Now go back to step 1 and iterate until convergence.

This is the basic scheme for one-shot DFT+DMFT calculations. Of
course, one has to make sure, that the chemical potential :math:`\mu`
is set such that the electron density is correct. This can be achieved
by adjusting it for the lattice Greens function such that the electron
count is fulfilled.

Full charge self-consistency
----------------------------

The feedback of the electronic correlations to the Kohn-Sham orbitals
is included by the interacting density matrix. With going into the
details, it basically consists of calculating the Kohn-Sham density
:math:`\rho(\mathbf{r})` in the presence of this interacting density
matrix. This new density now defines a new Kohn-Sham
exchange-correlation potential, which in turn leads to new
:math:`\varepsilon_{\nu\mathbf{k}}`,
:math:`\psi_{\nu\mathbf{k}}(\mathbf{r})`, and projectors
:math:`P_{m\nu}(\mathbf{k})`. The update of these
quantities can easily be included in the above
self-consistency cycle, for instance after
step 3, before the local lattice Green
function is downfolded again into orbital space.

How all these calculations can be done in practice with this
:program:`DFTTools` package is subject of the user guide in this documentation.
