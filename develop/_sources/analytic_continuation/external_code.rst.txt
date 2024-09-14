How to use an external AC code
===============================

You can use an external code to perform AC of the self-energy from the Matsubara frequency to the real frequency.
The only thing the code needs to do is to write the self-energy in the real frequency as a NumPy binary file, ``post/sigma_w.npz``
from the Matsubara frequency self-energy stored in ``seedname_sigma_iw.npz``.

Input file: ``seedname_sigma_iw.npz``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``seedname_sigma_iw.npz`` is a NumPy binary file, so you can load it with ``numpy.load`` function.

.. code-block:: python3

    import numpy as np
    npz = np.load("./seedname_sigma_iw.npz")

The returned object, ``npz``, is a dictionary. The keys are as follows:

.. csv-table::
    :header: Key, Type, Description
    :widths: 5, 10, 20

    "beta", float, "Inverse temperature"
    "iwn", array of complex, "Matsubara frequency"
    "data#", array of complex, "Self-energy of #-th inequivalent shell"
    "hartree\_fock#", array of complex, "Hartree-Fock term of #-th inequivalent shell"

Here, "#" is 0, 1, 2, ... .

"data#" is a :math:`N_{i\omega} \times N_\text{orb} \times N_\text{orb}` array, where :math:`N_{i\omega}` is the number of Matsubara frequencies, and :math:`N_\text{orb}` is the number of orbitals.

"hartree\_fock#" is a Hartree-Fock term of "#"-th inequivalent shell, :math:`H^\text{f}`.

.. math::
  
    H^\text{f}_{ik} = \sum_{jl} U_{ijkl} \left\langle c^\dagger_j c_l \right\rangle 

The data format is a :math:`N_\text{orb} \times N_\text{orb}` array, where :math:`N_\text{orb}` is the number of orbitals.

Output file: ``post/sigma_w.npz``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The output file, ``post/sigma_w.npz``, is also a NumPy binary file storing one dictionary with the following keys:

.. csv-table::
    :header: Key, Type, Description
    :widths: 5, 10, 20

    "omega", array of real, "frequency"
    "data#", array of complex, "Self-energy of #-th inequivalent shell"

Here, "#" is 0, 1, 2, ... .

"data#" is a :math:`N_{\omega} \times N_\text{orb} \times N_\text{orb}` array, where :math:`N_{\omega}` is the number of frequencies, and :math:`N_\text{orb}` is the number of orbitals.