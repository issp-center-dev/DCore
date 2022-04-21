Non-interacting limit: ``null``
===============================

This impurity solver returns zero self-energy, which allows to check the non-interacting limit.

How to invoke:

::

    [impurity_solver]
    name = null

In material calculations, we recommend users to use null solver first to check agreement with the DFT calculation.
Note that, since self-energy is zero, one should introduce an artificial broadening to the excitation spectra to make delta functions visible. This can be done by

::

    [tool]
    broadening = 0.001
