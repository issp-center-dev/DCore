.. _howtocthyb:

How to choose parameters for TRIQS/cthyb
========================================


Basic parameters for QMC
------------------------

The accuracy of the TRIQS/cthyb calcuration depends on mainly the following
three parameters.

 * n_cycles{int}

   The number of QMC cycles for comupting any quantities.
   The numerical noise can be reduced by increasing this parameter.
 
 * n_warmup_cycles{int}

   The number of QMC cycles for the thermalization before the above main calculation.
   If it is insufficient, the state to be computed may be different from the
   equiribrium state.

 * length_cycle{int}

   Each QMC cycles have subcycles.
   If
