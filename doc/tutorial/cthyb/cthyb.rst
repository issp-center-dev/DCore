.. _howtocthyb:

How to choose parameters for TRIQS/cthyb
========================================

Basic parameters for QMC
------------------------

The accuracy of the TRIQS/cthyb calcuration depends on mainly the following
three parameters.

 * ``n_cycles{int}``

   The number of QMC cycles for comupting any quantities.
   The numerical noise can be reduced by increasing this parameter.
 
   When we use the MPI parallelism, the total number of QMC cycle
   ( ``number-of-processes * n_cycles``) affects the accuracy.
   Therefore,
   **n_cycles in inverse propotion to the number of MPI processes keeps the accuracy**.
 
 * ``n_warmup_cycles{int}``

   The number of QMC cycles for the thermalization before the above main calculation.
   If it is insufficient, the state to be computed may be different from the
   equiribrium state.

   **This parameter is independent from the MPI parallelism.**
 
 * ``length_cycle{int}``

   Each QMC cycles have subcycles.
   The length of this subcycle should be long enough to escape from the auto-correlation.

   **This parameter is independent from the MPI parallelism.**

The computational time is propotional to ``length_cycle*(n_cycles+n_warmup_cycles)``.
The following script is an example to describe the procedure for searching appropriate QMC parameter. 

.. code-block:: bash

   #!/bin/bash

   cat > pre.ini <<EOF
   [model]
   seedname = bethe
   lattice = bethe
   norb = 1
   nelec = 1.0
   t = 1.0
   kanamori = [(4.0, 0.0, 0.0)]

   [system]
   beta = 40.0
   nk = 100

   [impurity_solver]
   name = TRIQS/cthyb
   n_warmup_cycles{int} = 10000
   n_cycles{int} = 10000
   length_cycle{int} = 100

   [control]
   max_step = 1
   EOF

   dcore_pre pre.ini

   n_cycles[0]=3000;   n_warmup_cycles[0]=5000; length_cycle[0]=50
   n_cycles[1]=10000;  n_warmup_cycles[1]=5000; length_cycle[1]=50
   n_cycles[2]=30000;  n_warmup_cycles[2]=5000; length_cycle[2]=50
   n_cycles[3]=100000; n_warmup_cycles[3]=5000; length_cycle[3]=50
   n_cycles[4]=300000; n_warmup_cycles[4]=5000; length_cycle[4]=50
   n_cycles[5]=300000; n_warmup_cycles[5]=5000; length_cycle[5]=100

   for i in `seq 0 5`
   do
       sed -e "/n_cycles/c n_cycles{int} = ${n_cycles[i]}" \
           -e "/n_wamup_cycles/c n_warmup_cycles{int} = ${n_wamup_cycles[i]}" \
           -e "/length_cycle/c length_cycle{int} = ${length_cycle[i]}" \
       pre.ini > dmft.ini
       mpiexec -np 4 dcore dmft.ini
       dcore_check dmft.ini --output bethe.pdf
       mv bethe_sigma.dat cyc${n_cycles[i]}_warm${n_warmup_cycles[i]}_len${length_cycle[i]}.dat
   done

Then, we use GnuPlot as

.. code-block:: gnuplot

   gnuplot> plot [0:10][-0.6:0] \
   "cyc3000_warm5000_len50.dat" u 1:3 w l lw 3, \
   "cyc10000_warm5000_len50.dat" u 1:3 w l lw 3, \
   "cyc30000_warm5000_len50.dat" u 1:3 w l lw 3, \
   "cyc100000_warm5000_len50.dat" u 1:3 w l lw 3, \
   "cyc300000_warm5000_len50.dat" u 1:3 w l lw 3, \
   "cyc300000_warm5000_len100.dat" u 1:3 w l lw 3

and obtain

.. image:: QMCparam.png
   :width: 500
   :align: center

From this plot, we can see that both parameter settings are insufficient and
we have to inclease ``n_cycles`` or ``length_cycle`` or both of them
(In almose cases, ``n_warmup_cycles`` has minor effect).

   
High-frequency tail fit
-----------------------

The self energy computed with QMC becomes noisy at the high frequency region.
This high-frequency tail can be fitted by using the following function:

.. math::

   \Sigma_{\rm tail}(i \omega) \approx \frac{a_1}{\omega} + \frac{a_2}{\omega^2} +
   \frac{a_3}{\omega^3} + \cdots

We will show the procedure to apply this technique.
The original input file (without tail-fit) is as follwos (dmft.ini):

.. code-block:: ini

   [model]
   seedname = bethe
   lattice = bethe
   norb = 1
   nelec = 1.0
   t = 1.0
   kanamori = [(4.0, 0.0, 0.0)]

   [system]
   beta = 40.0
   nk = 100

   [impurity_solver]
   name = TRIQS/cthyb
   n_warmup_cycles{int} = 10000
   n_cycles{int} = 100000
   length_cycle{int} = 50
   
   [control]
   max_step = 1
   
   [tool]
   omega_check = 30.0

Then we execute

.. code-block:: bash

   $ dcore_pre dmft.ini
   $ mpiexec -np 4 dcore dmft.ini
   $ dcore_check dmft.ini --output bethe.pdf

Then, by looking the plot in bethe.pdf,
we chose the energy range where the tail-fit is performed.

.. _tailfit:

.. image:: tailfit.png
   :align: center

In this case, we choose energy range from 6 to 15 [See (a) in the above figure], and
add the following parameters at the ``[system]`` block in the input file
(not the
`solver parameter <https://triqs.ipht.cnrs.fr/applications/cthyb/reference/solve_parameters.html>`_
for the tail fit):

.. code-block:: ini

   [system]
   beta = 40.0
   nk = 100
   perform_tail_fit = True
   fit_max_moment = 2
   fit_min_w = 6.0
   fit_max_w = 15.0

We run ``dcore_check`` again, and we obtain the result as (b) in the above figure.
If we need, we repeat editing the input file and running ``dcore_check`` to tune the input parameters.
After we finish to tuning parameters, we run ``dcore`` again and obtain the result as (c) in the
above figure.


Multi-band system
-----------------

For the multi-band systems, **we have to include the two-pairs insertion/removal move**
in the QMC cycles as

.. code-block:: ini

   [impurity_solver]
   name = TRIQS/cthyb
   move_double{bool} = True

because these moves are disabled in the default setting.

Pade approximation for DOS and spectrum function
------------------------------------------------

DCore currently support the Pade approximation for the numerical continuation
in the calculation of the spectrum function.

The reference grid points for the numerical continuation are specified with the parameter
``omega_pade`` in the ``[tool]`` block.
The good choise of ``omega_pade`` is the maximum frequency
before the self energy becomes noisy.
For example, we can find that ``omega_pade=4.0`` may be a good choice from (a) of the above figure.  
