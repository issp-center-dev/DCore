CT-QMC with segment implementation: ``JO/cthyb-seg``
===========================================================

`JO/CTHYB solver <https://github.com/j-otsuki/cthyb>`_ is an implementation of the hybrization-expansion continuous-time quantum Monte Carlo (CT-HYB) solver with the segment algorithm. Only the density-density interactions can be treared.

A minimal set of parameters are the following.

.. code-block:: ini

    [model]
    density_density = True  # required for segment algorithm

    [system]
    no_tail_fit = True  # recommended for CT-HYB

    [impurity_solver]
    name = JO/cthyb-seg
    exec_path{str} = hybqmc
    MC.n_msr{int} = 100000

There are many optional parameters. User can change them as follows.

.. code-block:: ini

    [impurity_solver]
    control.rand_seed{int} = 0
    control.max_order{int} = 1024
    control.flag_tp{bool} = False
    control.n_tp{int} = 32
    control.n_tp2{int} = 256
    control.flag_vx{bool} = False
    control.n_vx1{int} = 10
    control.n_vx2{int} = 1
    MC.n_warmup{int} = 1000000
    MC.n_bin{int} = 10
    MC.r_add{float} = 0.4
    MC.r_shift{float} = 0.4

This solver support MPI parallel computations. For `np` processes, the number of measurements `n_msr` is divided by `np` and each process performs `n_msr/np` measurements. Computation time is ruduced accordingly.

See the README and samples in `JO/CTHYB solver <https://github.com/j-otsuki/cthyb>`_ for more details.
