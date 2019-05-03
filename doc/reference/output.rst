.. _outputformat:

Output-file format
==================

``dcore_pre``
~~~~~~~~~~~~~

-   *seedname*.h5

    .. todo:: Update dcore_pre.txt

    It has two groups, ``dft_input`` and ``Dcore``.
    See `DFTTools <https://triqs.ipht.cnrs.fr/applications/dft_tools/reference/h5structure.html>`_ for the details of the data structure in ``dft_input`` group.
    The data included in ``Dcore`` group is list below:

    .. include:: dcore_pre.txt

``dcore``
~~~~~~~~~

-   *seedname*.out.h5

    .. todo:: Update dcore_out.txt

    All data are stored in ``dmft_out`` group.
    The following list summarizes data structure in the ``dmft_out`` group:

    .. include:: dcore_out.txt

    See
    `TRIQS manual <https://triqs.ipht.cnrs.fr/1.x/reference/gfs/py/full.html#hdf5>`_
    for the data structure of the Green's function and self-energy.

-   solver_dependent_output

    All solver-dependent output are stored in the working directory such as **work/imp_shell0_iter1**.

``dcore_check``
~~~~~~~~~~~~~~~

-   **sigma.dat**

    The local self energy at the final step.

    ::

       # Local self energy at imaginary frequency
       # [Column] Data
       # [1] Frequency
       # [2] Re(Sigma_{shell=0, spin=up, 0, 0})
       # [3] Im(Sigma_{shell=0, spin=up, 0, 0})
       # [4] Re(Sigma_{shell=0, spin=down, 0, 0})
       # [5] Im(Sigma_{shell=0, spin=down, 0, 0})
       -157.001093 0.994751 0.006358 0.994751 0.006358
       -156.844013 0.994751 0.006365 0.994751 0.006365
       -156.686934 0.994751 0.006371 0.994751 0.006371
       :

-   **iter_mu.dat**

    The chemical potential as a function of the iteration number.

    ::

        1 0.0000000000e+00
        2 1.3397270680e-01
        3 4.5709763936e-01
        4 6.2124557444e-01
        5 6.3750111472e-01
        6 6.7087331832e-01
        7 6.9841342338e-01

-   **iter_sigma.dat**

    The average self-energy as a function of the iteration number.

    ::

        1 0.6674359500130874 0.6674359500130874
        2 0.5244344115963627 0.5244344115963644
        3 0.3210245784417104 0.3210245784417107
        4 0.18521138795848135 0.1852113879584814
        5 0.1463266063561347 0.1463266063561347
        6 0.12071320930740709 0.12071320930740712
        7 0.10214335854450268 0.1021433585445027

``dcore_post``
~~~~~~~~~~~~~~

-   *seedname*\_dos.dat

    .. todo:: paste data

-   *seedname*\_dos0.dat

    .. todo:: paste data

-   *seedname*\_akw.dat

    .. todo:: paste data

-   *seedname*\_akw0.dat

    .. todo:: paste data
