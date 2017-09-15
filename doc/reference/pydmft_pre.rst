pydmft_pre
==========

Prepare model for pyDMFT.

Usage
-----

::

   $ pydmft_pre input-file

Input-file format
-----------------

* ``t``, 1.0, Transfer integral (Nearest neighbor) 
* ``tp``, 0.0, Transfer integral (Second nearest)
* ``U``, 0.0, On-site Coulomb potential
* ``J``, 0.0, On-site Hund potential
* ``model``, "single", Chosen from "single", "eg", "t2g", "full-d"
* ``nk``, 8, Number of k along each line
* ``nk0``, 0, Number of k (Only wannier90)
* ``nk1``, 0, Number of k (Only wannier90)
* ``nk2``, 0, Number of k (Only wannier90)
* ``ncor``, 1, Number of correlation shell.
* ``lattice``, "chain", Chosen from "chain", "square", "cubic", "bethe", and "wannier90"
* ``nelec``, 1.0, Number of electrons per unit cell.
* ``seedname``, "pydmft", Name of the system.
  It should be the same as the seedname of wannier90.
* ``l-0, l-1, l-2, ...``, 0, Anguler momentum of each shell. 
* ``norb-0, norb-1, norb-2 ...``, 1, Number of orbitals in each shell.
    
.. * ``equiv-0, equiv-1, equiv-2``, differ at each shell, Equivalence of shell.

