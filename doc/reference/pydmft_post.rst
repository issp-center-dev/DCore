pydmft_post
===========

Plot Green's function, self energy, etc.

Usage
-----

::

   $ pydmft_post dmft.ini seedname

Input-file format
-----------------

Example
~~~~~~~

.. literalinclude:: ../tutorial/qe/post.ini

Details
~~~~~~~

* [system] block

============= ================= ===================== ======================================================
Name          Type              Default               Description
============= ================= ===================== ======================================================
t             Float             1.0                   Transfer integral (Nearest neighbor) 
tp            Float             0.0                   Transfer integral (Second nearest)
orbital_model String            single                Chosen from "single", "eg", "t2g", "full-d"
ncor          Integer           1                     Number of correlation shell.
lattice       String            chain                 Chosen from "chain", "square", "cubic", "bethe", and
                                                      "wannier90"
seedname      String            pydmft                Name of the system.
                                                      It should be the same as the seedname of wannier90.
cshell        Integer array     [(0,1),...]           Anguler momentum, and the number of orbitals of each
                                                      correlation shell.
bvec          Float array       [(0.0,0.0,0.0),...]   Reciplocal lattice vectors
nnode         Integer           2                     Number of node for the *k* path
nk            Integer           8                     Number of *k* along each line
knode         Sting Float array [(G,0.0,0.0,0.0),...] The name and the fractional coordinate of each k-node.
============= ================= ===================== ======================================================
