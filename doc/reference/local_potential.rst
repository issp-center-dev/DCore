:orphan:

``local_potential_*`` parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An arbitrary local potential can be implemented using parameters ``local_potential_*``.
The format looks like

::

    [model]
    local_potential_matrix = {0: 'pot0.txt', 1: 'pot1.txt'}
    local_potential_factor = 0.01

Here, ``local_potential_matrix`` describes, in the python dictionary format, a set of the inequivalent shell index *ish* and the filename which defines the local potential matrix.
The parameter ``local_potential_factor`` defines a prefactor to the potential matrix.

For example, the Zeeman term along z-axis for S=1/2 is represented by

::

    $ cat pot0.txt
    # spin orb1 orb2  Re Im
    0 0 0   0.5 0.
    1 0 0  -0.5 0.

and the magnetic field is specified by ``local_potential_factor``.
