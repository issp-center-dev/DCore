:orphan:

``slater_basis`` parameter
================================================

The basis for the Slater interactions can be specified by ``slater_basis`` parameter. The format is ``slater_basis=[(basis, order), (basis, order),]``. One needs to specify *basis* and *order* for each inequivalent shell as in ``slater_f`` and ``slater_uj`` parameters. The details of *basis* and *order* are described below.

- *basis*

  - ``'cubic'`` : cubic harmonics (default)
  - ``'spherical'`` : spherical harmonics in L-S basis
  - ``'spherical_j'`` : spherical harmonics in J-Jz basis

- *order*

  - The order of basis is specified by index numbers.
  - The basis can be truncted if the number of indices is fewer than `2*l+1`
  - Special word such as 'eg' and 't2g' can be used.
  - If not specified, all bases are used in the default order


Examples
-----------------------------------------------

General examples::

    slater_basis = [('spherical',), ('spherical', 4, 3, 2, 1, 0), ('spherical', 2, 1, 0)]

    slater_basis = [('cubic',), ('cubic', 2, 4), ('cubic', 1, 0, 2)]
    slater_basis = [('cubic',), ('cubic', 'eg'), ('cubic', 1, 0, 2)] # Equivalent

    slater_basis = [('spherical_j',), ('spherical_j',), ('spherical',)]

    slater_basis = 'cubic'  # Default
    slater_basis = [('cubic',), ('cubic',), ('cubic',)] # Equivalent

More specific examples are shown below.

.. contents::
    :local:


cubic harmonics for p orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 3
    interaction = slater_uj
    slater_uj = [(1, 4.0, 0.9)]
    slater_basis = 'cubic'

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(1, 4.0, 0.9)]
    slater_basis(basis) = ['cubic']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 1
        | F_2m = [4.  4.5]
        | basis/sp = ['x' 'y' 'z']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['x' 'y' 'z']
        | basis(dn) = ['x' 'y' 'z']


cubic harmonics for d orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 5
    interaction = slater_uj
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = 'cubic'

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis(basis) = ['cubic']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 2
        | F_2m = [4.         7.73006135 4.86993865]
        | basis/sp = ['xy' 'yz' 'z^2' 'xz' 'x^2-y^2']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['xy' 'yz' 'z^2' 'xz' 'x^2-y^2']
        | basis(dn) = ['xy' 'yz' 'z^2' 'xz' 'x^2-y^2']


cubic harmonics for f orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 7
    interaction = slater_uj
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis = 'cubic'

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis(basis) = ['cubic']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 3
        | F_2m = [ 4.         10.72759974  7.1676259   5.3028777 ]
        | basis/sp = ['x(x^2-3y^2)' 'z(x^2-y^2)' 'xz^2' 'z^3' 'yz^2' 'xyz' 'y(3x^2-y^2)']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['x(x^2-3y^2)' 'z(x^2-y^2)' 'xz^2' 'z^3' 'yz^2' 'xyz' 'y(3x^2-y^2)']
        | basis(dn) = ['x(x^2-3y^2)' 'z(x^2-y^2)' 'xz^2' 'z^3' 'yz^2' 'xyz' 'y(3x^2-y^2)']


d-eg orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 2
    interaction = slater_uj
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = [('cubic', 'eg'),]
    # slater_basis = [('cubic', 2, 4),]  equivalent

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = [('cubic', 'eg')]
    slater_basis(basis) = ['cubic']
    slater_basis(order) = [['eg']]

    Slater interactions
    ish = 0
        | l = 2
        | F_2m = [4.         7.73006135 4.86993865]
        | basis/sp = ['xy' 'yz' 'z^2' 'xz' 'x^2-y^2']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['z^2' 'x^2-y^2']
        | basis(dn) = ['z^2' 'x^2-y^2']


d-t2g orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 3
    interaction = slater_uj
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = [('cubic', 't2g'),]
    # slater_basis = [('cubic', 0, 1, 3),]  equivalent

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = [('cubic', 't2g')]
    slater_basis(basis) = ['cubic']
    slater_basis(order) = [['t2g']]

    Slater interactions
    ish = 0
        | l = 2
        | F_2m = [4.         7.73006135 4.86993865]
        | basis/sp = ['xy' 'yz' 'z^2' 'xz' 'x^2-y^2']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['xy' 'yz' 'xz']
        | basis(dn) = ['xy' 'yz' 'xz']


d-(`xy`, `x^2-y^2`) orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 2
    interaction = slater_uj
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = [('cubic', 0, 4),]

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = [('cubic', 0, 4)]
    slater_basis(basis) = ['cubic']
    slater_basis(order) = [[0, 4]]

    Slater interactions
    ish = 0
        | l = 2
        | F_2m = [4.         7.73006135 4.86993865]
        | basis/sp = ['xy' 'yz' 'z^2' 'xz' 'x^2-y^2']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['xy' 'x^2-y^2']
        | basis(dn) = ['xy' 'x^2-y^2']


spherical harmonics for p orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 3
    interaction = slater_uj
    slater_uj = [(1, 4.0, 0.9)]
    slater_basis = 'spherical'

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(1, 4.0, 0.9)]
    slater_basis(basis) = ['spherical']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 1
        | F_2m = [4.  4.5]
        | basis/sp = ['p-1' 'p+0' 'p+1']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['p-1' 'p+0' 'p+1']
        | basis(dn) = ['p-1' 'p+0' 'p+1']


spherical harmonics for d orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 5
    interaction = slater_uj
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = 'spherical'

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis(basis) = ['spherical']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 2
        | F_2m = [4.         7.73006135 4.86993865]
        | basis/sp = ['d-2' 'd-1' 'd+0' 'd+1' 'd+2']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['d-2' 'd-1' 'd+0' 'd+1' 'd+2']
        | basis(dn) = ['d-2' 'd-1' 'd+0' 'd+1' 'd+2']


spherical harmonics for f orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 7
    interaction = slater_uj
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis = 'spherical'

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis(basis) = ['spherical']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 3
        | F_2m = [ 4.         10.72759974  7.1676259   5.3028777 ]
        | basis/sp = ['f-3' 'f-2' 'f-1' 'f+0' 'f+1' 'f+2' 'f+3']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['f-3' 'f-2' 'f-1' 'f+0' 'f+1' 'f+2' 'f+3']
        | basis(dn) = ['f-3' 'f-2' 'f-1' 'f+0' 'f+1' 'f+2' 'f+3']


j-jz basis for p orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 3
    interaction = slater_uj
    slater_uj = [(1, 4.0, 0.9)]
    slater_basis = 'spherical_j'
    spin_orbit = True

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(1, 4.0, 0.9)]
    slater_basis(basis) = ['spherical_j']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 1
        | F_2m = [4.  4.5]
        | basis/sp = ['p-1' 'p+0' 'p+1']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['j1/2-1/2' 'j3/2-3/2' 'j3/2-1/2']
        | basis(dn) = ['j1/2+1/2' 'j3/2+3/2' 'j3/2+1/2']


j-jz basis for d orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 5
    interaction = slater_uj
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis = 'spherical_j'
    spin_orbit = True

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(2, 4.0, 0.9)]
    slater_basis(basis) = ['spherical_j']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 2
        | F_2m = [4.         7.73006135 4.86993865]
        | basis/sp = ['d-2' 'd-1' 'd+0' 'd+1' 'd+2']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['j3/2-3/2' 'j3/2-1/2' 'j5/2-5/2' 'j5/2-3/2' 'j5/2-1/2']
        | basis(dn) = ['j3/2+3/2' 'j3/2+1/2' 'j5/2+5/2' 'j5/2+3/2' 'j5/2+1/2']


j-jz basis for f orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 7
    interaction = slater_uj
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis = 'spherical_j'
    spin_orbit = True

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis(basis) = ['spherical_j']
    slater_basis(order) = [None]

    Slater interactions
    ish = 0
        | l = 3
        | F_2m = [ 4.         10.72759974  7.1676259   5.3028777 ]
        | basis/sp = ['f-3' 'f-2' 'f-1' 'f+0' 'f+1' 'f+2' 'f+3']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['j5/2-5/2' 'j5/2-3/2' 'j5/2-1/2' 'j7/2-7/2' 'j7/2-5/2' 'j7/2-3/2' 'j7/2-1/2']
        | basis(dn) = ['j5/2+5/2' 'j5/2+3/2' 'j5/2+1/2' 'j7/2+7/2' 'j7/2+5/2' 'j7/2+3/2' 'j7/2+1/2']


`j=5/2` for f orbitals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    [model]
    norb = 3
    interaction = slater_uj
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis = [('spherical_j', 0, 1, 2),]
    spin_orbit = True

Output of ``dcore_pre`` :

::

    Generating U-matrix
    slater_uj = [(3, 4.0, 0.9)]
    slater_basis = [('spherical_j', 0, 1, 2)]
    slater_basis(basis) = ['spherical_j']
    slater_basis(order) = [[0, 1, 2]]

    Slater interactions
    ish = 0
        | l = 3
        | F_2m = [ 4.         10.72759974  7.1676259   5.3028777 ]
        | basis/sp = ['f-3' 'f-2' 'f-1' 'f+0' 'f+1' 'f+2' 'f+3']
        |
        | in SO rep (after transformed, reordered, or truncated)
        | basis(up) = ['j5/2-5/2' 'j5/2-3/2' 'j5/2-1/2']
        | basis(dn) = ['j5/2+5/2' 'j5/2+3/2' 'j5/2+1/2']

