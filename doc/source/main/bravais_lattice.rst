Bravais lattice
===============

A lattice defined in the class :py:class:`pystructrans.Lattice` is also called a Bravais lattice.
In 2D and 3D, these lattices have been systematically classified.
Using this classification, a 2D or 3D lattice can be simply defined
by its *Bravais lattice type* and *lattice parameters*.
5 2D Bravais lattices and 14 3D Bravais lattices are defined as following two tables.

- *2D Bravais lattices:*

   .. tabularcolumns:: L|L|L

   ======  ==================== ============================
   **ID**  **Name**             **Lattice Parameters**
   ======  ==================== ============================
   1       Square P             *a*
   2       Rectangular P        *a*, *b*
   3       Rectangular I        *a*, *b*
   4       Hexagonal P          *a*
   5       Oblique P            *a*, *b*, *γ*
   ======  ==================== ============================

- *3D Bravais lattices:*

   .. tabularcolumns:: L|L|L

   ======  ==================== ============================
   **ID**  **Name**             **Lattice Parameters**
   ======  ==================== ============================
   1       Cubic P              *a*
   2       Cubic F              *a*
   3       Cubic I              *a*
   4       Hexagonal P          *a*, *c*
   5       Rhombohedral P       *a*, *α*
   6       Tetragonal P         *a*, *c*
   7       Tetragonal I         *a*, *c*
   8       Orthorhombic P       *a*, *b*, *c*
   9       Orthorhombic C       *a*, *b*, *c*
   10      Orthorhombic F       *a*, *b*, *c*
   11      Orthorhombic I       *a*, *b*, *c*
   12      Monoclinic P         *a*, *b*, *c*, *β*
   13      Monoclinic C         *a*, *b*, *c*, *β*
   14      Triclinic  P         *a*, *b*, *c*, *α*, *β*, *γ*
   ======  ==================== ============================

:py:class:`pystructrans.BravaisLattice` is a subclass of :py:class:`pystructrans.Lattice`.
It is defined based on this classification.

.. note:: Mathematically, all :py:class:`pystructrans.Lattice` objects are Bravais lattices.
          Our naming is only for the convenience of implementation.

BravaisLattice object
---------------------

.. testsetup:: *

    from pystructrans import *

.. autoclass:: pystructrans.BravaisLattice
    :members:

    .. doctest::

        >>> Lat = BravaisLattice(12, [2, 3., 4, 89.75])
        >>> print(Lat)
        3D Bravais Lattice - Monoclinic P:    a = 2,    b = 3,    c = 4,    beta = 89.75

        >>> Lat2 = BravaisLattice(4, 2, N=2)
        >>> print(Lat2)
        2D Bravais Lattice - Hexagonal P:    a = 2