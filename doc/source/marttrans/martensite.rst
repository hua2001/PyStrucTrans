Martensite
==========

The class :py:class:`pystructrans.Martensite` implements the concept of martensitic phase transformation.

.. autoclass:: pystructrans.Martensite
    :members: getU, getLaue, setU, setLaue, isreversible, getvariant, getvariants

Examples
--------
Import the class

.. doctest::

    >>> from pystructrans import Martensite

Let's define a cubic to orthorhombic transformation

.. doctest::

    >>> M = Martensite().setU(1.1, 0.05, 0.98)

It will have the cubic Laue group

.. doctest::

    >>> M.getLaue().order()
    24

Note that :py:meth:`pystructrans.Martensite.setU` returns a new Martensite object and leave the invoking one untouched.
Same thing is applied to :py:meth:`pystructrans.Martensite.setLaue`

.. doctest::

    >>> A = Martensite()
    >>> A.setU(1.1, 0.05, 0.98) # doctest: +SKIP
    >>> A.getU()
    Traceback (most recent call last):
        ...
    AttributeError: U has not been initialized

In order to have the concept of `variants`, the Laue group of the low symmetry phase
must be a `proper` subgroup of that of the high symmetry one

.. doctest::

    >>> from pystructrans import Lattice
    >>> U = M.getU()
    >>> lg_low = Lattice(U).getLaueGroup()
    >>> lg_high = M.getLaue()
    >>> lg_high.hassubgroup(lg_low)
    True

Now, we are sure that this transformation is reversible and therefore has variants.

.. doctest::

    >>> M.isreversible()
    True
    >>> for U in M.getvariants():
    ...     print(U.tolist())
    [[1.1, 0.05, 0.0], [0.05, 1.1, 0.0], [0.0, 0.0, 0.98]]
    [[1.1, -0.05, 0.0], [-0.05, 1.1, 0.0], [0.0, 0.0, 0.98]]
    [[0.98, 0.0, 0.0], [0.0, 1.1, -0.05], [0.0, -0.05, 1.1]]
    [[0.98, 0.0, 0.0], [0.0, 1.1, 0.05], [0.0, 0.05, 1.1]]
    [[1.1, 0.0, -0.05], [0.0, 0.98, 0.0], [-0.05, 0.0, 1.1]]
    [[1.1, 0.0, 0.05], [0.0, 0.98, 0.0], [0.05, 0.0, 1.1]]

The number of variants follows `Lagrange's theorem <http://en.wikipedia.org/wiki/Lagrange's_theorem_(group_theory)>`_.

.. doctest::

    >>> lg_high.order()/lg_low.order()
    6.0
    >>> len(M.getvariants())
    6

If the group-subgroup relation does not hold for the transformation, trying to define variants will raise an error

.. doctest::

    >>> M2 = Martensite().setU(1)
    >>> Lattice(M2.getU()).getLaueGroup().order()
    24
    >>> M2.isreversible()
    False
    >>> M2.getvariants()
    Traceback (most recent call last):
        ...
    AttributeError: irreversible martensitic transformations have no variants