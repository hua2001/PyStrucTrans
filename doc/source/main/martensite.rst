Martensite
==========

The class :py:class:`pystructrans.Martensite` implements the concept of martensitic phase transformation.

Martensite object
-----------------

.. autoclass:: pystructrans.Martensite
    :members: getU, getLaue, isreversible, getvariant, getvariants

    .. automethod:: pystructrans.Martensite.setU

        It can be defined through one of the following ways

        * :param U: 3 x 3 :py:class:`numpy.ndarray`

            Directly define the transformation stretch tensor

        * :param |alpha|: :py:class:`int` or :py:class:`float`

            Define a transformation stretch tensor in the form of `diag(|alpha|, |alpha|, |alpha|)`

        * :param |alpha|: :py:class:`int` or :py:class:`float`
          :param |beta|: :py:class:`int` or :py:class:`float`

            Define a transformation stretch tensor in the form of `diag(|alpha|, |alpha|, |beta|)`

        * :param |alpha|: :py:class:`int` or :py:class:`float`
          :param |beta|: :py:class:`int` or :py:class:`float`
          :param |gamma|: :py:class:`int` or :py:class:`float`

            Define a transformation stretch tensor in the form

            .. math::

                \begin{pmatrix}
                \alpha & \beta & 0 \\
                \beta & \alpha & 0 \\
                0 & 0 & \gamma
                \end{pmatrix}

        * :param |alpha|: :py:class:`int` or :py:class:`float`
          :param |beta|: :py:class:`int` or :py:class:`float`
          :param |gamma|: :py:class:`int` or :py:class:`float`
          :param |delta|: :py:class:`int` or :py:class:`float`

            Define a transformation stretch tensor in the form

            .. math::

                \begin{pmatrix}
                \alpha & \beta & 0 \\
                \beta & \gamma & 0 \\
                0 & 0 & \delta
                \end{pmatrix}


        :return: a new Martensite object with the current `Laue` and newly assigned `U`
        :rtype: :py:class:`pystructrans.Martensite`
        :raises ValueError:

            if the input does not match any of above, or
            the resulting tensor is not `positive definite` and `symmetric`

        .. |alpha| unicode:: 0x03B1 .. alpha
        .. |beta| unicode:: 0x03B2 .. beta
        .. |gamma| unicode:: 0x03B3 .. gamma
        .. |delta| unicode:: 0x03B4 .. delta


    Note that :py:meth:`pystructrans.Martensite.setU` returns a new Martensite object and leave the invoking one untouched.

    .. testsetup::

        from pystructrans import Martensite

    .. doctest::

        >>> A = Martensite()
        >>> A.setU(1.1, 0.05, 0.98) # doctest: +SKIP
        >>> A.getU()
        Traceback (most recent call last):
            ...
        AttributeError: U has not been initialized

    Same caution should be applied to :py:meth:`pystructrans.Martensite.setLaue`

    .. automethod:: pystructrans.Martensite.setLaue

        There are two ways to set the Laue group for the transformation.

        * :param arg: an integer between 1 and 14
            type of the Bravais lattice of the high symmetry phase

        * :param group: :py:class:`pystructrans.MatrixGroup`
            directly assign the Laue group

        :return: a new Martensite object with the current `U` and newly assigned `Laue`
        :rtype: :py:class:`pystructrans.Martensite`
        :raises ValueError:

            unrecognizable input


Examples
--------
Import the class

.. doctest::

    >>> from structrans import Martensite

    Let's define a cubic to orthorhombic transformation



Let's define a cubic to orthorhombic transformation

.. doctest::

    >>> M = Martensite().setU(1.1, 0.05, 0.98)

It will have the cubic Laue group

.. doctest::

    >>> M.getLaue().order()
    24

    In order to have the concept of
    24

In order to have the concept of `variants`, the Laue group of the low symmetry phase
must be a `proper` subgroup of that of the high symmetry one

.. doctest::

    >>> from structrans import Lattice
        >>> U = M.getU()
        >>> lg_low = Lattice(U).getLauegroup()
        >>> lg_high = M.getLaue()
        >>> lg_high.hassubgroup(lg_low)
        True

        Now, we are sure that this transformation is reversible and therefore has variants.


        >>> U = M.getU()
        >>> lg_low = Lattice(U).getLauegroup()
        >>> lg_high = M.getLaue()
        >>> lg_high.hassubgroup(lg_low)
        True

    Now, we are sure that this transformation is reversible and therefore has variants.


    >>> U = M.getU()
    >>> lg_low = Lattice(U).getLauegroup()
    >>> lg_high = M.getLaue()
    >>> lg_high.hassubgroup(lg_low)
    True

    Now, we are sure that this transformation is reversible and therefore has variants.


    >>> U = M.getU()
    >>> lg_low = Lattice(U).getLauegroup()
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
    >>> Lattice(M2.getU()).getLauegroup().order()
    24
    >>> M2.isreversible()
    False
    >>> M2.getvariants()
    Traceback (most recent call last):
        ...
    AttributeError: irreversible martensitic transformations have no variants