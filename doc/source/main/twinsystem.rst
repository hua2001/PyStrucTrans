Twin System
===========

TwinSystem object
-----------------
The class :py:class:`pystructrans.TwinSystem` supports various operations on
a group of transformation stretch tensors, which are the variants of a martensitic phase transformation.

.. autoclass:: pystructrans.TwinSystem
    :members: getUlist, getLaue

    To construct, a list of positive definite symmetric matrices
    and a Laue group are required. Or, use a :py:class:`pystructrans.Martensite` object.

    :param args: If two arguments are given, they must be (``Ulist``, ``Laue``).
                 If only one is given, it could be ``Ulist`` or ``Martensite``.
                 In the case of the former, ``Laue`` is the cubic Laue group.
    :raise ValueError: Illegal construction parameters

    .. automethod:: pystructrans.TwinSystem.gettwintable

        .. testsetup:: *

            from pystructrans import *

        .. doctest::

            >>> ts = TwinSystem(Martensite().setU(0.9, 1.1))
            >>> ts.gettwintable()
            [(0, 1), (0, 2), (1, 2)]

    .. automethod:: pystructrans.TwinSystem.gettwinpairs

        The returned list of TwinPair objects one-to-one correspond
        to the return of :py:meth:`pystructrans.TwinSystem.gettwintable`

        .. doctest::

            >>> tps = ts.gettwinpairs()
            >>> tp = tps[0]

            >>> import numpy as np
            >>> np.array_equal(tp.getUi(), ts.getUlist()[0])
            True
            >>> np.array_equal(tp.getUj(), ts.getUlist()[1])
            True

    .. automethod:: pystructrans.TwinSystem.getconventional

        Returned indices corresponding to the return of
        :py:meth:`pystructrans.TwinSystem.gettwintable` or
        :py:meth:`pystructrans.TwinSystem.gettwinpairs`

        .. doctest::

            >>> ts.getconventional()
            [0, 1, 2]
            >>> [tp.isconventional(ts.getLaue()) for tp in tps]
            [True, True, True]

    .. automethod:: pystructrans.TwinSystem.getcompound

        See the explanation of :py:meth:`pystructrans.TwinSystem.getconventional`

    .. automethod:: pystructrans.TwinSystem.gettypeI

        See the explanation of :py:meth:`pystructrans.TwinSystem.getconventional`

    .. automethod:: pystructrans.TwinSystem.gettypeII

        See the explanation of :py:meth:`pystructrans.TwinSystem.getconventional`
