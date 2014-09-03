Twin System
===========
.. testsetup:: *

    from pystructrans import *

TwinSystem object
-----------------
The class :py:class:`pystructrans.TwinSystem` supports various operations on
a group of transformation stretch tensors, which are the variants of a martensitic phase transformation.

.. autoclass:: pystructrans.TwinSystem
    :members: 

    To construct, a list of positive definite symmetric matrices
    and a Laue group are required. Or, use a :py:class:`pystructrans.Martensite` object.

    :param args: If two arguments are given, they must be (``Ulist``, ``Laue``).
                 If only one is given, it could be ``Ulist`` or ``Martensite``.
                 In the case of the former, ``Laue`` is the cubic Laue group.
    :raise TypeError: Illegal construction parameters

    