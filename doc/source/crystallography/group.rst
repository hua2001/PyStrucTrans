Matrix representation of groups
===============================

The :py:class:`pystructrans.MatrixGroup` (or :py:class:`pystructrans.crystallography.MatrixGroup`)
class supports operations on a group using matrix representation.

.. autoclass:: pystructrans.MatrixGroup

    It supports the following instance methods:

    .. automethod:: pystructrans.MatrixGroup.matrices
    .. automethod:: pystructrans.MatrixGroup.multable
    .. automethod:: pystructrans.MatrixGroup.order
    .. automethod:: pystructrans.MatrixGroup.isabelian
    .. automethod:: pystructrans.MatrixGroup.hassubgroup

Class methods
-----------------------

.. automethod:: pystructrans.MatrixGroup.isgroup

Examples
--------

The class :py:class:`MatrixGroup` can be directly imported from
:py:mod:`pystructrans`.

.. doctest::

    >>> import numpy as np
    >>> from pystructrans import MatrixGroup

Define a list of matrices and check if it forms a group.

.. doctest::

    >>> B = [
    ... [[1, 0, 0,], [0, 1, 0], [0, 0, 1]],     # identity
    ...  [[np.cos(np.pi / 3), np.sin(np.pi / 3), 0],
    ...  [-np.sin(np.pi / 3), np.cos(np.pi / 3), 0],
    ...  [0, 0, 1]]     # three-fold rotation about [0, 0, 1]
    ... ]
    >>> MatrixGroup.isgroup(B)
    False

If it is a group, we get the multiplication table.

.. doctest::

    >>> A = [[[1, 0, 0,], [0, 1, 0], [0, 0, 1]]]    # only identity
    >>> MatrixGroup.isgroup(A)
    array([[0]])

    >>> A = [
    ...      [[1, 0, 0,], [0, 1, 0], [0, 0, 1]],     # identity
    ...      [[1, 0, 0], [0, -1, 0], [0, 0, -1]]     # two-fold rotation about [1, 0, 0]
    ...     ]
    >>> MatrixGroup.isgroup(A)
    array([[0, 1],
           [1, 0]])

Trying to construct a MatrixGroup with non-group matrix list will raise a :py:class:`ValueError`

.. doctest::

    >>> g = MatrixGroup(B)
    Traceback (most recent call last):
        ...
    ValueError: input does not form a group

Constructing a MatrixGroup with a propitiate input will generate the multiplication table during the construction.

.. doctest::

    >>> g = MatrixGroup(A)
    >>> g.multable()
    array([[0, 1],
           [1, 0]])
    >>> g.isabelian()
    True

Let's try a larger example

.. doctest::

    >>> from pystructrans.crystallography.matrix_group import CUBIC_LAUE_GROUP
    >>> for m in CUBIC_LAUE_GROUP.matrices():
    ...   print(m.tolist())
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
    [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    [[0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]]
    [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]]
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
    [[0, -1, 0], [0, 0, 1], [-1, 0, 0]]
    [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]
    [[0, 0, -1], [1, 0, 0], [0, -1, 0]]
    [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
    [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
    >>> CUBIC_LAUE_GROUP.order()
    24
    >>> CUBIC_LAUE_GROUP.isabelian()
    False
    >>> print(CUBIC_LAUE_GROUP.multable())
    [[ 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23]
     [ 1  0  3  2 15 14 12 13  9  8 11 10  6  7  5  4 20 22 21 23 16 18 17 19]
     [ 2  3  0  1 14 15  7  6 11 10  9  8 13 12  4  5 19 21 22 16 23 17 18 20]
     [ 3  2  1  0  5  4 13 12 10 11  8  9  7  6 15 14 23 18 17 20 19 22 21 16]
     [ 4 14 15  5  0  3 19 23 22 18 21 17 16 20  1  2 12 11  9  6 13 10  8  7]
     [ 5 15 14  4  3  0 20 16 21 17 22 18 23 19  2  1  7  9 11 13  6  8 10 12]
     [ 6 13  7 12 18 21  0  2 20 19 16 23  3  1 22 17 10 15  4  9  8  5 14 11]
     [ 7 12  6 13 22 17  2  0 23 16 19 20  1  3 18 21  9  5 14 10 11 15  4  8]
     [ 8  9 10 11 23 20 21 22  0  1  2  3 18 17 16 19 14 13 12 15  5  6  7  4]
     [ 9  8 11 10 19 16 18 17  1  0  3  2 21 22 20 23  5  7  6  4 14 12 13 15]
     [10 11  8  9 16 19 22 21  3  2  1  0 17 18 23 20 15  6  7 14  4 13 12  5]
     [11 10  9  8 20 23 17 18  2  3  0  1 22 21 19 16  4 12 13  5 15  7  6 14]
     [12  7 13  6 21 18  1  3 16 23 20 19  2  0 17 22 11  4 15  8  9 14  5 10]
     [13  6 12  7 17 22  3  1 19 20 23 16  0  2 21 18  8 14  5 11 10  4 15  9]
     [14  4  5 15  2  1 16 20 18 22 17 21 19 23  3  0 13  8 10  7 12  9 11  6]
     [15  5  4 14  1  2 23 19 17 21 18 22 20 16  0  3  6 10  8 12  7 11  9 13]
     [16 23 20 19 10  9 14  5 12  7 13  6 15  4 11  8 17  0  2 22 18  1  3 21]
     [17 21 18 22 13  7 11  9 15  5  4 14  8 10  6 12  0 16 20  3  2 23 19  1]
     [18 22 17 21  6 12  9 11 14  4  5 15 10  8 13  7  3 23 19  0  1 16 20  2]
     [19 20 23 16  9 10  4 15 13  6 12  7  5 14  8 11 21  2  0 18 22  3  1 17]
     [20 19 16 23 11  8  5 14  6 13  7 12  4 15 10  9 22  1  3 17 21  0  2 18]
     [21 17 22 18 12  6  8 10  5 15 14  4 11  9  7 13  2 19 23  1  0 20 16  3]
     [22 18 21 17  7 13 10  8  4 14 15  5  9 11 12  6  1 20 16  2  3 19 23  0]
     [23 16 19 20  8 11 15  4  7 12  6 13 14  5  9 10 18  3  1 21 17  2  0 22]]

We can check the group-subgroup relations between MatrixGroups
.. doctest::

    >>> [CUBIC_LAUE_GROUP.hassubgroup(g), g.hassubgroup(CUBIC_LAUE_GROUP), CUBIC_LAUE_GROUP.hassubgroup(CUBIC_LAUE_GROUP), g.hassubgroup(g)]
    [True, False, True, True]
