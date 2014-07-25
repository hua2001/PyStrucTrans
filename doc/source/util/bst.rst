Binary search tree
==================

.. autoclass:: pystructrans.BST
    :members: add, preorder, postorder, inorder, reverseorder, levelorder
    :show-inheritance:

    .. py:attribute:: elem

        the element at the node

    .. py:attribute:: left

        the reference to the left subtree. ``None`` if not exists

    .. py:attribute:: right

        the reference to the right subtree. ``None`` if not exists

    .. py:attribute:: size

        number of elements in the tree

Examples
--------

.. testsetup:: *

    from pystructrans import BST

Let's first define a helper function

    >>> def bst2str(node):
    ...     msgleft = bst2str(node.left) if node.left is not None else '#'
    ...     msgright = bst2str(node.right) if node.right is not None else '#'
    ...     return "(" + msgleft + ", " + str(node.elem) + ", " + msgright + ")"

Then we can start test BST.

    >>> t = BST(1)
    >>> print(bst2str(t))
    (#, 1, #)


Operator ``+`` can be used to merge to BSTs if their elements are comparable.

    >>> t2 = BST(2)
    >>> t3 = t + t2

    >>> print(t3.size)
    2
    >>> print(bst2str(t3))
    (#, 1, (#, 2, #))

The original trees remain no change

    >>> print(bst2str(t))
    (#, 1, #)

    >>> print(bst2str(t2))
    (#, 2, #)

We can also directly apply the operator on an element. However, this is not encouraged, because this operation creates a new BST object.

    >>> print(bst2str(t + 3))
    (#, 1, (#, 3, #))

In-place insertion is simply

    >>> t.add(3)
    >>> print(bst2str(t))
    (#, 1, (#, 3, #))

We cannot concatenate two BSTs with incomparable elements.

    >>> t2 = BST("A")
    >>> t + t2
    Traceback (most recent call last):
        ...
    TypeError: concatenate BST with incomparable elements

    >>> t + "B"
    Traceback (most recent call last):
        ...
    TypeError: add an incomparable element

To remove an element from the BST

    >>> t3.remove(2)
    >>> print(bst2str(t3))
    (#, 1, #)

A BST can be iterated (through inorder traversal)

    >>> t = BST(5)
    >>> for i in range(10):
    ...     t.add(i)
    >>> print(bst2str(t))
    ((#, 0, (#, 1, (#, 2, (#, 3, (#, 4, #))))), 5, (#, 6, (#, 7, (#, 8, (#, 9, #)))))

    >>> print([num for num in t])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

With iterator we can do element check

    >>> print(4 in t)
    True

There are also postorder, preorder and levelorder traversal iterators

    >>> print([num for num in t.postorder()])
    [4, 3, 2, 1, 0, 9, 8, 7, 6, 5]

    >>> print([num for num in t.preorder()])
    [5, 0, 1, 2, 3, 4, 6, 7, 8, 9]

    >>> print([num for num in t.levelorder()])
    [5, 0, 6, 1, 7, 2, 8, 3, 9, 4]

For BST, there is an extra traversal order: reverseorder, which will return elements in the reverse order.
In contrast, inorder always return elements in ascendant order.

    >>> import random
    >>> randtree = BST(random.randint(0, 14))
    >>> seq = list(range(15))
    >>> random.shuffle(seq)
    >>> for num in seq:
    ...     randtree.add(num)

Let's check the ordering of this randomly generated BST

    >>> print([num for num in randtree])
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    >>> print([num for num in randtree.reverseorder()])
    [14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]