from ..general_imports import *
from .binarytree import BinaryTree

class BST(BinaryTree):
    """
    Comparison-based binary search tree.
    Support classes with __cmp__ comparison operator
    """
    def add(self, elem):
        """
        add ``elem`` into the tree

        :param elem: an element comparable to other elements in the tree
        """
        if elem is None:
            pass
        elif elem < self.elem:
            if self.left is None:
                self.left = BST(elem)
            else:
                self.left.add(elem)
        elif elem > self.elem:
            if self.right is None:
                self.right = BST(elem)
            else:
                self.right.add(elem)
        self.updatesize()

    def remove(self, elem):
        """
        remove an element if it is in the tree,
        and return a **new** BST.

        :param elem: the element to be removed
        """
        tree = self.copy()
        if elem == self.elem:
            if self.left is None and self.right is None:
                self = None
            elif self.right is None:
                self = self.left
            elif self.left is None:
                self = self.right
            else:
                self = self.left + self.right
        elif elem < self.elem and self.left is not None:
            self.left = self.left.remove(elem)
        elif elem > self.elem and self.right is not None:
            self.right = self.right.remove(elem)

    def copy(self):
        tree = BST(self.elem)
        if self.left is not None:
            tree.left = self.left.copy()
        if self.right is not None:
            tree.right = self.right.copy()
        return tree

    def __add__(self, other):
        if other is None:
            return self.copy()

        if isinstance(other, BST):
            try:
                self.elem < other.elem
            except:
                raise TypeError("concatenate BST with incomparable elements")

            tree = self.copy()
            tree += other.left
            tree += other.right
            tree.add(other.elem)
            return tree

        else:
            try:
                self.elem < other
            except:
                raise TypeError("add an incomparable element")
            tree = self.copy()
            tree.add(other)
            return tree

    def __contains__(self, item):
        """
        check if item is in the tree
        """
        try:
            self.elem < item
        except:
            raise TypeError("incomparable item")

        if self.elem == item:
            return True
        if self.elem > item and self.left is not None:
            return item in self.left
        elif self.elem < item and self.right is not None:
            return item in self.right
        else:
            return False

    class ReverseorderIter(object):
        """
        Reverseorder traversal iterator.
        """
        def __init__(self, tree):
            self.node = tree
            if tree is not None:
                self.iter = BST.ReverseorderIter(self.node.right)
                self.current = "right"

        def __iter__(self):
            return self

        def next(self):
            return next(self)

        def __next__(self):
            if self.node is None:
                raise StopIteration
            elif self.current is "right":
                try:
                    return next(self.iter)
                except StopIteration:
                    self.current = "left"
                    self.iter = BST.ReverseorderIter(self.node.left)
                    return self.node.elem
            else:
                return next(self.iter)

    def reverseorder(self):
        """

        :return: reverse-order iterator
        """
        return BST.ReverseorderIter(self)

