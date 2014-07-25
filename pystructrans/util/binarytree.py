from ..general_imports import *

class BinaryTree(object):
    """
    Comparison-based binary search tree.
    Support classes with __cmp__ comparison operator
    """
    def __init__(self, elem):
        self.elem = elem
        self.left = None
        self.right = None
        self.size = 1

    def add(self, elem):
        """
        add ``elem`` into the tree

        :param elem: an element comparable to other elements in the tree
        """
        if self.left is None:
            self.left = BinaryTree(elem)
        if self.right is None:
            self.right = BinaryTree(elem)
        if self.left.size < self.right.size:
            self.left.add(elem)
        else:
            self.right.add(elem)
        self.updatesize()

    def updatesize(self):
        # update size
        self.size = 1
        self.size += 0 if self.left is None else self.left.size
        self.size += 0 if self.right is None else self.right.size

    def copy(self):
        tree = BinaryTree(self.elem)
        if self.left is not None:
            tree.left = self.left.copy()
        if self.right is not None:
            tree.right = self.right.copy()
        return tree

    # TODO: implement iterators by stack/queue
    class InorderIter(object):
        """
        Inorder traversal iterator.
        """
        def __init__(self, tree):
            self.node = tree
            if tree is not None:
                self.iter = BinaryTree.InorderIter(self.node.left)
                self.current = "left"

        def __iter__(self):
            return self

        def next(self):
            return next(self)

        def __next__(self):
            if self.node is None:
                raise StopIteration
            else:
                if self.current is "left":
                    try:
                        return next(self.iter)
                    except StopIteration:
                        elem = self.node.elem
                        self.current = "right"
                        self.iter = BinaryTree.InorderIter(self.node.right)
                        return elem
                else:
                    return next(self.iter)

    class PreorderIter(object):
        """
        Preorder traversal iterator. This is the default iterator of binary trees.
        """
        def __init__(self, tree):
            self.node = tree
            self.current = "elem"

        def __iter__(self):
            return self

        def next(self):
            return next(self)

        def __next__(self):
            if self.node is None:
                raise StopIteration
            else:
                if self.current is "elem":
                    elem = self.node.elem
                    self.current = "left"
                    self.iter = BinaryTree.PreorderIter(self.node.left)
                    return elem
                elif self.current is "left":
                    try:
                        return next(self.iter)
                    except StopIteration:
                        self.current = "right"
                        self.iter = BinaryTree.PreorderIter(self.node.right)
                return next(self.iter)

    class PostorderIter(object):
        """
        Postorder traversal iterator.
        """
        def __init__(self, tree):
            self.node = tree
            if tree is not None:
                self.iter = BinaryTree.PostorderIter(self.node.left)
                self.current = "left"

        def __iter__(self):
            return self

        def next(self):
            return next(self)

        def __next__(self):
            if self.node is None or self.current is "elem":
                raise StopIteration
            elif self.current is "left":
                try:
                    return next(self.iter)
                except StopIteration:
                    self.current = "right"
                    self.iter = BinaryTree.PostorderIter(self.node.right)
                    return self.next()
            elif self.current is "right":
                try:
                    return next(self.iter)
                except StopIteration:
                    elem = self.node.elem
                    self.current = "elem"
                    return elem

    class LevelIter(object):
        """
        level order traversal iterator
        """
        def __init__(self, tree):
            self.lev = 0
            self.queue = [tree]

        def __iter__(self):
            return self

        def next(self):
            return next(self)

        def __next__(self):
            if len(self.queue) == 0:
                raise StopIteration
            else:
                head = self.queue.pop(0)
                if head.left is not None:
                    self.queue.append(head.left)
                if head.right is not None:
                    self.queue.append(head.right)
                return head.elem

    def __iter__(self):
        """
        inorder traversal iterator (because it is sorted)
        """
        return self.inorder()

    def preorder(self):
        """

        :return: preorder iterator
        """
        return BinaryTree.PreorderIter(self)

    def postorder(self):
        """

        :return: postorder iterator
        """
        return BinaryTree.PostorderIter(self)

    def inorder(self):
        """

        :return: inorder iterator; This is the default iterator of binary trees.
        """
        return BinaryTree.InorderIter(self)

    def levelorder(self):
        """

        :return: by level iterator
        """
        return BinaryTree.LevelIter(self)

    def isbalanced(self):
        """

        :return: if the tree is balanced
        """
