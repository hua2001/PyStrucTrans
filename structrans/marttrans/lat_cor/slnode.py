from __future__ import absolute_import
from structrans.general_imports import *

class SLNode:
    """Tree structure of GL(n,Z) group"""

    dim = 3
    # all 12 transvectives in the form of a 6x2 array of matrices
    T1 = np.array([np.eye(3, dtype="int")
                   + np.tensordot(np.eye(3, dtype="int")[i], np.eye(3, dtype="int")[j], axes=0)
                   for i in xrange(3) for j in xrange(3) if i != j])
    T2 = np.array([np.eye(3, dtype="int")
                   - np.tensordot(np.eye(3, dtype="int")[i], np.eye(3, dtype="int")[j], axes=0)
                   for i in xrange(3) for j in xrange(3) if i != j])
    _T = np.append(T1, T2, axis=0)
    onetotwo = [(i, j) for i in xrange(3) for j in xrange(3) if i != j]

    def __init__(self, elem, cache, parent=None, grandpa=None):
        """
        constructed by node element
        """
        self.elem = elem
        self.parent = parent
        self.grandpa = grandpa
        self.elem_dist = SLNode.cache_dist(self.elem, cache)
        self.dist = cache['dist']
        self.cache = cache
        self.children = self.generate_children()
        self.children_dist = None

    def generate_children(self):
        """
        generate children nodes
        """
        nobacktrans = []
        for i in xrange(len(SLNode._T)):
            if SLNode.dim != 3 or (
                not SLNode.gobacktoparent(i, self.parent) and
                not SLNode.gobacktouncle(i, self.parent, self.grandpa) and
                not SLNode.gobacktosibling(i, self.parent)
            ):
                nobacktrans.append(i)
        return ((j, self.elem.dot(SLNode._T[j])) for j in nobacktrans)


    @classmethod
    def cache_dist(cls, elem, cache):
        return cache['dist'](elem)
        # key = tuple(elem.flatten())
        # if not key in cache:
        #     equiv = (m1.dot(elem).dot(m2) for m1 in cache['lg1'] for m2 in cache['lg2'])
        #     d = cache['dist'](elem)
        #     for E in equiv:
        #         cache[tuple(E.flatten())] = d
        # return cache[key]


    @classmethod
    def gobacktoparent(cls, child, parent):
        N = SLNode.dim * (SLNode.dim - 1)
        return False if parent is None else child == (parent + N) % N

    @classmethod
    def gobacktouncle(cls, child, parent, grandpa):
        """
        [Tij, Tjk] = Tik => Tij Tjk = Tik Tjk Tij
        """
        if parent is None or grandpa is None:
            return False
        else:
            N = SLNode.dim * (SLNode.dim - 1)
            ic, kc = SLNode.onetotwo[child % N]
            jp, kp = SLNode.onetotwo[parent % N]
            ig, jg = SLNode.onetotwo[grandpa % N]
            return ic == ig and jg == jp and kp == kc

    @classmethod
    def gobacktosibling(cls, child, parent):
        """
        [Tij Tkl] = 1 =>
        """
        if parent is None:
            return False
        else:
            N = SLNode.dim * (SLNode.dim - 1)
            i, j = SLNode.onetotwo[child % N]
            k, l = SLNode.onetotwo[parent % N]
            return (i, j) > (k, l)

    def calc_neighbor_dist(self):
        """distance values of neighbors"""
        nb = (self.elem.dot(t) for t in SLNode._T)
        return (SLNode.cache_dist(c, self.cache) for c in nb)

    def steepest_des(self):
        """direction of the steepest descendant"""
        dmin = 1E10
        nbmin = None
        for i, nbdist in enumerate(self.calc_neighbor_dist()):
            if nbdist - self.elem_dist < dmin:
                nbmin = i
                dmin = nbdist - self.elem_dist
        return SLNode(self.elem.dot(SLNode._T[nbmin]), self.cache) if dmin < 0 else None

    def is_min(self):
        return all(list(self.calc_neighbor_dist()) >= self.elem_dist)

    def copy(self, that):
        self.elem = that.elem
        self.elem_dist = that.elem_dist
        self.children = that.generate_children()

    def __str__(self):
        return "GLTree - {:s}".format(str(self.elem.flatten()))