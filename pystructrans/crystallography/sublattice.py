from pystructrans.general_imports import *
import itertools
from operator import mul

def _divisors(n):
    '''
    Find all divisors of an integer n
    
    :param n: the integer to be factorized
    :type n: integer
    :return: vector - the list of all the divisors of n
    '''

    dvrs = [ i + 1
        for i in xrange(int(np.floor(n/2)))
        if n % (i + 1) == 0 ]
    dvrs.append(n)
    return tuple(dvrs)

def hnf_from_diag(diag):
    '''
    Find all the hermite normal forms from the given diagonal elements
    
    :param diag: diagonal elments [d1, ..., dN]
    :type diag: list or 1D array
    :return: [n x N x N] matrix, where n is the number of HNFs with the given diagonal.
    '''
    N = len(diag)
    def set_ij(hnfs, modif):
        newH = []
        for H in hnfs:
            if H[modif[0], modif[1]] == 0:
                M = H.copy()
                M[modif[0], modif[1]] = modif[2]
                newH.append(M)
        hnfs.extend(newH)
        return hnfs
    hnfs = [np.diag(diag)]
    modifs = ((i, j, -d) for i in xrange(N) for j in xrange(i) for d in xrange(1, int(diag[i])))
    reduce(set_ij, modifs, hnfs)
    return np.array(hnfs)

def hnf_from_det(n, N=3):
    '''    
    Find all Hermite Normal Forms with determinant n and dimension N.    
    
    :param n: the integer to be factorized
    :type n: integer
    :return: [9 x N] matrix, where N is the number of HNFs with the given diagonal.
    '''
    if n < 0:
        raise ValueError("negative determinant of HNF")
    # get all the divisors of n
    diags = (d for d in itertools.product(_divisors(n), repeat=N) if reduce(mul, d) == n)
    def ext(x, y):
        if x is None:
            return hnf_from_diag(y)
        return np.append(x, hnf_from_diag(y), axis=0)
    return np.array(reduce(ext, diags, None))

import unittest

class TestSublattice(unittest.TestCase):
    def test_divisor(self):
        self.assertEqual(_divisors(1), (1,))
        self.assertEqual(_divisors(14), (1, 2, 7, 14))
        self.assertEqual(_divisors(18), (1, 2, 3, 6, 9, 18))

    def test_hnf_from_diag(self):
        self.assertIn([[2, 0, 0], [0, 2, 0], [0, 0, 1]], hnf_from_diag([2, 2, 1]))
        hnfs = [
            [[6, 0, 0], [0, 3, 0], [0, 0, 2]],
            [[6, 0, 0], [-1, 3, 0], [0, 0, 2]],
            [[6, 0, 0], [-2, 3, 0], [0, 0, 2]],
            [[6, 0, 0], [0, 3, 0], [-1, 0, 2]],
            [[6, 0, 0], [-1, 3, 0], [-1, 0, 2]],
            [[6, 0, 0], [-2, 3, 0], [-1, 0, 2]],
            [[6, 0, 0], [0, 3, 0], [0, -1, 2]],
            [[6, 0, 0], [-1, 3, 0], [0, -1, 2]],
            [[6, 0, 0], [-2, 3, 0], [0, -1, 2]],
            [[6, 0, 0], [0, 3, 0], [-1, -1, 2]],
            [[6, 0, 0], [-1, 3, 0], [-1, -1, 2]],
            [[6, 0, 0], [-2, 3, 0], [-1, -1, 2]]
        ]
        self.assertEqual(len(hnf_from_diag([6, 3, 2])), len(hnfs))
        for h in hnf_from_diag([6, 3, 2]):
            self.assertIn(h.tolist(), hnfs)

    def test_hnf_from_det(self):
        self.assertListEqual(hnf_from_det(0).tolist(), [[[0, 0, 0], [0, 0, 0], [0, 0, 0]]])
        hnfs = [
            [[1, 0, 0], [0, 1, 0], [0, 0, 2]],
            [[1, 0, 0], [0, 1, 0], [-1, 0, 2]],
            [[1, 0, 0], [0, 1, 0], [0, -1, 2]],
            [[1, 0, 0], [0, 1, 0], [-1, -1, 2]],
            [[1, 0, 0], [0, 2, 0], [0, 0, 1]],
            [[1, 0, 0], [-1, 2, 0], [0, 0, 1]],
            [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
        ]
        self.assertEqual(len(hnf_from_det(2)), len(hnfs))
        for h in hnf_from_det(2):
            self.assertIn(h.tolist(), hnfs)

        self.assertRaises(ValueError, hnf_from_det, -2)