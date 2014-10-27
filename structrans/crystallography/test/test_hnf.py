from __future__ import absolute_import
from structrans.general_imports import *
import unittest

from structrans.crystallography.hnf import *


class TestSublattice(unittest.TestCase):
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
        self.assertListEqual(hnf_from_det(1, 1).tolist(), [[[1]]])

    def test_hnfd(self):
        dim = np.random.randint(1, 5)
        det = np.random.randint(1, 6)
        hnfs = hnf_from_det(det, N=dim)
        for H in hnfs:
            A = np.random.randint(-5, 5, size=(dim, dim))
            while abs(la.det(A)) != 1:
                A = np.random.randint(-5, 5, size=(dim, dim))
            M = H.dot(A)
            H2, L = hnfd(M)
            self.assertTrue(np.array_equal(H2.dot(L), M))
            self.assertTrue(np.array_equal(L, A))
            self.assertTrue(np.array_equal(H, H2))