from ..general_imports import *
import unittest
from .. import util
from ..util import Euler, rotation

class TestUtilities(unittest.TestCase):
    def test_pos_def_sym(self):
        U = 1
        self.assertFalse(util._pos_def_sym(U))
        U = [[1, 0, 0], [0, 2, 0], [0 ,0, 3]]
        self.assertFalse(util._pos_def_sym(U))
        U = np.array([[1, 0, 0], [0, 2, 0], [0 ,0, 3]])
        self.assertTrue(util._pos_def_sym(U))
        U = np.array([[1, 0, 0], [0, -2, 0], [0 ,0, 3]])
        self.assertFalse(util._pos_def_sym(U))

    def test_divisor(self):
        self.assertEqual(util.divisors(1), (1,))
        self.assertEqual(util.divisors(14), (1, 2, 7, 14))
        self.assertEqual(util.divisors(18), (1, 2, 3, 6, 9, 18))

    def test_rotations(self):
        t = [90, 90, 90]
        self.assertTrue((abs(Euler(t) - np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])) < 1.0E-15).all())

        t = 60
        z = [0, 0, 1]
        diff = rotation(t, z) - np.array([[0.5, -np.sqrt(3)*0.5, 0], [np.sqrt(3)*0.5, 0.5, 0], [0, 0, 1]])
        self.assertTrue(np.max(np.abs(diff)) < 1.0E-12)

        t = 120
        z = [1, 1, 1]
        v = [1, 0, 0]
        self.assertTrue(np.max(np.abs(rotation(t, z).dot(v) - np.array([0, 1, 0]))) < 1.0E-12)