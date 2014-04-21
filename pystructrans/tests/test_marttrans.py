import unittest
import numpy as np
try:
    import scipy.linalg as la
except:
    import numpy.linalg as la

from pystructrans.marttrans.twins import *

class TestTwins(unittest.TestCase):

    def test_variants(self):
        with self.assertRaises(TwinsError):
            L = [3]
            getVariants(L)
        with self.assertRaises(TwinsError):
            L = [[1,2,3],[3,4,5],[5,6,7]]
            getVariants(L)
        L = BravaisLattice(6,[2,3])
        vrns = getVariants(L)
        self.assertTrue((vrns[1]==[ 0,  0,  3,  0, -2,  0,  2,  0,  0]).all())