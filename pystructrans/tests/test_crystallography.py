import unittest
import numpy as np
try:
    import scipy.linalg as la
except:
    import numpy.linalg as la


from pystructrans.crystallography.lattice import Lattice, LatticeError
from pystructrans.crystallography.bravais_lattice import BravaisLattice, BravaisLatticeError
from pystructrans.crystallography.sublattice import *

class TestLattice(unittest.TestCase):

    def test_lattice(self):
        E = np.eye(4)
        L = Lattice(E)
        self.assertEqual(L.getDimension(), 4)
        self.assertTrue(min(L.getBase()[:,0] == np.array([1, 0, 0, 0])))
        del L

    def test_lattice_eq(self):
        E1 = np.random.random((3,3))
        L1 = Lattice(E1)
        E2 = np.random.random((4,4))
        L2 = Lattice(E2)
        
        with self.assertRaises(LatticeError):
            L1 == E2
        with self.assertRaises(LatticeError):
            L1 == L2
        del L2
        A = np.array([[0, 1, 0],[1, 0 ,0],[0, 0, 1]])
        L2 = Lattice(np.dot(E1, A))
        self.assertTrue(L1==L2)
        del L1
        del L2
        
class TestBravaisLattice(unittest.TestCase):

    def test_bravais_lattice(self):
        L = BravaisLattice(2,2)
        self.assertEqual(L.getLatID(), 2)
        del L

    def test_point_group(self):
        L = BravaisLattice(2,2)
        E = L.getBase()
        LG = L.getLaueGroup()
        LGZ = L.getSpecialLatticeGroup()
        Q = LG[-1]
        A = LGZ[-1]
        self.assertTrue((np.dot(Q, E)==np.dot(E, A)).all())
        del L

class TestSublattice(unittest.TestCase):

    def test_divisors(self):
        pass
        