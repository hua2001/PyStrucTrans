import unittest
import numpy as np
try:
    import scipy.linalg as la
except:
    import numpy.linalg as la

from pystructrans.crystallography.lattice import Lattice, LatticeError
from pystructrans.crystallography.bravais_lattice import BravaisLattice, BravaisLatticeError
import pystructrans.crystallography as cry

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
    
class TestOrientationRelationship(unittest.TestCase):
    def test_vec_trans(self):
        # set a lattice correspondence
        L = np.array([[1,0,-1],[0,1,0],[1,0,1]]).T
        self.assertTrue(la.det(L)>0)
        # one index
        n_A = [1,0,0]
        n_M = cry.direc_trans(L, n_A)
        self.assertTrue((n_M == np.array([0.5,0,0.5])).all())
        # multiple index
        n_A = [[1,0,0],[1,1,0],[1,1,1]]
        n_M = cry.direc_trans(L, n_A)
        self.assertTrue((n_M == np.array([[ 0.5, 0., 0.5], [ 0.5, 1., 0.5], [ 0., 1., 1. ]])).all())
        
    def test_Rvec_trans(self):
        # set a lattice correspondence
        L = np.array([[1,0,-1],[0,1,0],[1,0,1]]).T
        # one index
        n_M = [0.5,0,0.5]
        n_A = cry.plane_trans(L, n_M)
        self.assertTrue((n_A == np.array([0,0,1.])).all())
        # multiple index
        n_M = [[ 0.5, 0., 0.5], [ 0.5, 1., 0.5], [ 0., 1., 1. ]]
        n_A = cry.plane_trans(L, n_M)
        self.assertTrue((n_A == np.array([[0,0,1.],[0,1.,1.],[-1.,1.,1.]])).all())
        