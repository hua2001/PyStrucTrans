from pystructrans.general_imports import *
from pystructrans import Lattice, BravaisLattice
import unittest

class TestLattice(unittest.TestCase):
    def test_construction(self):
        E = np.eye(4)
        L = Lattice(E)
        self.assertEqual(L.getDimension(), 4)
        self.assertTrue(all(L.getBase()[:, 0] == np.array([1, 0, 0, 0])))
        E = "adfda"
        self.assertRaises(ValueError, Lattice, E)
        E = [[1, 0, 0], [0, 0, 0], [0, 0, 1]]
        self.assertRaises(ValueError, Lattice, E)

    def test_equality(self):
        E1 = np.random.random((3, 3))
        L1 = Lattice(E1)
        E2 = np.random.random((4, 4))
        L2 = Lattice(E2)

        self.assertFalse(L1 == L2)
        self.assertFalse(L2 == E1)

        A = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
        L2 = Lattice(np.dot(E1, A))
        self.assertTrue(L1 == L2)
        self.assertFalse(L1 == E2)

    def test_point_group(self):
        L = Lattice(np.eye(4))
        self.assertRaises(AttributeError, L.getLaueGroup)
        self.assertRaises(AttributeError, L.getPointGroup)

        L = BravaisLattice(2, 2)
        Q = [[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]]
        self.assertTrue(L.inPointGroup(Q))
        self.assertEqual(len(L.getLaueGroup()), 24)
        self.assertEqual(len(L.getPointGroup()), 48)
        lg = L.getLaueGroup().tolist()
        for Q in L.getPointGroup():
            self.assertEqual(L, Lattice(Q.dot(L.getBase())))
            self.assertListEqual(Q.T.dot(Q).tolist(), np.eye(3).tolist())
            if Q.tolist() in lg:
                self.assertTrue(la.det(Q) == 1)
            else:
                self.assertTrue(la.det(Q) == -1)

        L = BravaisLattice(12, [2, 3, 4, 87])
        self.assertFalse(L.inPointGroup(Q))
        self.assertEqual(len(L.getLaueGroup()), 2)
        self.assertEqual(len(L.getPointGroup()), 4)

    def test_lattice_group(self):
        L = Lattice(np.eye(4))
        self.assertRaises(AttributeError, L.getSpecialLatticeGroup)
        self.assertRaises(AttributeError, L.getLatticeGroup)

        L = BravaisLattice(2, 2)
        self.assertEqual(len(L.getSpecialLatticeGroup()), 24)
        self.assertEqual(len(L.getLatticeGroup()), 48)
        slg = L.getSpecialLatticeGroup().tolist()
        for M in L.getLatticeGroup():
            self.assertEqual(L, Lattice(L.getBase().dot(M)))
            if M.tolist() in slg:
                self.assertTrue(la.det(M) == 1)
            else:
                self.assertTrue(la.det(M) == -1)

class TestBravaisLattice(unittest.TestCase):
    def test_construction(self):
        L = BravaisLattice(2, 2)
        self.assertEqual(L.getLatID(), 2)
        self.assertEqual(L.getDimension(), 3)
        msg = "3D Bravais Lattice - Cubic F (face centered cubic):    a = 2"
        self.assertEqual(str(L), msg)

        L = BravaisLattice(2, [3, 2], 2)
        self.assertEqual(L.getLatID(), 2)
        self.assertEqual(L.getDimension(), 2)
        msg = "2D Bravais Lattice - Rectangular P:    a = 3,    b = 2"
        self.assertEqual(str(L), msg)

    def test_conversion(self):
        L = BravaisLattice(2, 2)
        T = L.getConventionalTrans()
        self.assertTrue((T == np.array([[1, 1, -1], [-1, 1, 1], [1, -1, 1]])).all())
        E = L.getBase()
        self.assertListEqual(E.dot(T).tolist(), [[2, 0, 0], [0, 2, 0], [0, 0, 2]])

        self.assertRaises(ValueError, L.toConventional, 1)
        self.assertRaises(ValueError, L.toConventional, (1, 0, 0, 1))

        idx = [1, 1, 0] # in conventional
        self.assertListEqual(L.toPrimitive(idx).tolist(), [2, 0, 0])
        idx = [(1, 1, 0), (1, 0, 0), (1, 1, 1)] # in conventional
        self.assertListEqual(L.toPrimitive(idx).tolist(), [[2, 0, 0], [1, -1, 1], [1, 1, 1]])

        idx = [2, 0, 0] # in primitive
        self.assertListEqual(L.toConventional(idx).tolist(), [1, 1, 0])
        idx = [[2, 0, 0], [1, -1, 1], [1, 1, 1]]
        self.assertListEqual(L.toConventional(idx).tolist(), [[1, 1, 0], [1, 0, 0], [1, 1, 1]])

# class TestOrientationRelationship(unittest.TestCase):
#     def test_vec_trans(self):
#         # set a lattice correspondence
#         L = np.array([[1,0,-1],[0,1,0],[1,0,1]]).T
#         self.assertTrue(la.det(L)>0)
#         # one index
#         n_A = [1,0,0]
#         n_M = cry.direc_trans(L, n_A)
#         self.assertTrue((n_M == np.array([0.5,0,0.5])).all())
#         # multiple index
#         n_A = [[1,0,0],[1,1,0],[1,1,1]]
#         n_M = cry.direc_trans(L, n_A)
#         self.assertTrue((n_M == np.array([[ 0.5, 0., 0.5], [ 0.5, 1., 0.5], [ 0., 1., 1. ]])).all())
#         
#     def test_Rvec_trans(self):
#         # set a lattice correspondence
#         L = np.array([[1,0,-1],[0,1,0],[1,0,1]]).T
#         # one index
#         n_M = [0.5,0,0.5]
#         n_A = cry.plane_trans(L, n_M)
#         self.assertTrue((n_A == np.array([0,0,1.])).all())
#         # multiple index
#         n_M = [[ 0.5, 0., 0.5], [ 0.5, 1., 0.5], [ 0., 1., 1. ]]
#         n_A = cry.plane_trans(L, n_M)
#         self.assertTrue((n_A == np.array([[0,0,1.],[0,1.,1.],[-1.,1.,1.]])).all())
        