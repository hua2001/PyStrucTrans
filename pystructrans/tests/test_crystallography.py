from pystructrans.general_imports import *
import unittest

from pystructrans import Lattice, BravaisLattice

class TestLattice(unittest.TestCase):
    def test_construction(self):
        E = np.eye(4)
        L = Lattice(E)
        self.assertEqual(L.getdimension(), 4)
        self.assertListEqual(L.getbase()[:, 0].tolist(), [1, 0, 0, 0])
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
        self.assertRaises(AttributeError, L.getLauegroup)
        self.assertRaises(AttributeError, L.getpointgroup)

        L = BravaisLattice(2, 2)
        Q = [[0., 1., 0.], [1., 0., 0.], [0., 0., -1.]]
        self.assertTrue(L.inpointgroup(Q))
        self.assertEqual(L.getLauegroup().order(), 24)
        self.assertEqual(L.getpointgroup().order(), 48)
        lg = L.getLauegroup().matrices().tolist()
        for Q in L.getpointgroup().matrices():
            self.assertEqual(L, Lattice(Q.dot(L.getbase())))
            self.assertListEqual(Q.T.dot(Q).tolist(), np.eye(3).tolist())
            if Q.tolist() in lg:
                self.assertTrue(la.det(Q) == 1)
            else:
                self.assertTrue(la.det(Q) == -1)

        L = BravaisLattice(12, [2, 3, 4, 87])
        self.assertFalse(L.inpointgroup(Q))
        self.assertEqual(len(L.getLauegroup().matrices()), 2)
        self.assertEqual(len(L.getpointgroup().matrices()), 4)

    def test_lattice_group(self):
        L = Lattice(np.eye(4))
        self.assertRaises(AttributeError, L.getspeciallatticegroup)
        self.assertRaises(AttributeError, L.getlatticegroup)

        L = BravaisLattice(2, 2)
        self.assertEqual(L.getspeciallatticegroup().order(), 24)
        self.assertEqual(L.getlatticegroup().order(), 48)
        slg = L.getspeciallatticegroup().matrices().tolist()
        for M in L.getlatticegroup().matrices():
            self.assertEqual(L, Lattice(L.getbase().dot(M)))
            if M.tolist() in slg:
                self.assertTrue(la.det(M) == 1)
            else:
                self.assertTrue(la.det(M) == -1)

class TestBravaisLattice(unittest.TestCase):
    def test_construction(self):
        L = BravaisLattice(2, 2)
        self.assertEqual(L.getLatID(), 2)
        self.assertEqual(L.getdimension(), 3)
        msg = "3D Bravais Lattice - Cubic F (face centered cubic):    a = 2"
        self.assertEqual(str(L), msg)

        L = BravaisLattice(2, [3, 2], 2)
        self.assertEqual(L.getLatID(), 2)
        self.assertEqual(L.getdimension(), 2)
        msg = "2D Bravais Lattice - Rectangular P:    a = 3,    b = 2"
        self.assertEqual(str(L), msg)

    def test_conversion(self):
        L = BravaisLattice(2, 2)
        T = L.getConventionalTrans()
        self.assertListEqual(T.tolist(), [[1, 1, -1], [-1, 1, 1], [1, -1, 1]])
        E = L.getbase()
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


from ..crystallography.matrix_group import MatrixGroup, CUBIC_LAUE_GROUP


class TestMatrixGroup(unittest.TestCase):
    def test_is_group(self):
        A = []
        self.assertFalse(MatrixGroup.isgroup(A))
        A = [[1, 2, 3]]
        self.assertFalse(MatrixGroup.isgroup(A))
        A = [[[1, 2, 3], [2, 3, 4]]]
        self.assertFalse(MatrixGroup.isgroup(A))
        A = [[[1, 2, 3], [2, 3, 4], [3, 4, 5]]]
        self.assertFalse(MatrixGroup.isgroup(A))
        A = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]
        self.assertTrue(MatrixGroup.isgroup(A) is not False)
        self.assertListEqual(MatrixGroup.isgroup(A).tolist(), [[0]])

    def test_construction(self):
        A = [[[1, 2, 3], [2, 3, 4]]]
        self.assertRaises(ValueError, MatrixGroup, A)
        g = CUBIC_LAUE_GROUP
        t = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
                       17, 18, 19, 20, 21, 22, 23],
                      [1, 0, 3, 2, 15, 14, 12, 13, 9, 8, 11, 10, 6, 7, 5, 4, 20,
                       22, 21, 23, 16, 18, 17, 19],
                      [2, 3, 0, 1, 14, 15, 7, 6, 11, 10, 9, 8, 13, 12, 4, 5, 19,
                       21, 22, 16, 23, 17, 18, 20],
                      [3, 2, 1, 0, 5, 4, 13, 12, 10, 11, 8, 9, 7, 6, 15, 14, 23,
                       18, 17, 20, 19, 22, 21, 16],
                      [4, 14, 15, 5, 0, 3, 19, 23, 22, 18, 21, 17, 16, 20, 1, 2, 12,
                       11, 9, 6, 13, 10, 8, 7],
                      [5, 15, 14, 4, 3, 0, 20, 16, 21, 17, 22, 18, 23, 19, 2, 1, 7,
                       9, 11, 13, 6, 8, 10, 12],
                      [6, 13, 7, 12, 18, 21, 0, 2, 20, 19, 16, 23, 3, 1, 22, 17, 10,
                       15, 4, 9, 8, 5, 14, 11],
                      [7, 12, 6, 13, 22, 17, 2, 0, 23, 16, 19, 20, 1, 3, 18, 21, 9,
                       5, 14, 10, 11, 15, 4, 8],
                      [8, 9, 10, 11, 23, 20, 21, 22, 0, 1, 2, 3, 18, 17, 16, 19, 14,
                       13, 12, 15, 5, 6, 7, 4],
                      [9, 8, 11, 10, 19, 16, 18, 17, 1, 0, 3, 2, 21, 22, 20, 23, 5,
                       7, 6, 4, 14, 12, 13, 15],
                      [10, 11, 8, 9, 16, 19, 22, 21, 3, 2, 1, 0, 17, 18, 23, 20, 15,
                       6, 7, 14, 4, 13, 12, 5],
                      [11, 10, 9, 8, 20, 23, 17, 18, 2, 3, 0, 1, 22, 21, 19, 16, 4,
                       12, 13, 5, 15, 7, 6, 14],
                      [12, 7, 13, 6, 21, 18, 1, 3, 16, 23, 20, 19, 2, 0, 17, 22, 11,
                       4, 15, 8, 9, 14, 5, 10],
                      [13, 6, 12, 7, 17, 22, 3, 1, 19, 20, 23, 16, 0, 2, 21, 18, 8,
                       14, 5, 11, 10, 4, 15, 9],
                      [14, 4, 5, 15, 2, 1, 16, 20, 18, 22, 17, 21, 19, 23, 3, 0, 13,
                       8, 10, 7, 12, 9, 11, 6],
                      [15, 5, 4, 14, 1, 2, 23, 19, 17, 21, 18, 22, 20, 16, 0, 3, 6,
                       10, 8, 12, 7, 11, 9, 13],
                      [16, 23, 20, 19, 10, 9, 14, 5, 12, 7, 13, 6, 15, 4, 11, 8, 17,
                       0, 2, 22, 18, 1, 3, 21],
                      [17, 21, 18, 22, 13, 7, 11, 9, 15, 5, 4, 14, 8, 10, 6, 12, 0,
                       16, 20, 3, 2, 23, 19, 1],
                      [18, 22, 17, 21, 6, 12, 9, 11, 14, 4, 5, 15, 10, 8, 13, 7, 3,
                       23, 19, 0, 1, 16, 20, 2],
                      [19, 20, 23, 16, 9, 10, 4, 15, 13, 6, 12, 7, 5, 14, 8, 11, 21,
                       2, 0, 18, 22, 3, 1, 17],
                      [20, 19, 16, 23, 11, 8, 5, 14, 6, 13, 7, 12, 4, 15, 10, 9, 22,
                       1, 3, 17, 21, 0, 2, 18],
                      [21, 17, 22, 18, 12, 6, 8, 10, 5, 15, 14, 4, 11, 9, 7, 13, 2,
                       19, 23, 1, 0, 20, 16, 3],
                      [22, 18, 21, 17, 7, 13, 10, 8, 4, 14, 15, 5, 9, 11, 12, 6, 1,
                       20, 16, 2, 3, 19, 23, 0],
                      [23, 16, 19, 20, 8, 11, 15, 4, 7, 12, 6, 13, 14, 5, 9, 10, 18,
                       3, 1, 21, 17, 2, 0, 22]])
        logging.debug(g.multable())
        self.assertTrue(np.array_equal(t, g.multable()))
        self.assertEqual(g.order(), 24)



from ..crystallography.sublattice import hnf_from_det, hnf_from_diag, hnfd


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