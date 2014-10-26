from __future__ import absolute_import, print_function, division
import unittest
from structrans.marttrans.lat_cor import *
from structrans.general_imports import *
import structrans as pst
from structrans.crystallography import HNFDecomposition as __hnfd
from structrans.marttrans.lat_cor.lat_opt import lat_opt
from structrans.marttrans.lat_cor.dist import Cauchy_dist as dist


def splitHL(Ltot, bravi, pi, bravf, pf):
    lati = pst.BravaisLattice(bravi, pi)
    latf = pst.BravaisLattice(bravf, pf)
    convi = lati.getConventionalTrans()
    convf = latf.getConventionalTrans()
    Lp = convi.dot(Ltot).dot(la.inv(convf))
    return __hnfd(Lp)


class TestLatCorOnNiTi(unittest.TestCase):
    def test_lat_opt_NiTi(self):
        # Kudoh's data
        bravi = 3
        pi = 3.015
        bravf = 12
        pf = [2.898, 4.108, 4.646, 97.78]

        Ltot = np.array([
            [0.5, -0.5, -0.5],
            [1, 1, 0],
            [0.5, -0.5, 1.5]
        ]).T
        H, L = splitHL(Ltot, bravi, pi, bravf, pf)
        lati = pst.BravaisLattice(bravi, pi)
        latf = pst.BravaisLattice(bravf, pf)
        Ei = lati.getbase()
        E1 = Ei.dot(H)
        E2 = latf.getbase()
        E2r = pst.crystallography.lattice_reduction.LLL(E2)
        chi = np.array(np.rint(la.inv(E2r).dot(E2)), dtype='int')
        E2rinv = la.inv(E2r)
        sols = lat_opt(E1, E2r, nsol=3)
        # check dist
        for d, l in zip(sols['dopt'], sols['lopt']):
            self.assertEqual(dist(l, E1, E2rinv), d)

        LGi = pst.Lattice(lati.getConventionalBase()).getspeciallatticegroup().matrices()
        LGf = latf.getspeciallatticegroup().matrices()
        Ci = lati.getConventionalTrans()
        # best solution should be the given one
        has_same_L = False
        lopt = la.inv(Ci).dot(H.dot(sols['lopt'][0]).dot(chi))
        equiv = (mi.dot(lopt).dot(mf) for mi in LGi for mf in LGf)
        for E in equiv:
            if np.array_equal(E, Ltot):
                has_same_L = True
                break
        self.assertTrue(has_same_L)
        best = sols['dopt'][0]

        guess = np.array([0.070472746999312816, 0.60553079720431169, 0.95558248659191025])
        self.assertTrue(np.allclose(sols['dopt'], guess))

        # double the size of the cell
        Ltot = np.array([
            [1, -1, -1],
            [1, 1, 0],
            [1, -1, 3]
        ]).T
        pf[0] *= 2
        pf[2] *=2
        H, L = splitHL(Ltot, bravi, pi, bravf, pf)
        latf = pst.BravaisLattice(bravf, pf)
        E1 = Ei.dot(H)
        E2 = latf.getbase()
        E2r = pst.crystallography.lattice_reduction.LLL(E2)
        sols = lat_opt(E1, E2r, nsol=3)
        self.assertAlmostEqual(guess[0], sols['dopt'][0])

from numpy import array


class TestMartTrans(unittest.TestCase):

    def test_lat_opt(self):
        E1 = np.array([[1, -1, 2], [1, 1, 0], [0, 0, 2]], dtype="float")
        E2 = np.array([[1.414, 0, 0], [0, 1.414, 0], [0, 0, 2]])
        res = lat_opt(E1, E2, disp=0)
        self.assertListEqual(res['lopt'][0].tolist(), [[1, 0, -1], [0, 1, 1], [0, 0, 1]])
        self.assertEqual(1.8251822436107899e-07, res['dopt'][0])

    def test_lat_cor(self):
        res = lat_cor(2, 2, 1, 2, nsol=3, disp=0, miniter=2, maxiter=3)
        target = [{'cor': array([[ 0.,  0.,  1.],
           [ 0.,  1.,  0.],
           [-1.,  0.,  0.]]), 'd': 0.0, 'h': array([[ 1,  0,  0],
           [-1,  2,  0],
           [-1,  0,  2]]), 'l': array([[1, 1, 1],
           [0, 1, 0],
           [0, 0, 1]]), 'lams': array([ 1.,  1.,  1.]), 'U': array([[ 1.,  0.,  0.],
           [ 0.,  1.,  0.],
           [ 0.,  0.,  1.]])}, {'cor': array([[ 0. ,  1. ,  0.5],
           [-0.5,  0.5, -0.5],
           [-0.5, -0.5,  1. ]]), 'd': 1.0, 'h': array([[ 1,  0,  0],
           [ 0,  1,  0],
           [-2,  0,  4]]), 'l': array([[ 0,  2, -1],
           [-1, -1,  0],
           [ 0,  1,  0]]), 'lams': array([ 0.75043998,  0.88161645,  1.51148678]), 'U': array([[ 0.92037457, -0.16109394, -0.04426511],
           [-0.16109394,  1.33109267,  0.27792276],
           [-0.04426511,  0.27792276,  0.89207597]])}, {'cor': array([[ 0.5, -0.5,  0.5],
           [ 0.5,  0.5, -0.5],
           [ 0. ,  1. ,  1. ]]), 'd': 1.25, 'h': array([[1, 0, 0],
           [0, 2, 0],
           [0, 0, 2]]), 'l': array([[ 1, -1, -1],
           [ 0,  1,  0],
           [ 0,  0,  1]]), 'lams': array([ 0.70710678,  1.        ,  1.41421356]), 'U': array([[ 1.20710678,  0.20710678,  0.        ],
           [ 0.20710678,  1.20710678,  0.        ],
           [ 0.        ,  0.        ,  0.70710678]])}, {'cor': array([[ 1. , -0.5,  0. ],
           [ 0. ,  0.5,  1. ],
           [ 0. , -1. ,  0. ]]), 'd': 1.25, 'h': array([[ 1,  0,  0],
           [ 0,  1,  0],
           [-2, -3,  4]]), 'l': array([[ 1,  1,  1],
           [-1,  0,  1],
           [ 0,  0,  1]]), 'lams': array([ 0.70710678,  1.        ,  1.41421356]), 'U': array([[ 0.97140452,  0.02859548, -0.23570226],
           [ 0.02859548,  0.97140452,  0.23570226],
           [-0.23570226,  0.23570226,  1.1785113 ]])}, {'cor': array([[ 0. ,  1. ,  0. ],
           [ 0.5,  0. ,  1. ],
           [ 0.5,  0. , -1. ]]), 'd': 1.25, 'h': array([[ 1,  0,  0],
           [ 0,  1,  0],
           [-3,  0,  4]]), 'l': array([[ 0,  1,  2],
           [ 1, -1,  0],
           [ 0,  1,  1]]), 'lams': array([ 0.70710678,  1.        ,  1.41421356]), 'U': array([[ 1.        ,  0.        ,  0.        ],
           [ 0.        ,  1.06066017,  0.35355339],
           [ 0.        ,  0.35355339,  1.06066017]])}]

        for i in xrange(len(res)):
            sol = res[i]
            tar = target[i]
            self.assertEqual(sol['d'], tar['d'])
            self.assertTrue((np.abs(sol['lams']-tar['lams'])<0.001).all())