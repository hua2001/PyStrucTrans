from pystructrans.general_imports import *
import unittest
from pystructrans.marttrans.lat_opt import lat_opt
from pystructrans import lat_cor
import logging
from numpy import array


from ..marttrans import orientation_relationship as OR


class TestOrientationRelationship(unittest.TestCase):
    def test_vec_trans(self):
        # set a lattice correspondence
        L = np.array([[1, 0, -1], [0, 1, 0], [1, 0, 1]]).T

        n_A = [1, 0, 0]
        n_M = OR.direc_trans(L, n_A)
        self.assertListEqual(n_M, [0.5, 0, 0.5])

        n_A = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
        n_M = OR.direc_trans(L, n_A)
        self.assertListEqual(n_M, [[0.5, 0., 0.5], [0.5, 1., 0.5], [0., 1., 1.]])

    def test_plan_trans(self):
        # set a lattice correspondence
        L = np.array([[1, 0, -1], [0, 1, 0], [1, 0, 1]]).T

        p_A = [1, 0, 0]
        self.assertListEqual(OR.plane_trans(L, p_A), [1, 0, 1])

        p_A = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
        p_M = OR.plane_trans(L, p_A)
        self.assertListEqual(p_M, [[1, 0, 1], [1, 1, 1], [0, 1, 2]])

# class TestMartTrans(unittest.TestCase):
#
#     def test_lat_opt(self):
#         E1 = np.array([[1, -1, 2],[1, 1, 0],[0, 0, 2]], dtype="float")
#         E2 = np.array([[1.414, 0, 0],[0, 1.414, 0],[0, 0, 2]])
#         res = lat_opt(E1, E2, disp=0)
#         self.assertListEqual(res['lopt'][0].tolist(), [[1, 0, -1], [0, 1, 1], [0, 0, 1]])
#         self.assertEqual(1.8251822436107899e-07, res['dopt'][0])
#
#     def test_lat_cor(self):
#         res = lat_cor(2, 2, 1, 2, nsol=3, disp=0)
#         target = [{'cor': array([[ 0.,  0.,  1.],
#            [ 0.,  1.,  0.],
#            [-1.,  0.,  0.]]), 'd': 0.0, 'h': array([[ 1,  0,  0],
#            [-1,  2,  0],
#            [-1,  0,  2]]), 'l': array([[1, 1, 1],
#            [0, 1, 0],
#            [0, 0, 1]]), 'lams': array([ 1.,  1.,  1.]), 'U': array([[ 1.,  0.,  0.],
#            [ 0.,  1.,  0.],
#            [ 0.,  0.,  1.]])}, {'cor': array([[ 0. ,  1. ,  0.5],
#            [-0.5,  0.5, -0.5],
#            [-0.5, -0.5,  1. ]]), 'd': 1.0, 'h': array([[ 1,  0,  0],
#            [ 0,  1,  0],
#            [-2,  0,  4]]), 'l': array([[ 0,  2, -1],
#            [-1, -1,  0],
#            [ 0,  1,  0]]), 'lams': array([ 0.75043998,  0.88161645,  1.51148678]), 'U': array([[ 0.92037457, -0.16109394, -0.04426511],
#            [-0.16109394,  1.33109267,  0.27792276],
#            [-0.04426511,  0.27792276,  0.89207597]])}, {'cor': array([[ 0.5, -0.5,  0.5],
#            [ 0.5,  0.5, -0.5],
#            [ 0. ,  1. ,  1. ]]), 'd': 1.25, 'h': array([[1, 0, 0],
#            [0, 2, 0],
#            [0, 0, 2]]), 'l': array([[ 1, -1, -1],
#            [ 0,  1,  0],
#            [ 0,  0,  1]]), 'lams': array([ 0.70710678,  1.        ,  1.41421356]), 'U': array([[ 1.20710678,  0.20710678,  0.        ],
#            [ 0.20710678,  1.20710678,  0.        ],
#            [ 0.        ,  0.        ,  0.70710678]])}, {'cor': array([[ 1. , -0.5,  0. ],
#            [ 0. ,  0.5,  1. ],
#            [ 0. , -1. ,  0. ]]), 'd': 1.25, 'h': array([[ 1,  0,  0],
#            [ 0,  1,  0],
#            [-2, -3,  4]]), 'l': array([[ 1,  1,  1],
#            [-1,  0,  1],
#            [ 0,  0,  1]]), 'lams': array([ 0.70710678,  1.        ,  1.41421356]), 'U': array([[ 0.97140452,  0.02859548, -0.23570226],
#            [ 0.02859548,  0.97140452,  0.23570226],
#            [-0.23570226,  0.23570226,  1.1785113 ]])}, {'cor': array([[ 0. ,  1. ,  0. ],
#            [ 0.5,  0. ,  1. ],
#            [ 0.5,  0. , -1. ]]), 'd': 1.25, 'h': array([[ 1,  0,  0],
#            [ 0,  1,  0],
#            [-3,  0,  4]]), 'l': array([[ 0,  1,  2],
#            [ 1, -1,  0],
#            [ 0,  1,  1]]), 'lams': array([ 0.70710678,  1.        ,  1.41421356]), 'U': array([[ 1.        ,  0.        ,  0.        ],
#            [ 0.        ,  1.06066017,  0.35355339],
#            [ 0.        ,  0.35355339,  1.06066017]])}]
#         for i in xrange(len(res)):
#             sol = res[i]
#             tar = target[i]
#             self.assertEqual(sol['d'], tar['d'])
#             self.assertTrue((np.abs(sol['lams']-tar['lams'])<0.001).all())