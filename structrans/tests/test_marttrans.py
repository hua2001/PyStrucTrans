# import unittest
#
# from ..general_imports import *
#
# # from ..marttrans.twin import _solvetwin
# from ..marttrans import Martensite
# from ..marttrans import orientation_relationship as OR
#
#
# class TestOrientationRelationship(unittest.TestCase):
#     def test_vec_trans(self):
#         # set a lattice correspondence
#         L = np.array([[1, 0, -1], [0, 1, 0], [1, 0, 1]]).T
#
#         n_A = [1, 0, 0]
#         n_M = OR.direc_trans(L, n_A)
#         self.assertListEqual(n_M, [0.5, 0, 0.5])
#
#         n_A = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
#         n_M = OR.direc_trans(L, n_A)
#         self.assertListEqual(n_M, [[0.5, 0., 0.5], [0.5, 1., 0.5], [0., 1., 1.]])
#
#     def test_plan_trans(self):
#         # set a lattice correspondence
#         L = np.array([[1, 0, -1], [0, 1, 0], [1, 0, 1]]).T
#
#         p_A = [1, 0, 0]
#         self.assertListEqual(OR.plane_trans(L, p_A), [1, 0, 1])
#
#         p_A = [[1, 0, 0], [1, 1, 0], [1, 1, 1]]
#         p_M = OR.plane_trans(L, p_A)
#         self.assertListEqual(p_M, [[1, 0, 1], [1, 1, 1], [0, 1, 2]])
#
# # class TestTwin(unittest.TestCase):
# #     def test_solve_twin(self):
# #         #set U and 2-fold rotation
# #         U = np.array([[2.2,0.1,0.0],
# #                       [0.1,2.1,0.0],
# #                       0.0, 0.0,1.7])
# #         e = np.array([1,1,2])/sqrt(6)
# #         (Q1,a1,n1),(Q2,a2,n2) =Martensite.twin._solvetwin(U, e)
# #         self.assertListEqual(Q1.T.dot(Q1), np.eye(3))
# #         self.assertListEqual(Q2.T.dot(Q2), np.eye(3))
#
#
# class TestMartensite(unittest.TestCase):
#     def test_set_U(self):
#         M = Martensite()
#         U = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
#         M = M.setU(U)
#         self.assertListEqual(M.getU().tolist(), U)
#         self.assertRaises(ValueError, M.setLaue, 0)
#         self.assertRaises(ValueError, M.setLaue, U)
#
#         Us = M.getvariants()
#         self.assertEqual(len(Us), 6)
#         idx_list = M.getLaueidx()
#         lg = M.getLaue().matrices()
#         for i, idx in enumerate(idx_list):
#             for j in idx:
#                 Q = lg[j]
#                 V = Us[i]
#                 self.assertListEqual(Q.dot(U).dot(Q.T).tolist(), V.tolist())