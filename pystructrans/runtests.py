# import unittest
# 
# # from tests.test_crystallography import TestLattice, TestBravaisLattice,TestOrientationRelationship
# from tests.test_marttrans import TestMartTrans
# 
# def test():
#     '''
#     test the package
#     '''
#     # load the test suites
# #     suite_lattice = unittest.TestLoader().loadTestsFromTestCase(TestLattice, TestBravaisLattice)
# #     suite_bravais_lattice = unittest.TestLoader().loadTestsFromTestCase(TestBravaisLattice)
# #     suite_orientation_relationship = unittest.TestLoader().loadTestsFromTestCase(TestOrientationRelationship)
#     suite_marttrans = unittest.TestLoader().loadTestsFromTestCase(TestMartTrans)
#     
#     # run the tests
#     # alltest = unittest.TestSuite([suite_lattice, suite_bravais_lattice,suite_orientation_relationship])
#     alltest = unittest.TestSuite([
#                 suite_marttrans
#                 ])
#     unittest.TextTestRunner(verbosity=1).run(alltest)
#     
# if __name__ == "__main__":
# #     import os
# #     import sys
# #     BASE_DIR = os.path.dirname(os.path.dirname(__file__))
# #     sys.path.insert(0, os.path.join(BASE_DIR))
# #     print sys.path
#     
#     test()
#     