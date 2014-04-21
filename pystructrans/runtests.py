import unittest

from tests.test_crystallography import TestLattice, TestBravaisLattice, TestSublattice
from tests.test_marttrans import TestTwins

def test():
    '''
    test the package
    '''
    # load the test suites
    suite_lattice = unittest.TestLoader().loadTestsFromTestCase(TestLattice)
    suite_bravais_lattice = unittest.TestLoader().loadTestsFromTestCase(TestBravaisLattice)
    suite_sublattice = unittest.TestLoader().loadTestsFromTestCase(TestSublattice)
    suite_twins = unittest.TestLoader().loadTestsFromTestCase(TestTwins)
    
    # run the tests
    alltest = unittest.TestSuite([suite_lattice, suite_bravais_lattice,suite_sublattice,suite_twins])
    unittest.TextTestRunner(verbosity=1).run(alltest)
    
if __name__ == "__main__":
    test()