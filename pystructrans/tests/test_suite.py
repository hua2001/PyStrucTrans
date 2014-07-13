from __future__ import print_function, division

import unittest
import pystructrans

class TestLattice(unittest.TestCase):

    def test_lattice(self):
        lat = pystructrans.BravaisLattice(2, 3)
        print(lat)
        pass