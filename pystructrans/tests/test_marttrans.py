from __future__ import print_function

import unittest
import numpy as np
import sys
print(sys.path)
from pystructrans.marttrans.lat_opt import lat_opt

class TestMartTrans(unittest.TestCase):

    def test_lat_opt(self):
        E1 = np.array([[1,-1,2],[1,1,0],[0,0,2]], dtype="float")
        E2 = np.array([[1.414,0,0],[0,1.414,0],[0,0,2]])
          
        # pst.lat_opt(E1, E2, num_sol = 10) # t = 0.330821990967 sec
        # Lo =  [ 1.  0. -1.  0.  1.  1.  0.  0.  1.]
        # Lopt =  [ 0.  1. -1. -1.  0.  1.  0.  0.  1.]
        # min_dist =  1.82518224361e-07
        # dopt = [  1.82518224e-07   1.25030227e+00   1.25120882e+00   2.00151091e+00
        #    3.00060436e+00   3.00241737e+00   3.25271955e+00   4.25271946e+00
        #    5.00271946e+00   5.25241719e+00]
          
        lopt, dopt, _ = lat_opt(E1, E2)
        print(lopt[0])
        self.assertTrue((lopt[0]==[1.,  0., -1.,  0.,  1.,  1.,  0.,  0.,  1.]).all())