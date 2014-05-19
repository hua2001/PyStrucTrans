import numpy as np
from numpy.linalg import inv, det, eig, norm
from math import sqrt
from pystructrans.mat_math import unique_rows, eigSort
from Martensite import Martensite
import inspect

class TwinSystemError(Exception):
    pass

class TwinSystem():
    r'''
    To generate a twin system for a Martensite class.
    
    Methods
    '''
    
    def __init__(self, *args):
        '''
        Constructor
        '''
        if len(args) == 1:
            if isinstance(args[0], np.ndarray):
                self._uList = args[0]
            elif isinstance(args[0], Martensite):
                self._uList = args[0].getVariants()
        else:
            raise TwinSystemError('Please input a list of variants or a Martensite class.')
         
        self._N = len(self._uList) # Number of variants
        
    
    
    def getTwinpairs(self):
        t_idx = [[0,0]]
        for i in xrange(self._N-1):
            for j in np.arange(i+1, self._N):
                tw_temp = np.dot(inv(self._uList[i]).dot(self._uList[j]), self._uList[j].dot(inv(self._uList[i])))
                e, v = eigSort(tw_temp)
#                 print [i+1, j+1], e[1]
                if abs(e[1] - 1) < 1e-9:
                    t_idx = np.append(t_idx, [[i+1,j+1]], axis=0)
                
        return t_idx[1:]
    
    def isTwinnable(self, i, j):
        
        tw_temp = np.dot(inv(self._uList[i-1]).dot(self._uList[j-1]), self._uList[j-1].dot(inv(self._uList[i-1])))
        e, v = eigSort(tw_temp)
        if abs(e[1] - 1) < 1e-9:
            return True
        else:
            return False













