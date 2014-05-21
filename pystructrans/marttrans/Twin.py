import numpy as np
from numpy.linalg import inv, det, eig, norm
from math import sqrt
from pystructrans.mat_math import unique_rows, eigSort
from Martensite import Martensite
from pystructrans.crystallography import CUBIC_LAUE_GROUP

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
                if abs(e[1] - 1) < 1e-9:
                    t_idx = np.append(t_idx, [[i+1,j+1]], axis=0)
                
        return t_idx[1:]
    
    def getConventional(self):
        con_idx = []
        for i in xrange(self._N-1):
            for j in np.arange(i+1, self._N):
                tp = TwinPair(self._uList[i], self._uList[j])
                if tp.isConventional():
                    con_idx.append([i+1,j+1])
        return np.array(con_idx)
    def getCompound(self):
        con_idx = []
        for i in xrange(self._N-1):
            for j in np.arange(i+1, self._N):
                tp = TwinPair(self._uList[i], self._uList[j])
                if tp.isCompound():
                    con_idx.append([i+1,j+1])
        return np.array(con_idx)
    def getTypeI(self):
        con_idx = []
        for i in xrange(self._N-1):
            for j in np.arange(i+1, self._N):
                tp = TwinPair(self._uList[i], self._uList[j])
                if tp.isTypeI():
                    con_idx.append([i+1,j+1])
        return np.array(con_idx)
    def getTypeII(self):
        con_idx = []
        for i in xrange(self._N-1):
            for j in np.arange(i+1, self._N):
                tp = TwinPair(self._uList[i], self._uList[j])
                if tp.isTypeII():
                    con_idx.append([i+1,j+1])
        return np.array(con_idx)
        
    
    
        
class TwinPair():
    r'''
    
    It is about all properties associated with a martensite twin.
    A TwinPair class is constructed by two symmetry-related
    martensite variants U_i and U_j.
    
    '''
    def __init__(self, ui, uj):
        
        if isinstance(ui, np.ndarray) and isinstance(uj, np.ndarray):
            if ui.shape == (3,3) and uj.shape == (3,3):
                self._Ui = ui
                self._Uj = uj
            else:
                'Please input a 3x3 positive symmetry matrix as ui, uj.'
        else:
            'Please input a 3x3 positive symmetry matrix as ui, uj.'
        
        
        self._C = np.dot(inv(self._Ui).dot(self._Uj), self._Uj.dot(inv(self._Ui)))
        self._e, self._v = eigSort(self._C)
        self._ax = [[1,0,0],[0,1,0],[0,0,1],
                    [1/sqrt(2), 1/sqrt(2), 0], [1/sqrt(2), -1/sqrt(2), 0],
                    [1/sqrt(2), 0, 1/sqrt(2)], [1/sqrt(2), 0, -1/sqrt(2)],
                    [0, 1/sqrt(2), 1/sqrt(2)], [0, 1/sqrt(2), 1/sqrt(2)]]
        
    def getUi(self):
        return self._Ui
    def getUj(self):
        return self._Uj
    def getAx(self):
        As = self._Ui.dot(self._Ui)
        A_11 = np.dot(self._v[0], As.dot(self._v[0]))
        A_23 = np.dot(self._v[1], As.dot(self._v[2]))
        A_13 = np.dot(self._v[2], As.dot(self._v[0]))
        A_21 = np.dot(self._v[1], As.dot(self._v[0]))
        ehat = []
        for s in [-1, 1]:
            if abs(s*sqrt(self._e[2])*A_23 + A_21) < 1e-6:
                del1 = sqrt(2*(A_11 + s*sqrt(self._e[2])*A_13))
                del3 = s*sqrt(self._e[2])*del1
                _ehat = del1*self._Ui.dot(self._v[0])+del3*self._Ui.dot(self._v[2])
                ehat.append(_ehat/norm(_ehat))
       
        return np.array(ehat)
        
    
    def isTwinnable(self):
        if abs(self._e[1] - 1) < 1e-9:
            return True
        else:
            return False
        
    def isConventional(self):
        twofold = CUBIC_LAUE_GROUP[1:10]
#         print len(twofold)
        flag = False
        for tf in twofold:
            du = (np.dot(tf, self._Ui.dot(tf.T)) - self._Uj).reshape(1,9)
            if norm(du)<1e-6:
                flag = True
        return flag
    
    def isCompound(self):
        n = len(self.getAx())
        if n>1:
            return True
        else:
            return False
        
    def isTypeI(self):
        n = len(self.getAx())
        if n==1:
            return True
        else:
            return False
        
    def isTypeII(self):
        n = len(self.getAx())
        if n==1:
            return True
        else:
            return False

    
    
        
    
                    
        
            
        
        
        
def TwinSolver(U, e, type):
    
    e = np.array(e)
    
    if len(U) - 3 > 1e-9:
        print 'The input should be a 3x3 symmetric matrix!'
        return
      
    if norm((U.T - U).reshape(9)) > 1e-9:
        print 'The input should be a symmetric matrix!'
        return
    if abs(type - 1) < 1e-9:
        n = e
        denominator = np.dot(inv(U).dot(e), inv(U).dot(e))
#         print denominator
        a = 2*(np.dot(inv(U),e)/denominator - U.dot(e))
    else:
        if abs(type - 2) < 1e-9:
            n = 2*(e - np.dot(U, U.dot(e))/np.dot(U.dot(e), U.dot(e)))
            a = U.dot(e)
        else:
            print 'Please input the type of twin system: 1 for Type I twin; 2 for Type II twin'
            return
    
    return np.array([a, n])















