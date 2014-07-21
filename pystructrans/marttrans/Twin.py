from ..general_imports import *
import numpy as np
from numpy.linalg import inv, det, eig, norm
from math import sqrt
from ..mat_math import unique_rows, eigSort
from ..crystallography import CUBIC_LAUE_GROUP
from .. import mat_math as _math
from .martensite import Martensite
from . import Compatibility as cp

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
    
        
    
    
class TwinPairError(Exception):
    pass     
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
        if norm((self._Ui - self._Uj).reshape(1,9))>1e-6:
            if abs(self._e[1] - 1) < 1e-9:
                return True
            else:
                return False
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
    def satisfyCofI(self):
        if abs(self.XI()-1)<1e-4:
            return True
        else:
            return False
    def satisfyCofII(self):
        if abs(self.XII()-1)<1e-4:
            return True
        else:
            return False
    
    def getTwinParam(self):
        r'''
        It returns twinning parameters for
        the TwinPair class if it's twinnable.
        The twinning parameters include
        (a, n) where
        q Uj - Ui = a \otimes n
        Therefore there are two solutions usually.
        
        '''
        t_ax = self.getAx()
        tParam = [np.zeros(6)]
        if len(t_ax)>0:
            for ehat in t_ax:
                t1 = TwinSolver(self._Ui, ehat, 1).flatten()
                t2 = TwinSolver(self._Ui, ehat, 2).flatten()
                t1IsNew = True
                t2IsNew = not isSameAN(t1, t2)
                for t in tParam:
                    if isSameAN(t,t1):
                        t1IsNew = False
                    if isSameAN(t,t2):
                        t2IsNew = False
                if t1IsNew:
                    tParam.append(t1)
                if t2IsNew:
                    tParam.append(t2)
            
            return np.array(tParam[1:]).reshape(len(tParam)-1, 2,3)
        else:
            print('The pair of variants is not twinnable!')
            
    def QI(self):
        r'''
        return the rotation matrix of Uj relative to Ui for the first solution
        '''
        anlist = self.getTwinParam()
        if len(anlist) == 2:
            return (self._Ui + np.outer(anlist[0][0], anlist[0][1])).dot(inv(self._Uj))
        else:
            raise TwinPairError('Not a valid twin pair. Please check if it is twinnable.')
        
    def QII(self):
        r'''
        return the rotation matrix of Uj relative to Ui for the second solution
        '''
        anlist = self.getTwinParam()
        if len(anlist) == 2:
            return (self._Ui + np.outer(anlist[1][0], anlist[1][1])).dot(inv(self._Uj))
        else:
            raise TwinPairError('Not a valid twin pair. Please check if it is twinnable.')
    def XI(self):
        if self.isTypeI():
            t_ax = self.getAx()[0]
            return norm(inv(self._Ui).dot(t_ax))
    def XII(self):
        if self.isTypeII():
            t_ax = self.getAx()[0]
            return norm(self._Ui.dot(t_ax))
        
    def volume_frac(self, *args):
        anlist = self.getTwinParam()
        if len(args)>0:
            if args[0] == 1:
                a = anlist[args[0]-1][0]
                n = anlist[args[0]-1][1]
                cofU = _math.cofactor(self._Ui.dot(self._Ui) - np.eye(3))
                minus_alpha = 2*np.dot(self._Ui.dot(a), cofU.dot(n))    
                beta = det(self._Ui.dot(self._Ui)-np.eye(3))+minus_alpha/4
                if abs(minus_alpha)<1e-5:
                    return np.array([0,1])
                elif beta/minus_alpha > 1e-6:
                    f = 0.5+sqrt(beta/minus_alpha)
                    return np.array([f, 1-f])
                else:
                    return None
                
            elif args[0] == 2:
                a = anlist[args[0]-1][0]
                n = anlist[args[0]-1][1]
                cofU = _math.cofactor(self._Ui.dot(self._Ui) - np.eye(3))
                minus_alpha = 2*np.dot(self._Ui.dot(a), cofU.dot(n))    
                beta = det(self._Ui.dot(self._Ui)-np.eye(3))+minus_alpha/4
                if abs(minus_alpha)<1e-5:
                    return np.array([0,1])
                elif beta/minus_alpha > 1e-6:
                    f = 0.5+sqrt(beta/minus_alpha)
                    return np.array([f, 1-f])
                else:
                    return None
            else:
                raise TwinPairError('Please input n = 1 or 2, which means type 1 or 2 twin systems. Or leave it as blank.')
        else:
            vol = []
            for an in anlist:
                a = an[0]
                n = an[1]
                cofU = _math.cofactor(self._Ui.dot(self._Ui) - np.eye(3))
                minus_alpha = 2*np.dot(self._Ui.dot(a), cofU.dot(n))    
                beta = det(self._Ui.dot(self._Ui)-np.eye(3))+minus_alpha/4
#                 print("minus_alpha=%.5f"%minus_alpha)
#                 print("g(0)g(1/2)=%.5f"%(beta*(-minus_alpha/4+beta)))
#                 print("g(0)=%.5f"%(-minus_alpha/4+beta))
#                 print("g(1/2)=%.5f"%(beta))
                if abs(minus_alpha)<1e-5:
                    vol.append([0,1])
                elif beta/minus_alpha > 1e-6:
                    f = 0.5+sqrt(beta/minus_alpha)
                    vol.append([f, 1-f])
                else:
                    pass
            if len(vol)>0:
                return np.array(vol)
            else:
                return None
            
    def habit_planes(self,n):
        anlist = self.getTwinParam()
        habit = []
        rotation = []
        if n==1 or n==2:
            f = self.volume_frac(n)
            if not f==None:
                a = anlist[n-1][0]
                n = anlist[n-1][1]
                Uf = [self._Ui + f[0]*np.outer(a, n),
                  self._Ui + f[1]*np.outer(a, n)]
                for uf in Uf:
                    if cp.isCompatible(uf.T.dot(uf), np.eye(3)):
                        [[b1,m1], [b2,m2]] = cp.AM_Solver(uf.T.dot(uf))
                        R = [(np.eye(3) + np.outer(b1, m1)).dot(inv(uf)), 
                             (np.eye(3) + np.outer(b2, m2)).dot(inv(uf))]
                        habit.append([[b1, m1], [b2, m2]])
                        rotation.append(R)
        else:
            raise TwinPairError('Please input n = 1 or 2.')
        if len(habit)>0:
            return np.array(habit), np.array(rotation)
        else:
            return None
                
            
def isSameAN(t1, t2):
    M1 = np.outer(t1[0:3], t1[3:6])
    M2 = np.outer(t2[0:3], t2[3:6])   
    return np.max(np.abs(M1-M2)) < 1E-3  
        
        
def TwinSolver(U, e, type):
    
    e = np.array(e)
    
    if len(U) - 3 > 1e-9:
        print('The input should be a 3x3 symmetric matrix!')
        return
      
    if norm((U.T - U).reshape(9)) > 1e-9:
        print('The input should be a symmetric matrix!')
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
            print('Please input the type of twin system: 1 for Type I twin; 2 for Type II twin')
            return
    rho = norm(n)
    a = np.round(rho*a, 6)
    n = np.round(n/rho, 6)
    return np.array([a, n])















