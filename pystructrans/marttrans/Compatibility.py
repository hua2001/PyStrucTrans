from ..general_imports import *

from numpy.linalg import inv, det, eig, norm
from math import sqrt

from ..crystallography import BravaisLattice, CUBIC_LAUE_GROUP
from .. import mat_math as _math
from .martensite import Martensite

class CompatibilityError(Exception):
    pass

def isCompatible(A, B):
    r'''
    check rank1 connection between A and B
    A - B = rank1
    '''
    if isinstance(A, np.ndarray) and isinstance(B, np.ndarray):
        C = np.dot(inv(B.T).dot(A.T), A.dot(inv(B)))
        e, v = _math.eigSort(C)
        if abs(e[1])-1<1e-5:
            return True
        else:
            return False
    else:
        print('The input should be two arrays.')
        return False


def AM_Solver(A):
    
    kappa = 1
    
    if len(A) - 3 > 1e-9:
        raise CompatibilityError('The input should be a positive symmetric matrix!')
      
    if norm((A.T - A).reshape(9)) > 1e-9:
        raise CompatibilityError('The input should be a positive symmetric matrix!')
    
    e = eig(A)
    e = [np.append(e[0][i],e[1][:,i]).real for i in xrange(3)]
    e = np.array(sorted(e, key=lambda x:x[0]))
    eval = e[:,0]
    evec = np.array([[1.,0,0],[0,1.,0],[0,0,1.]]).dot(e[:,1:4])
    
    
    c = sqrt(eval[2]-eval[0])
    c1 = sqrt(1-eval[0])
    c3 = sqrt(eval[2]-1)
    
    if c < 1e-6:
        print('solution is b = e, m = - 2 e where |e| = 1.')
        return
    else:
        if abs(eval[1] - 1) <1e-4:
            m1 = ((sqrt(eval[2])-sqrt(eval[0]))/c)*(-c1*evec[0] + kappa*c3*evec[2])
            m2 = ((sqrt(eval[2])-sqrt(eval[0]))/c)*(-c1*evec[0] - kappa*c3*evec[2])
            rho1 = norm(m1)
            rho2 = norm(m2)
            m1 = m1/rho1
            m2 = m2/rho2
            b1 = rho1*((sqrt(eval[2])*c1/c)*evec[0] + kappa*(sqrt(eval[0])*c3/c)*evec[2])
            b2 = rho2*((sqrt(eval[2])*c1/c)*evec[0] - kappa*(sqrt(eval[0])*c3/c)*evec[2])
        
            return np.array([[b1,m1], [b2,m2]])
        

class Microstructure():
    r'''
    
    Given a positive symmetric 3x3 U matrix and the point group of austenite,
    it generates microstructure between austenite and single/twinned martensite.
    
    Anaglous to class Martensite, it uses the same constructor.
    
    '''
    def __init__(self, *args):
        '''
        Constructor
        '''
        if len(args) > 0:
            self._U = args[0]
        if len(args) > 1:
            self._Laue = args[1]
        else:
            self._Laue = CUBIC_LAUE_GROUP
        self._uList = None
        self._laueIdx = None
        
    def getU(self):
        return self._U
    
    def setU(self, *args):
        
        if len(args) == 1:
            if isinstance(args[0], np.ndarray):
                self._U = args[0]
            else:
                self._U = np.diag([args[0],args[0],args[0]])
        elif len(args) == 2:
            self._U = np.diag([args[0], args[0], args[1]])
        elif len(args) == 3:
            self._U = np.array([[args[0],args[1],0],[args[1],args[0],0],[0,0,args[2]]])
        elif len(args) == 4:
            self._U = np.array([[args[0], args[1], 0],[args[1], args[2], 0],[0, 0, args[3]]])
        else:
            raise MartensiteError('You can input a 3x3 U matrix or a vector [a_i] where i = 1, 2, 3, 4 to construct a variant.')
                                    
            
    def getLaue(self):
        return self._Laue
    def setLaue(self, arg):
        latParam = [2,2,2,[2,3],[2,120],[1.41, 2], [1.41, 2],[1.41, 1.42, 2], [1.41, 1.42, 2],[1.41, 1.42, 2], [1.41, 1.42, 2], [1.41, 1.42, 2, 89], [1.41, 1.42, 2, 89], [1.41, 1.42, 2, 89, 88, 92]]
        
        if isinstance(arg, int) and 0 < arg <= 14:            
            Lat_init = BravaisLattice(arg, latParam[arg-1])
            self._Laue = Lat_init.getLaueGroup()
        else:
            raise MartensiteError("Please input Bravais lattice ID that is an integer between 1 and 14. ")
    
