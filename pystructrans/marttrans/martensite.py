import numpy as np
from numpy.linalg import inv, det, eig, norm
from math import sqrt
from pystructrans.crystallography import BravaisLattice, CUBIC_LAUE_GROUP
from pystructrans.mat_math import unique_rows
import itertools

class MartensiteError(Exception):
    pass

class Martensite():
    r'''
    To generate all variants and the associated
    Laue group indecies for a given 3x3 U matrix.
    
    Martensite class can be constructed by the Laue group
    of the initial phase and the transformation stretch tensor U matrix.
    
    Methods
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
    
    def calcVariants(self):
        U1 = np.array(self._U)
        ulist = U1.reshape(1,9)
#         print ulist
        laue_idx = [np.array([0])]
        i_exist = [0];
        for i, lg in enumerate(self._Laue):
            utemp = (np.dot(lg, U1.dot(lg.T))).reshape(1, 9)
            flag = True
            
            for j in xrange(len(ulist)):
                du = ulist[j] - utemp
                if norm(du) < 1e-6:
                    flag = False
                    if not i in i_exist:
                        i_exist.append(i);
                        laue_idx[j] = np.append(laue_idx[j], i);            
            
            if flag:
                ulist = np.append(ulist, np.array(utemp), axis = 0)
                i_exist.append(i);
                laue_idx.append(np.array([i]))
                
        laue_idx = np.array(laue_idx)
        
        self._uList = ulist.reshape(len(ulist),3,3);
        self._laueIdx = laue_idx;
            
    def getVariants(self):
        if self._uList == None:
            self.calcVariants()
            
        return self._uList
    
    def getVariant(self,n):
        if self._uList == None:
            self.calcVariants()
        return self._uList[n-1]

    def getLaueIdx(self, *args):
        if self._laueIdx == None:
            self.calcVariants()
        if len(args) == 0:
            return self._laueIdx
        else:
            return self._laueIdx[args[0]-1]
        
    def getLaue_M(self):
        idx = self._laueIdx[0]
        return np.array([self._Laue[i] for i in idx])
    
    def getCor_list(self, E, L):
        if isinstance(E, np.ndarray) and isinstance(E, np.ndarray):
            if E.shape == (3,3) and E.shape == (3,3):
                laue_idx = self.getLaueIdx()
                lg_a = self.getLaue()
                lcor = []
                for i in xrange(len(laue_idx)):
                    ltemp = []
                    for j in xrange(len(laue_idx[0])):
                        ltemp.append((np.dot(inv(E), lg_a[laue_idx[i][j]].dot(E))).dot(L))
                    lcor.append(ltemp)
                return np.array(lcor)
            else:
                raise MartensiteError('Please input lattice base vectors E in R^{3x3} and correspondence matrix L.')
        else:
            raise MartensiteError('Please input lattice base vectors E in R^{3x3} and correspondence matrix L.')
        
        
        
        
        
