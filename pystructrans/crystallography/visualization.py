import numpy as np
from bravais_lattice import BravaisLattice

class visualError(Exception):
    pass

def vertex(o, E, n=3):
    dim = np.array(E).shape
    if n==3:
        if dim[0]==3 and dim[1]==3 and len(o)==3:
            return np.array([o,
                              o+E[0],o+E[1],o+E[2],
                              o+E[0]+E[1],o+E[1]+E[2],o+E[0]+E[2],
                              o+E[0]+E[1]+E[2]])
        else:
            print 'E has to be a 3x3 matrix.'
            return None
    elif n==2:
        if dim[0]==2 and dim[1]==2 and len(o)==2:
            return np.array([o,
                             o+E[0], o+E[1],
                             o+E[0]+E[1]])
        else:
            print 'E has to be a 2x2 matrix.'
            return None
    else:
        raise visualError('The dimension n has to be 2 or 3!')
            
            
class UnitCell():
    r'''
    unit_cell class gives the vertices
    it can be constructed by a lattice (2D or 3D).
    
    '''
    def __init__(self, origin, lat_ID, lat_Params, N=3):
        r'''
        Constructed by a Bravais lattice
        '''
        self._lattice = BravaisLattice(lat_ID, lat_Params, N=3)
        self._o = np.array(origin)
        self._dim = N
    def setOrigin(self, new_origin):
        self._o = new_origin
        
    
    def getPrimitive(self):
        E = self._lattice.getBase()
        return vertex(self._o, E, self._dim)
    
    def getConventional(self):
        C = self._lattice.getConventionalBase()
        return vertex(self._o, C, self._dim)
        
        
        
    
    