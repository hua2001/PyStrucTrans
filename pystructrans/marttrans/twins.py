import numpy as np
from pystructrans.crystallography.bravais_lattice import BravaisLattice
from pystructrans.crystallography.lattice import Lattice, CUBIC_LAUE_GROUP

class TwinsError(Exception):
    pass

def getVariants(L):
    '''
    get all the variants for a 3D lattice
    
    :return: [n x 9] matrix, where n is the number of variants
    '''
    if not isinstance(L, Lattice):
        raise TwinsError('this function only works with 3D Bravais lattice currently.')
    if L.getDimension()!=3:
        raise TwinsError('this function only works with 3D Bravais lattice currently.')
    
    CLG = CUBIC_LAUE_GROUP
    E = L.getBase()
    vrns = np.array([E.reshape(9)])
    
    for i in range(len(CLG)):
        Q = CLG[i]
        E_new = np.dot(Q, E)
        L_new = Lattice(E_new)
        included = False # by default, not included
        for row in vrns:
            E_exist = row.reshape(3,3)
            L_exist = Lattice(E_exist)
            if L_new == L_exist: # if the new lattice is not contained in the variants list
                included = True # flip the flag            
            del L_exist
        if not included:
            vrns = np.vstack((vrns, E_new.reshape(9)))
        del L_new
    return vrns
    

    
    
    