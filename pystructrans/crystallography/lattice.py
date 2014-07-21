from pystructrans.general_imports import *
from .matrix_group import MatrixGroup, CUBIC_LAUE_GROUP, SQUARE_GROUP, SQUARE_GROUP_EXT

class Lattice():
    r'''
    This is the base class of lattices. 
    
    Attributes:
    
    .. py:attribute:: __E (2D numpy array) - base matrix of the lattice, each column is a lattice vector
    .. py:attribute:: __N (integer) - the dimension of the vector space where the lattice lives
    
    '''   
    
    def __init__(self, E):
        '''
        Constructor
        '''
        # check if E is a square matrix
        try:
            if not isinstance(E, np.ndarray):
                E = np.array(E)
            if not (len(E.shape) == 2 and E.shape[0]==E.shape[1]):
                raise ValueError('Base is not a square matrix')
            if la.det(E)**2 <= 0.0:
                raise ValueError('Bae matrix is singular')
        except Exception:
            raise ValueError('Base is not a square matrix')

        self.__E = E
        self.__N = len(E)
        self.__LaueGroup = None
        self.__PointGroup = None
        self.__LatticeGroup = None
        self.__SLatticeGroup = None
    
    def getBase(self):
        '''
        Get the base matrix of the lattice.
        
        :return: (*2D numpy array*) - the base matrix 
        '''
        return self.__E
        
    def getDimension(self):
        '''
        Get the dimenstion of the lattice.
        
        :return:  (*integer*) - the dimenstion of the lattice
        '''
        return self.__N

    def getLaueGroup(self):
        '''
        get the Laue group in SO(3) of a 3D lattice

        :return: (*ndarray*) - matrix representation of Laue group
        '''
        if not self.getDimension() in [2, 3]:
            raise AttributeError('Only Laue groups for lattices in 2D and 3D have been implemented.')
        if self.__LaueGroup is None and self.getDimension() == 3:
            self.__LaueGroup = MatrixGroup(np.array([Q for Q in CUBIC_LAUE_GROUP.matrices() if self.inPointGroup(Q)]))
        if self.__LaueGroup is None and self.getDimension() == 2:
            self.__LaueGroup = MatrixGroup(np.array([Q for Q in SQUARE_GROUP if self.inPointGroup(Q)]))
        return self.__LaueGroup
    
    def getPointGroup(self):
        '''
        get the point group in O(3) of a 3D Bravais lattice

        :return: (*ndarray*) - matrix representation of Laue group
        '''
        if not self.getDimension() in [2, 3]:
            raise AttributeError('Only point groups for lattices in 2D and 3D have been implemented.')
        if self.__PointGroup is None and self.__N == 3:
            lg = self.getLaueGroup().matrices()
            self.__PointGroup = MatrixGroup(np.append(lg, -lg, axis=0))
        elif self.__PointGroup is None and self.__N == 2:
            self.__PointGroup = MatrixGroup(np.array([Q for Q in SQUARE_GROUP_EXT if self.inPointGroup(Q)]))
        return self.__PointGroup
    
    def getSpecialLatticeGroup(self):
        '''
        get the special lattice group
        '''
        if not self.getDimension() in [2, 3]:
            raise AttributeError('Only special lattice groups for lattices in 2D and 3D have been implemented.')
        if self.__SLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.getLaueGroup().matrices()], dtype='int')
            return MatrixGroup(mats)
    
    def getLatticeGroup(self):
        '''
        get the lattice group
        '''
        if not self.getDimension() in [2, 3]:
            raise AttributeError('Only lattice groups for lattices in 2D and 3D have been implemented.')
        if self.__SLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.getPointGroup().matrices()], dtype='int')
            return MatrixGroup(mats)

    def __eq__(self, other):
        if not isinstance(other, Lattice):
            return False

        E1 = self.getBase()
        E2 = other.getBase()
        if not self.getDimension() == other.getDimension():
            return False

        A = np.dot(la.inv(E1), E2)
        A_int = np.round(A)
        return np.max(np.abs(A-A_int)) < 1.0E-10

    def __str__(self):
        des_str = '{:d} dimensional lattice:'.format(self.__N)
        for i in range(self.__N):
            des_str = des_str + '\n    e_{:d} = {:s}'.format(i+1, self.__E[:, i])
        return des_str

    def inPointGroup(self, Q):
        '''
        :return: if Q is in the point group
        '''
        # if Q is not an orthogonal matrix, return false
        try:
            Q = Q if isinstance(Q, np.ndarray) else np.array(Q)
        except:
            return False
        if not _in_O3(Q):
            return False
        if len(Q) != self.__N:
            return False
        return self == Lattice(np.dot(Q, self.__E))

def _is_square_mat(M):
    if not isinstance(M, np.ndarray):
        M = np.array(M)
    return (M.ndim == 2) and (M.shape[0] == M.shape[1])

def _in_O3(Q):
    if not _is_square_mat(Q):
        return False
    else:
        return abs(la.det(Q)) == 1.0 and (np.dot(Q, Q.T) == np.eye(len(Q))).all()