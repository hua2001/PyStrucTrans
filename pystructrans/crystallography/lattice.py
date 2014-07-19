from pystructrans.general_imports import *

CUBIC_LAUE_GROUP = np.array([# identity
                             [[1,0,0],[0,1,0],[0,0,1]],
                             # two fold rotations
                             [[1,0,0],[0,-1,0],[0,0,-1]],   # [100]
                             [[-1,0,0],[0,1,0],[0,0,-1]],   # [010]
                             [[-1,0,0],[0,-1,0],[0,0,1]],   # [001]
                             [[0,1,0],[1,0,0],[0,0,-1]],    # [110]
                             [[0,-1,0],[-1,0,0],[0,0,-1]],  # [1-10]
                             [[0,0,1],[0,-1,0],[1,0,0]],    # [101]
                             [[0,0,-1],[0,-1,0],[-1,0,0]],  # [10-1]
                             [[-1,0,0],[0,0,1],[0,1,0]],    # [011]
                             [[-1,0,0],[0,0,-1],[0,-1,0]],  # [01-1]
                             # four fold rotations about [100]
                             [[1,0,0],[0,0,-1],[0,1,0]],
                             [[1,0,0],[0,0,1],[0,-1,0]],
                             # four fold rotations about [010]
                             [[0,0,1],[0,1,0],[-1,0,0]],
                             [[0,0,-1],[0,1,0],[1,0,0]],
                             # four fold rotations about [001]
                             [[0,-1,0],[1,0,0],[0,0,1]],
                             [[0,1,0],[-1,0,0],[0,0,1]],
                             # three-fold rotations about [111]
                             [[0,1,0],[0,0,1],[1,0,0]],
                             [[0,0,1],[1,0,0],[0,1,0]],
                             # three-fold rotations about [-111]
                             [[0,0,-1],[-1,0,0],[0,1,0]],
                             [[0,-1,0],[0,0,1],[-1,0,0]],
                             # three-fold rotations about [-1-11]
                             [[0,1,0],[0,0,-1],[-1,0,0]],
                             [[0,0,-1],[1,0,0],[0,-1,0]],
                             # three-fold rotations about [1-11]
                             [[0,0,1],[-1,0,0],[0,-1,0]],
                             [[0,-1,0],[0,0,-1],[1,0,0]]
                             ])

_CUBIC_LAUE_GROUP = np.array([
    # identity
    [[1,0,0],[0,1,0],[0,0,1]],
    # two fold rotations
    [[1,0,0],[0,-1,0],[0,0,-1]],   # [100]
    [[-1,0,0],[0,1,0],[0,0,-1]],   # [010]
    [[-1,0,0],[0,-1,0],[0,0,1]],   # [001]
    [[0,1,0],[1,0,0],[0,0,-1]],    # [110]
    [[0,-1,0],[-1,0,0],[0,0,-1]],  # [1-10]
    [[0,0,1],[0,-1,0],[1,0,0]],    # [101]
    [[0,0,-1],[0,-1,0],[-1,0,0]],  # [10-1]
    [[-1,0,0],[0,0,1],[0,1,0]],    # [011]
    [[-1,0,0],[0,0,-1],[0,-1,0]],  # [01-1]
    # four fold rotations about [100]
    [[1,0,0],[0,0,-1],[0,1,0]],
    [[1,0,0],[0,0,1],[0,-1,0]],
    # four fold rotations about [010]
    [[0,0,1],[0,1,0],[-1,0,0]],
    [[0,0,-1],[0,1,0],[1,0,0]],
    # four fold rotations about [001]
    [[0,-1,0],[1,0,0],[0,0,1]],
    [[0,1,0],[-1,0,0],[0,0,1]],
    # three-fold rotations about [111]
    [[0,1,0],[0,0,1],[1,0,0]],
    [[0,0,1],[1,0,0],[0,1,0]],
    # three-fold rotations about [-111]
    [[0,0,-1],[-1,0,0],[0,1,0]],
    [[0,-1,0],[0,0,1],[-1,0,0]],
    # three-fold rotations about [-1-11]
    [[0,1,0],[0,0,-1],[-1,0,0]],
    [[0,0,-1],[1,0,0],[0,-1,0]],
    # three-fold rotations about [1-11]
    [[0,0,1],[-1,0,0],[0,-1,0]],
    [[0,-1,0],[0,0,-1],[1,0,0]]
])

# this is not really a group, because it contains no identity. Also it does not contain the 2 fold rotation about z-axis.
HEX_LAUE_GROUP = np.array([
                [[np.cos(np.pi/3), np.sin(np.pi/3), 0], [-np.sin(np.pi/3), np.cos(np.pi/3), 0],[ 0, 0, 1.]],
                [[np.cos(-np.pi/3), np.sin(-np.pi/3), 0], [-np.sin(-np.pi/3), np.cos(-np.pi/3), 0],[ 0, 0, 1.]],
                [[np.cos(2*np.pi/3), np.sin(2*np.pi/3), 0], [-np.sin(2*np.pi/3), np.cos(2*np.pi/3), 0],[ 0, 0, 1.]],
                [[np.cos(-2*np.pi/3), np.sin(-2*np.pi/3), 0], [-np.sin(-2*np.pi/3), np.cos(-2*np.pi/3), 0],[ 0, 0, 1.]],
                [[np.cos(np.pi/3), np.sin(np.pi/3), 0], [np.sin(np.pi/3), -np.cos(np.pi/3), 0],[ 0, 0, -1.]],
                [[np.cos(-np.pi/3), np.sin(-np.pi/3), 0], [np.sin(-np.pi/3), -np.cos(-np.pi/3), 0],[ 0, 0, -1.]],
                [[np.cos(2*np.pi/3), np.sin(2*np.pi/3), 0], [np.sin(2*np.pi/3), -np.cos(2*np.pi/3), 0],[ 0, 0, -1.]],
                [[np.cos(-2*np.pi/3), np.sin(-2*np.pi/3), 0], [np.sin(-2*np.pi/3), -np.cos(-2*np.pi/3), 0],[ 0, 0, -1.]]
                ])

_HEX_LAUE_GROUP = np.array([
    [[np.cos(np.pi/3), np.sin(np.pi/3), 0], [-np.sin(np.pi/3), np.cos(np.pi/3), 0],[ 0, 0, 1.]],
    [[np.cos(-np.pi/3), np.sin(-np.pi/3), 0], [-np.sin(-np.pi/3), np.cos(-np.pi/3), 0],[ 0, 0, 1.]],
    [[np.cos(2*np.pi/3), np.sin(2*np.pi/3), 0], [-np.sin(2*np.pi/3), np.cos(2*np.pi/3), 0],[ 0, 0, 1.]],
    [[np.cos(-2*np.pi/3), np.sin(-2*np.pi/3), 0], [-np.sin(-2*np.pi/3), np.cos(-2*np.pi/3), 0],[ 0, 0, 1.]],
    [[np.cos(np.pi/3), np.sin(np.pi/3), 0], [np.sin(np.pi/3), -np.cos(np.pi/3), 0],[ 0, 0, -1.]],
    [[np.cos(-np.pi/3), np.sin(-np.pi/3), 0], [np.sin(-np.pi/3), -np.cos(-np.pi/3), 0],[ 0, 0, -1.]],
    [[np.cos(2*np.pi/3), np.sin(2*np.pi/3), 0], [np.sin(2*np.pi/3), -np.cos(2*np.pi/3), 0],[ 0, 0, -1.]],
    [[np.cos(-2*np.pi/3), np.sin(-2*np.pi/3), 0], [np.sin(-2*np.pi/3), -np.cos(-2*np.pi/3), 0],[ 0, 0, -1.]]
])
_FULL_LAUE_GROUP = np.append(_CUBIC_LAUE_GROUP, _HEX_LAUE_GROUP, axis=0)

# square group
SQUARE_GROUP = np.array([# identity
                         [[1,0],[0,1]],
                         # 90 rotations
                         [[0,-1],[1,0]],
                         [[0,1],[-1,0]],
                         # 180 rotation,
                         [[-1,0],[0,-1]]
                         ])

SQUARE_GROUP_ext = np.array([# mirror along x
                             [[1,0],[0,-1]],
                             # mirror along y 
                             [[-1,0],[0,1]],
                             # mirror along (1,1)
                             [[0,1],[1,0]],
                             # mirror along (-1,1)
                             [[0,-1],[-1,0]]
                             ])

_SQUARE_GROUP = np.array([
    # identity
    [[1,0],[0,1]],
    # 90 rotations
    [[0,-1],[1,0]],
    [[0,1],[-1,0]],
    # 180 rotation,
    [[-1,0],[0,-1]]
])
_SQUARE_GROUP_EXT = np.array([
    # mirror along x
    [[1,0],[0,-1]],
    # mirror along y
    [[-1,0],[0,1]],
    # mirror along (1,1)
    [[0,1],[1,0]],
    # mirror along (-1,1)
    [[0,-1],[-1,0]]
])
_SQUARE_GROUP_EXT = np.append(_SQUARE_GROUP, _SQUARE_GROUP_EXT, axis=0)

class LatticeError(Exception):
    pass

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
            self.__LaueGroup = np.array([Q for Q in _FULL_LAUE_GROUP if self.inPointGroup(Q)])
        if self.__LaueGroup is None and self.getDimension() == 2:
            self.__LaueGroup = np.array([Q for Q in _SQUARE_GROUP if self.inPointGroup(Q)])
        return self.__LaueGroup
    
    def getPointGroup(self):
        '''
        get the point group in O(3) of a 3D Bravais lattice

        :return: (*ndarray*) - matrix representation of Laue group
        '''
        if not self.getDimension() in [2, 3]:
            raise AttributeError('Only point groups for lattices in 2D and 3D have been implemented.')
        if self.__PointGroup is None and self.__N == 3:
            lg = self.getLaueGroup()
            self.__PointGroup = np.append(lg, -lg, axis=0)
        elif self.__PointGroup is None and self.__N == 2:
            self.__PointGroup = np.array([Q for Q in _SQUARE_GROUP_EXT if self.inPointGroup(Q)])
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
            return np.array([(Einv.dot(Q.dot(E))).astype('int') for Q in self.getLaueGroup()])
    
    def getLatticeGroup(self):
        '''
        get the lattice group
        '''
        if not self.getDimension() in [2, 3]:
            raise AttributeError('Only lattice groups for lattices in 2D and 3D have been implemented.')
        if self.__SLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            return np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.getPointGroup()], dtype='int')

    def __eq__(self, other):
        if not isinstance(other, Lattice):
            return False

        E1 = self.getBase()
        E2 = other.getBase()
        if not self.getDimension() == other.getDimension():
            return False

        A = np.dot(la.inv(E1), E2)
        A_int = np.round(A)
        return (np.max(np.abs(A-A_int)) < 1.0E-10)

    def __str__(self):
        des_str = '{:d} dimensional lattice:'.format(self.__N)
        for i in range(self.__N):
            des_str = des_str + '\n    e_{:d} = {:s}'.format(i+1, self.__E[:,i])
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
