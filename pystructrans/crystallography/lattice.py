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
        self.__setBase(E);
    
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
        
    def __setBase(self, E):
        '''
        Set the base matrix :math:`\mathbf E` of the lattice.
        
        :param E: base matrix, each column is a lattice vector.
        :type E: 2D numpy array, square
        :return: none
        '''
        # check if E is a square matrix 
        if not (len(E.shape) == 2 and E.shape[0]==E.shape[1]):
            raise LatticeError('The base matrix must be a square matrix.')
        if la.det(E)**2 < 1e-16:
            raise LatticeError('The base matrix cannot be singular.')
        self.__E = E;
        self.__N = E.shape[0];
    
    def transformBase(self, A):
        r'''
        Transform the base matrix of the lattice.
        
        .. warning:: The argument must be a unimodular integral matrix.
        
        :param A: base transformation matrix
        :type A: 2D numpy array, integers
        :return: none
        '''
        # check if A is a square matrix 
        if not (len(A.shape) == 2 and A.shape[0]==A.shape[1]):
            raise LatticeError('Base transformation matrix must be a square matrix.')
        # check if A is an integeral matrix
        if not A.dtype == 'int':
            raise LatticeError('Base transformation matrix must be an integeral matrix.')
        # check if A is unimodular
        if not (np.abs(la.det(A))-1)**2 < 1.e-10:
            raise LatticeError('Base transformation matrix must be a unimodular matrix.')
        # check if the dimension of A is the same as that of E
        if not A.shape[0] == self.__E.shape[0]:
            raise LatticeError('Base transformation matrix must have the same size of the base matrix.')
        
        self.__setBase(np.dot(self.__E, A))    
    
    def __eq__(self, other):
        E1 = self.getBase()
        
        try:
            E2 = other.getBase()
        except:
            raise LatticeError('the two objects are not comparable')
        
        if not self.getDimension()==other.getDimension():
            raise LatticeError('the two lattices have to have the same dimension before comparison')
        
        A = np.dot(la.inv(E1), E2)
        A_int = np.round(A)
        return (np.max(np.abs(A-A_int)) < 1e-12)
        
    def __str__(self):
        des_str = '{:d} dimensional lattice:'.format(self.__N)
        for i in range(self.__N):
            des_str = des_str + '\n    e_{:d} = {:s}'.format(i+1, self.__E[:,i])
        return des_str   
    
    def getLaueGroup(self):
        '''
        get the Laue group in SO(3) of a 3D lattice
        '''
        if not self.getDimension() in [2, 3]:
            raise LatticeError('only 2D or 3D lattices have Laue groups')
        
        if self.getDimension() == 3:
            CLG = CUBIC_LAUE_GROUP.copy()
            HLG = HEX_LAUE_GROUP.copy()
            CLG = np.append(CLG, HLG, axis=0)
            
            LG = [np.eye(3)]
            for i in range(1, len(CLG)):
                Q = CLG[i]
                if inPointGroup(Q, self):
                    LG = np.append(LG, [Q], axis = 0)
        
        if self.getDimension() == 2:
            SG = SQUARE_GROUP
            
            LG = [np.eye(2)]
            for i in range(1, len(SG)):
                Q = SG[i]
                if inPointGroup(Q, self):
                    LG = np.append(LG, [Q], axis = 0)                
        return np.array(LG)
    
    def getPointGroup(self):
        '''
        get the point group in O(3) of a 3D Bravais lattice
        '''
        if self.getDimension()==3:
            LG = self.getLaueGroup()
            return np.append(LG, -LG, axis=0)   
        else:
            SG = SQUARE_GROUP.copy()
            SG2 = SQUARE_GROUP_ext.copy()
            SLG = np.append(SG, SG2, axis=0)
            LG = [np.eye(2)]
            for i in range(1, len(SLG)):
                Q = SLG[i]
                if inPointGroup(Q, self):
                    LG = np.append(LG, [Q], axis = 0)
            return np.array(LG)
    
    def getSpecialLatticeGroup(self):
        '''
        get the special lattice group
        '''
        SPG = self.getLaueGroup()
        SLG = [np.int_(np.eye(self.getDimension()))]
        E = self.getBase()
        for i in range(1, len(SPG)):
            R = SPG[i]
            Et = np.dot(R,E)
            A = np.dot(la.inv(E), Et)
            SLG = np.append(SLG,[np.int_(np.rint(A))], axis=0)
        return np.array(SLG)
    
    def getLatticeGroup(self):
        '''
        get the lattice group
        '''
        PG = self.getPointGroup()
        LG = [np.int_(np.eye(self.getDimension()))]
        E = self.getBase()
        for i in range(1, len(PG)):
            R = PG[i]
            Et = np.dot(R,E)
            A = np.dot(la.inv(E), Et)
            LG = np.append(LG,[np.int_(np.rint(A))], axis=0)
        return np.array(LG)

def inPointGroup(Q, L):
    '''
    check if Q is in the point group of L
    '''
    # if Q is not an orthogonal matrix, return false
    if not (abs(la.det(Q)) == 1.0 and (np.dot(Q, Q.T) == np.eye(len(Q))).all()):
        return False
    
    if isinstance(L, Lattice):
        L_old = L
        E = L.getBase()
        L_new = Lattice(np.dot(Q, E))
    elif isinstance(L, np.ndarray):
        try:
            L_old = Lattice(L)
        except:
            raise LatticeError('cannot construct lattice from L. L has to be a base matrix or an instance of Lattice.')
        try:
            L_new = Lattice(np.dot(Q, L))
        except:
            raise LatticeError('cannot construct lattice from the given Q and L. Q.L has to be a base matrix or L has to be an instance of Lattice.')
    else:
        raise LatticeError('cannot tell if Q is in the point group of L. Check the types of input parameters.')    
    return L_new==L_old
