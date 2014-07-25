from pystructrans.general_imports import *
from .matrix_group import MatrixGroup, CUBIC_LAUE_GROUP, HEX_LAUE_GROUP, SQUARE_GROUP, SQUARE_GROUP_EXT

class Lattice():
    '''
    This is the base class of lattices.

    :param E: two-dimentional square matrix, any type can be converted to :py:class:`numpy.ndarray`
    '''   
    
    def __init__(self, E):
        '''
        Constructor
        '''
        # check if E is a square matrix
        try:
            if not isinstance(E, np.ndarray):
                E = np.array(E)
            if not _is_square_mat(E):
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

    def getbase(self):
        '''

        :return: the base matrix
        :rtype: :py:class:`numpy.ndarray`
        '''
        E = np.empty_like(self.__E)
        E[:] = self.__E
        return E

    def getdimension(self):
        '''
        Get the dimenstion of the lattice.

        :return:  (*integer*) - the dimenstion of the lattice
        '''
        return self.__N

    def inpointgroup(self, Q):
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

    def getLauegroup(self):
        '''

        :return: get the Laue group in SO(3) of a 3D  lattice
        :rtype: :py:class:`pystructrans.crystallography.MatrixGroup`
        '''
        if not self.getdimension() in [2, 3]:
            raise AttributeError('Only Laue groups for lattices in 2D and 3D have been implemented.')
        if self.__LaueGroup is None and self.getdimension() == 3:
            self.__LaueGroup = MatrixGroup(np.array([Q for Q in CUBIC_LAUE_GROUP.matrices() if self.inpointgroup(Q)]))
            hexgroup = MatrixGroup(np.array([Q for Q in HEX_LAUE_GROUP.matrices() if self.inpointgroup(Q)]))
            if hexgroup.order() > self.__LaueGroup.order():
                self.__LaueGroup = hexgroup
        if self.__LaueGroup is None and self.getdimension() == 2:
            self.__LaueGroup = MatrixGroup(np.array([Q for Q in SQUARE_GROUP if self.inpointgroup(Q)]))
        return self.__LaueGroup

    def getpointgroup(self):
        '''

        :return: get the point group in O(3) of a 3D lattice
        :rtype: :py:class:`pystructrans.crystallography.MatrixGroup`
        '''
        if not self.getdimension() in [2, 3]:
            raise AttributeError('Only point groups for lattices in 2D and 3D have been implemented.')
        if self.__PointGroup is None and self.__N == 3:
            lg = self.getLauegroup().matrices()
            self.__PointGroup = MatrixGroup(np.append(lg, -lg, axis=0))
        elif self.__PointGroup is None and self.__N == 2:
            self.__PointGroup = MatrixGroup(np.array([Q for Q in SQUARE_GROUP_EXT if self.inpointgroup(Q)]))
        return self.__PointGroup

    def getspeciallatticegroup(self):
        '''

        :return: get the special lattice group in SL(3, Z) of a 3D lattice
        :rtype: :py:class:`pystructrans.crystallography.MatrixGroup`
        '''
        if not self.getdimension() in [2, 3]:
            raise AttributeError('Only special lattice groups for lattices in 2D and 3D have been implemented.')
        if self.__SLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.getLauegroup().matrices()], dtype='int')
            return MatrixGroup(mats)

    def getlatticegroup(self):
        '''

        :return: get the lattice group in GL(3, Z) of a 3D lattice
        :rtype: :py:class:`pystructrans.crystallography.MatrixGroup`
        '''
        if not self.getdimension() in [2, 3]:
            raise AttributeError('Only lattice groups for lattices in 2D and 3D have been implemented.')
        if self.__SLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.getpointgroup().matrices()], dtype='int')
            return MatrixGroup(mats)

    def __eq__(self, other):
        if not isinstance(other, Lattice):
            return False

        E1 = self.getbase()
        E2 = other.getbase()
        if not self.getdimension() == other.getdimension():
            return False

        A = np.dot(la.inv(E1), E2)
        A_int = np.round(A)
        return nanmax(np.abs(A-A_int)) < 1.0E-10

    def __str__(self):
        des_str = '{:d} dimensional lattice:'.format(self.__N)
        for i in range(self.__N):
            des_str = des_str + '\n    e_{:d} = {:s}'.format(i+1, str(self.__E[:, i]))
        return des_str

def _is_square_mat(M):
    if not isinstance(M, np.ndarray):
        M = np.array(M)
    return (M.ndim == 2) and (M.shape[0] == M.shape[1])

def _in_O3(Q):
    if not _is_square_mat(Q):
        return False
    else:
        return abs(la.det(Q)) == 1.0 and (np.dot(Q, Q.T) == np.eye(len(Q))).all()