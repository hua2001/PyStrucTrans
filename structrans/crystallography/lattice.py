from structrans.general_imports import *
from structrans.crystallography.matrix_group import MatrixGroup, CUBIC_LAUE_GROUP, HEX_LAUE_GROUP
from structrans.crystallography.matrix_group import SQUARE_GROUP, SQUARE_EXT_GROUP, HEX2D_EXT_GROUP, HEX2D_GROUP


class Lattice():
    """
    This is the base class of lattices.

    :attribute E: two-dimensional square matrix (any type can be converted to :py:class:`numpy.ndarray`)
    """
    
    def __init__(self, E):
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
        self.__pointGroup = None
        self.__latticeGroup = None
        self.__specLatticeGroup = None

    def base(self):
        """
        :return: the base matrix
        :rtype: :py:class:`numpy.ndarray`
        """
        E = np.empty_like(self.__E)
        E[:] = self.__E
        return E

    def getbase(self):
        warn("getbase() is deprecated, please use base() instead.", DeprecationWarning)
        return self.base()

    def dimension(self):
        """
        Get the dimenstion of the lattice.

        :return:  (*integer*) - the dimenstion of the lattice
        """
        return self.__N

    def inpointgroup(self, Q):
        """
        Check if `Q` is in the point group of the lattice.
        This method works for lattices in any dimension.
        """
        # if Q is not an orthogonal matrix, return false
        try:
            Q = np.array(Q)
        except Exception:
            return False
        if not _in_O3(Q):
            return False
        if len(Q) != self.__N:
            return False
        return self == Lattice(np.dot(Q, self.__E))

    def inlatticegroup(self, M):
        """
        Check if `M` is in the lattice group of the lattice.
        This method works for lattices in any dimension.
        """
        # if Q is not an orthogonal matrix, return false
        try:
            M = np.array(M)
        except Exception:
            return False
        M_int = np.rint(M)
        if not np.allclose(M, M_int, atol=1E-8):
            return False
        # QE = EM
        E = self.base()
        Q = E.dot(M).dot(la.inv(E))
        return _in_O3(Q)

    def Laue_group(self):
        '''
        :return: get the Laue group in :math:`SO(3)` of a 3D lattice, in :math:`SO(2)` of a 2D lattice.
        :rtype: :py:class:`MatrixGroup <structrans.crystallography.MatrixGroup>`
        :raises AttributeError: Laue groups for other than 2D and 3D lattices are not implemented yet.
        '''
        if not self.dimension() in [2, 3]:
            raise AttributeError('Laue groups for other than 2D and 3D lattices are not implemented yet.')
        if self.__LaueGroup is None and self.dimension() == 3:
            self.__LaueGroup = MatrixGroup(np.array([Q for Q in CUBIC_LAUE_GROUP.matrices() if self.inpointgroup(Q)]))
            if self.__LaueGroup.order() == 4:
                hexgroup = MatrixGroup(np.array([Q for Q in HEX_LAUE_GROUP.matrices() if self.inpointgroup(Q)]))
                if hexgroup.order() > self.__LaueGroup.order():
                    self.__LaueGroup = hexgroup
        if self.__LaueGroup is None and self.dimension() == 2:
            self.__LaueGroup = MatrixGroup(np.array([Q for Q in SQUARE_GROUP.matrices() if self.inpointgroup(Q)]))
            if self.__LaueGroup.order() == 2:
                hexGroup = MatrixGroup(np.array([Q for Q in HEX2D_GROUP.matrices() if self.inpointgroup(Q)]))
                if hexGroup.order() > self.__LaueGroup.order():
                    self.__LaueGroup = hexGroup
        return self.__LaueGroup

    def getLauegroup(self):
        warn("getLauegroup() is deprecated, please use Laue_group() instead.", DeprecationWarning)
        return self.Laue_group()

    def point_group(self):
        """
        :return: get the point group in :math:`O(3)` of a 3D lattice, in :math:`O(2)` of a 2D lattice
        :rtype: :py:class:`MatrixGroup <structrans.crystallography.MatrixGroup>`
        :raises AttributeError: Point groups for other than 2D and 3D lattices are not implemented yet.
        """
        if not self.dimension() in [2, 3]:
            raise AttributeError('Point groups for other than 2D and 3D lattices are not implemented yet.')
        if self.__pointGroup is None and self.__N == 3:
            lg = np.array(self.Laue_group().matrices())
            self.__pointGroup = MatrixGroup(np.append(lg, -lg, axis=0))
        elif self.__pointGroup is None and self.__N == 2:
            self.__pointGroup = MatrixGroup(np.array([Q for Q in SQUARE_EXT_GROUP.matrices() if self.inpointgroup(Q)]))
            if self.__pointGroup.order() == 4:
                hexGroup = MatrixGroup(np.array([Q for Q in HEX2D_EXT_GROUP.matrices() if self.inpointgroup(Q)]))
                if hexGroup.order() > self.__pointGroup.order():
                    self.__pointGroup = hexGroup
        return self.__pointGroup

    def getpointgroup(self):
        warn("getpointgroup() is deprecated, please use point_group() instead.", DeprecationWarning)
        return self.point_group()

    def special_lattice_group(self):
        """
        :return: get the special lattice group in :math:`SL(3, \mathbb Z)` of a 3D lattice
        :rtype: :py:class:`MatrixGroup <structrans.crystallography.MatrixGroup>`
        :raises AttributeError: Special lattice groups for other than 2D and 3D lattices are not implemented yet.
        """
        if not self.dimension() in [2, 3]:
            raise AttributeError('Special lattice groups for other than 2D and 3D lattices are not implemented yet.')
        if self.__specLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.Laue_group().matrices()], dtype='int')
            return MatrixGroup(mats)

    def getspeciallatticegroup(self):
        warn("getspeciallatticegroup() is deprecated, please use special_lattice_group() instead.", DeprecationWarning)
        return self.special_lattice_group()

    def lattice_group(self):
        """
        :return: get the lattice group in :math:`GL(3, \mathbb Z)` of a 3D lattice
        :rtype: :py:class:`MatrixGroup <structrans.crystallography.MatrixGroup>`
        :raises AttributeError: Lattice groups for other than 2D and 3D lattices are not implemented yet.
        """
        if not self.dimension() in [2, 3]:
            raise AttributeError('Lattice groups for other than 2D and 3D lattices are not implemented yet.')
        if self.__specLatticeGroup is None:
            E = self.__E
            Einv = la.inv(E)
            mats = np.array([np.rint(Einv.dot(Q.dot(E))) for Q in self.point_group().matrices()], dtype='int')
            return MatrixGroup(mats)

    def getlatticegroup(self):
        warn("getlatticegroup() is deprecated, please use lattice_group() instead.", DeprecationWarning)
        return self.special_lattice_group()

    def __eq__(self, other):
        if not isinstance(other, Lattice):
            return False

        if not self.dimension() == other.dimension():
            return False

        E1 = self.base()
        E2 = other.base()
        # E1 A = E2
        A = np.dot(la.inv(E1), E2)
        A_int = np.rint(A)
        return np.allclose(A, A_int, atol=1E-8)

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
    Q = np.array(Q)
    if not _is_square_mat(Q):
        return False
    else:
        if Q.dtype == np.int:
            return abs(la.det(Q)) == 1 and np.array_equal(np.dot(Q, Q.T), np.eye(len(Q)))
        else:
            return np.isclose(abs(la.det(Q)), 1.0) and np.allclose(np.dot(Q, Q.T), np.eye(len(Q)))