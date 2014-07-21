from ..general_imports import *
import logging
from .util import divisors

CUBIC_LAUE_GROUP = np.array([
    # identity
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    # two fold rotations
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # [100]
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # [010]
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],  # [001]
    [[0, 1, 0], [1, 0, 0], [0, 0, -1]],  # [110]
    [[0, -1, 0], [-1, 0, 0], [0, 0, -1]],  # [1-10]
    [[0, 0, 1], [0, -1, 0], [1, 0, 0]],  # [101]
    [[0, 0, -1], [0, -1, 0], [-1, 0, 0]],  # [10-1]
    [[-1, 0, 0], [0, 0, 1], [0, 1, 0]],  # [011]
    [[-1, 0, 0], [0, 0, -1], [0, -1, 0]],  # [01-1]
    # four fold rotations about [100]
    [[1, 0, 0], [0, 0, -1], [0, 1, 0]],
    [[1, 0, 0], [0, 0, 1], [0, -1, 0]],
    # four fold rotations about [010]
    [[0, 0, 1], [0, 1, 0], [-1, 0, 0]],
    [[0, 0, -1], [0, 1, 0], [1, 0, 0]],
    # four fold rotations about [001]
    [[0, -1, 0], [1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    # three-fold rotations about [111]
    [[0, 1, 0], [0, 0, 1], [1, 0, 0]],
    [[0, 0, 1], [1, 0, 0], [0, 1, 0]],
    # three-fold rotations about [-111]
    [[0, 0, -1], [-1, 0, 0], [0, 1, 0]],
    [[0, -1, 0], [0, 0, 1], [-1, 0, 0]],
    # three-fold rotations about [-1-11]
    [[0, 1, 0], [0, 0, -1], [-1, 0, 0]],
    [[0, 0, -1], [1, 0, 0], [0, -1, 0]],
    # three-fold rotations about [1-11]
    [[0, 0, 1], [-1, 0, 0], [0, -1, 0]],
    [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
])

# this is not really a group,
# because it contains no identity.
# Also it does not contain the 2 fold rotation about z-axis.
HEX_LAUE_GROUP = np.array([
    [[np.cos(np.pi / 3), np.sin(np.pi / 3), 0], [-np.sin(np.pi / 3), np.cos(np.pi / 3), 0], [0, 0, 1.]],
    [[np.cos(-np.pi / 3), np.sin(-np.pi / 3), 0], [-np.sin(-np.pi / 3), np.cos(-np.pi / 3), 0], [0, 0, 1.]],
    [[np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3), 0], [-np.sin(2 * np.pi / 3), np.cos(2 * np.pi / 3), 0], [0, 0, 1.]],
    [[np.cos(-2 * np.pi / 3), np.sin(-2 * np.pi / 3), 0], [-np.sin(-2 * np.pi / 3), np.cos(-2 * np.pi / 3), 0],
     [0, 0, 1.]],
    [[np.cos(np.pi / 3), np.sin(np.pi / 3), 0], [np.sin(np.pi / 3), -np.cos(np.pi / 3), 0], [0, 0, -1.]],
    [[np.cos(-np.pi / 3), np.sin(-np.pi / 3), 0], [np.sin(-np.pi / 3), -np.cos(-np.pi / 3), 0], [0, 0, -1.]],
    [[np.cos(2 * np.pi / 3), np.sin(2 * np.pi / 3), 0], [np.sin(2 * np.pi / 3), -np.cos(2 * np.pi / 3), 0],
     [0, 0, -1.]],
    [[np.cos(-2 * np.pi / 3), np.sin(-2 * np.pi / 3), 0], [np.sin(-2 * np.pi / 3), -np.cos(-2 * np.pi / 3), 0],
     [0, 0, -1.]]
])
FULL_LAUE_GROUP = np.append(CUBIC_LAUE_GROUP, HEX_LAUE_GROUP, axis=0)

# square group
SQUARE_GROUP = np.array([
    # identity
    [[1, 0], [0, 1]],
    # 90 rotations
    [[0, -1], [1, 0]],
    [[0, 1], [-1, 0]],
    # 180 rotation,
    [[-1, 0], [0, -1]]
])
SQUARE_GROUP_EXT = np.array([
    # mirror along x
    [[1, 0], [0, -1]],
    # mirror along y
    [[-1, 0], [0, 1]],
    # mirror along (1,1)
    [[0, 1], [1, 0]],
    # mirror along (-1,1)
    [[0, -1], [-1, 0]]
])
SQUARE_GROUP_EXT = np.append(SQUARE_GROUP, SQUARE_GROUP_EXT, axis=0)

class MatrixGroup():
    """
    Matrix representation of a group, in particular point groups and Laue groups.

    ``matrices`` should be a datatype that can be converted to a 3D :py:class:`numpy.ndarray`
    and the last two dimensions define a square matrix.
    """

    def __init__(self, matrices):
        self.__mtable = self.isgroup(matrices)
        if self.__mtable is False:
            raise ValueError("input does not form a group")
        self.__mats = matrices

    def multable(self):
        """
        :return: multiplication table
        :rtype: :py:class:`numpy.ndarray`
        """
        return self.__mtable

    def matrices(self):
        """
        :return: all the elements
        :rtype: :py:class:`numpy.ndarray`
        """
        return self.__mats

    def order(self):
        """
        :return: the group order
        :rtype: integer
        """
        return len(self.__mats)

    def isabelian(self):
        """
        :return: if the group is abelian
        :rtype: boolean
        """
        N = self.order()
        E = self.__mats[self.__mtable[0, 0]]   # id
        for i in xrange(N):
            if not np.array_equal(self.__mats[self.__mtable[i, i]], E):
                return False
            for j in xrange(i+1, N):
                if not np.array_equal(self.__mats[self.__mtable[i, j]], self.__mats[self.__mtable[j, i]]):
                    return False
        return True

    def hassubgroup(self, g):
        """
        :param g: another MatrixGroup
        :return: if ``g`` is a subgroup of the invoking group
        :rtype: boolean
        :raises ValueError: if `g` is not an instance of :py:class:`pystructrans.MatrixGroup`

        """
        if not isinstance(g, MatrixGroup):
            raise ValueError("input must be a MatrixGroup")

        if g.order() not in divisors(self.order()):
            return False

        for m1 in g.matrices():
            found = False
            for m2 in self.matrices():
                if np.max(np.abs(m1 - m2)) < 1.0E-12:
                    found = True
                    break
            if not found:
                return False

        return True

    @classmethod
    def isgroup(cls, matrices):
        """
        check is a list of `square matrices` form a group

        :param matrices:  anything can be converted in to a [`M` x `N` x `N`] :py:class:`numpy.ndarray`
        :return: return the multiplication table if is a group, otherwise return False
        :rtype: boolean or 2D :py:class:`numpy.ndarray`
        """
        try:
            matrices = np.array(matrices)
        except Exception:
            logging.debug("not array")
            return False
        shape = matrices.shape
        if len(shape) != 3 or shape[1] != shape[2]:
            logging.debug("wrong shape")
            return False

        dim = shape[1]
        hashmap = {}
        # existance of identity
        id_exist = False
        for i, m in enumerate(matrices):
            m[m == -0] = 0
            m[m == -0.0] = 0.0
            hashmap[m.tostring()] = i
            if np.max(np.abs(m - np.eye(dim))) < 1.0e-10:
                id_exist = True
        if not id_exist:
            logging.debug("no id")
            for m in matrices:
                logging.debug(m)
            return False

        N = len(matrices)
        if len(hashmap) != N:
            # exist duplication
            logging.debug("duplicate")
            return False

        mtable = np.empty((N, N), dtype='int')
        for i in xrange(N):
            for j in xrange(N):
                m = matrices[i].dot(matrices[j])
                m[m == -0] = 0
                m[m == -0.0] = 0.0
                if m.tostring() in hashmap:
                    mtable[i, j] = hashmap[m.tostring()]
                else:
                    logging.debug("false multable")
                    logging.debug("mi = {:s}".format(str(matrices[i])))
                    logging.debug("mj = {:s}".format(str(matrices[j])))
                    logging.debug("m = {:s}".format(str(m)))
                    return False

        return mtable
