from structrans.general_imports import *
from itertools import chain
from structrans.util import divisors, rotation

class MatrixGroup():
    """
    Matrix representation of a group, in particular point groups and Laue groups.

    ``matrices`` should be a datatype that can be converted to a 3D :py:class:`numpy.ndarray`
    and the last two dimensions define a square matrix.
    """

    @classmethod
    def isgroup(cls, matrices):
        """
        check is a list of `square matrices` form a group

        :param matrices:  anything can be converted in to a [`M` x `N` x `N`] :py:class:`numpy.ndarray`
        :return: return the multiplication table if is a group, otherwise return False
        :rtype: boolean or 2D :py:class:`numpy.ndarray`
        """
        try:
            mats = np.round(np.array(matrices), 6)
        except Exception:
            return False
        shape = mats.shape
        if len(shape) != 3 or shape[1] != shape[2]:
            return False

        dim = shape[1]
        hashmap = {}
        # existance of identity
        id_exist = False
        for i, m in enumerate(mats):
            if tuple(m.flatten()) in hashmap:
                logging.debug("duplicate")
                logging.debug("i = {:d}, md = {:s}".format(i, str(m)))
                return False
            hashmap[tuple(m.flatten())] = i
            if np.allclose(m, np.eye(dim)):
                id_exist = True
        if not id_exist:
            logging.debug("no id")
            for m in matrices:
                logging.debug(m)
            return False

        N = len(mats)

        def ingroup(M):
            for i, m in enumerate(mats):
                if nanmax(np.abs(M - m)) < 1.0e-5:
                    return i
            return -1

        mtable = np.empty((N, N), dtype='int')
        for i in xrange(N):
            for j in xrange(N):
                m = np.round(mats[i].dot(mats[j]), 6)
                idx = ingroup(m)
                if idx > -1:
                    mtable[i, j] = idx
                else:
                    logging.debug("false multable")
                    logging.debug("mi = {:s}".format(str(mats[i])))
                    logging.debug("mj = {:s}".format(str(mats[j])))
                    logging.debug("m = {:s}".format(str(m)))
                    return False

        return mtable

    def __init__(self, matrices):
        self.__mtable = None
        # self.__mtable = self.isgroup(matrices)
        # if self.__mtable is False:
        #     raise ValueError("input does not form a group")
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
        return [m[:] for m in self.__mats]

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
        for i in xrange(N):
            for j in xrange(i + 1, N):
                if self.__mtable[i, j] != self.__mtable[j, i]:
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
            raise TypeError("input must be a MatrixGroup")

        if g.order() not in divisors(self.order()):
            return False

        for m1 in g.matrices():
            found = False
            for m2 in self.matrices():
                if np.allclose(m1, m2):
                    found = True
                    break
            if not found:
                return False

        return True


cubic_Laue_group_mats = np.array([
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
CUBIC_LAUE_GROUP = MatrixGroup(cubic_Laue_group_mats)

# this is not really a group,
# because it contains no identity.
# Also it does not contain the 2 fold rotation about z-axis.
__C1 = np.cos(np.pi / 3.0)
__S1 = np.sin(np.pi / 3.0)
hex_Laue_group_mats = np.array([
    # identity
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    # two fold rotations
    [[1, 0, 0], [0, -1, 0], [0, 0, -1]],  # [100]
    [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],  # [010]
    [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],  # [001]
    [[__C1, -__S1, 0], [__S1, __C1, 0],  [0, 0, 1]],
    [[__C1, __S1, 0], [-__S1, __C1, 0],  [0, 0, 1]],
    [[-__C1, -__S1, 0], [__S1, -__C1, 0],  [0, 0, 1]],
    [[-__C1, __S1, 0], [-__S1, -__C1, 0],  [0, 0, 1]],
    [[__C1, -__S1, 0], [-__S1, -__C1, 0],  [0, 0, -1]],
    [[__C1, __S1, 0], [__S1, -__C1, 0],  [0, 0, -1]],
    [[-__C1, -__S1, 0], [-__S1, __C1, 0],  [0, 0, -1]],
    [[-__C1, __S1, 0], [__S1, __C1, 0],  [0, 0, -1]]
])
HEX_LAUE_GROUP = MatrixGroup(hex_Laue_group_mats)
# HEX_LAUE_GROUP = CUBIC_LAUE_GROUP

# square group
SQUARE_GROUP = MatrixGroup(np.array([
    # identity
    [[1, 0], [0, 1]],
    # 90 rotations
    [[0, -1], [1, 0]],
    [[0, 1], [-1, 0]],
    # 180 rotation,
    [[-1, 0], [0, -1]]
]))
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
SQUARE_GROUP_EXT = MatrixGroup(np.append(SQUARE_GROUP.matrices(), SQUARE_GROUP_EXT, axis=0))