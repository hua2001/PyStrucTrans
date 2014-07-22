from ..general_imports import *
from ..crystallography import CUBIC_LAUE_GROUP, MatrixGroup, BravaisLattice, Lattice


class Martensite():
    r"""
    The class of martensitic phase transformations.

    * the first argument (:py:class:`numpy.ndarray`) is the transformation stretch tensor `U`
    * the second argument (:py:class:`pystructrans.MatrixGroup`) is the Laue group of the parent phase (usually chosen to be the high symmetry one),
      the default value is the full `cubic Laue group`

    Both arguments are optional, `i.e.` one can instantiate an empty Martensite object.

    """
    
    def __init__(self, *args):
        self.__U = args[0] if len(args) > 0 else None
        self.__Laue = args[1] if len(args) > 1 else CUBIC_LAUE_GROUP
        self.__Ulist = None
        self.__Laueidx = None

    def setU(self, *args):
        r"""
        Set the transformation stretch tensor of the transformation,
        and return a `new` Martensite object.
        """
        try:
            U = None
            if isinstance(args[0], np.ndarray) or isinstance(args[0], list):
                U0 = np.array(args[0])
                if U0.shape == (3, 3):
                   U = U0
                else:
                    raise ValueError
            elif len(args) == 1:
                U = np.diag([args[0], args[0], args[0]])
            elif len(args) == 2:
                U = np.diag([args[0], args[0], args[1]])
            elif len(args) == 3:
                U = np.array([[args[0], args[1], 0], [args[1], args[0], 0], [0, 0, args[2]]])
            elif len(args) == 4:
                U = np.array([[args[0], args[1], 0], [args[1], args[2], 0], [0, 0, args[3]]])
        except:
            raise ValueError("unrecognizable input")

        if not np.array_equal(U, U.T):
            raise ValueError("not symmetric")

        if not la.det(U) > 0:
            raise ValueError("not positive definite")

        return Martensite(U, self.getLaue())

    def getU(self):
        """
        :return: a `copy` of the transformation stretch matrix
        :rtype: :py:class:`numpy.ndarray`
        """
        if self.__U is None:
            raise AttributeError("U has not been initialized")
        U = np.empty_like(self.__U)
        U[:] = self.__U
        return U

    def getLaue(self):
        """
        :return: the Laue group of the high symmetry phase
        :rtype: :py:class:`pystructrans.MatrixGroup`
        """
        return self.__Laue

    def setLaue(self, arg):
        """
        set the Laue group, return a `new` Martensite object
        """
        latParam = [2, 2, 2,
                    [2, 3],
                    [2, 120],
                    [1.41, 2], [1.41, 2],
                    [1.41, 1.42, 2], [1.41, 1.42, 2], [1.41, 1.42, 2], [1.41, 1.42, 2],
                    [1.41, 1.42, 2, 89], [1.41, 1.42, 2, 89],
                    [1.41, 1.42, 2, 89, 88, 92]]

        if isinstance(arg, MatrixGroup):
            return Martensite(self.__U, arg)
        elif isinstance(arg, int) and 0 < arg <= 14:
            lattice = BravaisLattice(arg, latParam[arg-1])
            return Martensite(self.__U, lattice.getLauegroup())
        else:
            raise ValueError("unrecognizable input")

    def getvariants(self):
        """

        :return: all the variants
        :rtype: [N x 3 x 3] :py:class:`numpy.ndarray`
        :raises AttributeError:

            if the transformation is not reversible

        """
        if self.__Ulist is None:
            self.calcvariants()
        return self.__Ulist

    def getvariant(self, n):
        """

        :param n: index of the variant, an integer between 1 and # of variants
        :return: the `n`-th variant
        :rtype: :py:class:`numpy.ndarray`
        """
        if self.__Ulist is None:
            self.calcvariants()
        return self.__Ulist[n - 1]

    def isreversible(self):
        """
        It is reversible if the symmetry group of ``U``
        is a **proper subgroup** of the object's ``Laue`` group

        :return: whether the transformation is reversible
        :rtype: boolean
        :raises AttributeError:

                if `U` has not been assigned
        """
        if self.__U is None:
            raise AttributeError("U has not been initialized")
        lg = Lattice(self.getU()).getLauegroup()
        lg0 = self.getLaue()
        return lg.order() < lg0.order() and lg0.hassubgroup(lg)

    def getLaueidx(self):
        if self.__Laueidx is None:
            self.calcvariants()
        return self.__Laueidx

    def calcvariants(self):
        """
        (lazy) calculate all the variants of the transformation.
        store the matrices associated to each variant
        and the indices of elements in the Laue group that maps
        the original `U` to them.
        """
        if self.__Ulist is None:
            if not self.isreversible():
                raise AttributeError("irreversible martensitic transformations have no variants")

            U1 = self.getU()
            ulist = [U1]
            idx = [[]]
            for i, Q in enumerate(self.__Laue.matrices()):
                V = Q.dot(U1).dot(Q.T)
                newU = True
                for j, U in enumerate(ulist):
                    if np.max(np.abs(U - V)) < SMALL:
                        newU = False
                        idx[j].append(i)

                if newU:
                    ulist.append(V)
                    idx.append([i])

            self.__Ulist = np.array(ulist)
            self.__Laueidx = np.array(idx)

#     def getLaueIdx(self, *args):
#         if self._laueIdx == None:
#             self.calcVariants()
#         if len(args) == 0:
#             return self._laueIdx
#         else:
#             return self._laueIdx[args[0]-1]
#
#     def getLaue_M(self):
#         idx = self._laueIdx[0]
#         return np.array([self._Laue[i] for i in idx])
#
#     def getCor_list(self, E, L):
#         r'''
#         E should be the primitive lattice base
#         '''
#         if isinstance(E, np.ndarray) and isinstance(E, np.ndarray):
#             if E.shape == (3,3) and E.shape == (3,3):
#                 laue_idx = self.getLaueIdx()
#                 lg_a = self.getLaue()
#                 lcor = []
#                 for i in xrange(len(laue_idx)):
#                     ltemp = []
#                     for j in xrange(len(laue_idx[0])):
#                         ltemp.append((np.dot(inv(E), lg_a[laue_idx[i][j]].dot(E))).dot(L))
#                     lcor.append(ltemp)
#                 return np.array(lcor)
#             else:
#                 raise MartensiteError('Please input lattice base vectors E in R^{3x3} and correspondence matrix L.')
#         else:
#             raise MartensiteError('Please input lattice base vectors E in R^{3x3} and correspondence matrix L.')
#     def getLc_list(self, C, Lc):
#         r'''
#         C should be the conventional lattice base
#         '''
#         if isinstance(C, np.ndarray) and isinstance(C, np.ndarray):
#             if C.shape == (3,3) and C.shape == (3,3):
#                 laue_idx = self.getLaueIdx()
#                 lg_a = self.getLaue()
#                 lcor = []
#                 for i in xrange(len(laue_idx)):
#                     ltemp = []
#                     for j in xrange(len(laue_idx[0])):
#                         ltemp.append((np.dot(inv(C), lg_a[laue_idx[i][j]].dot(C))).dot(Lc))
#                     lcor.append(ltemp)
#                 return np.array(lcor)
#             else:
#                 raise MartensiteError('Please input lattice base vectors E in R^{3x3} and correspondence matrix L.')
#         else:
#             raise MartensiteError('Please input lattice base vectors E in R^{3x3} and correspondence matrix L.')
        
        
        
        
        
