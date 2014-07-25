from ..general_imports import *


class HashArray:
    """
    A class of hashable arrays.
    Main extra features are hashability and comparability.
    """
    def __init__(self, array, tol=1.E-12):
        self.array = np.array(array)
        # self.__key = list(self.array.shape)
        # self.__key.extend(self.array.flatten().tolist())
        self.__key = tuple(self.array.tolist())
        self.__tol = tol

    def key(self):
        return self.__key

    def __lt__(self, other):
        if not isinstance(other, HashArray):
            raise TypeError("unorderable tpyes: {:s} < {:s}".format(type(self).__name__, type(other).__name__))
        else:
            return self.array.tolist() < other.array.tolist()

    def __gt__(self, other):
        if not isinstance(other, HashArray):
            raise TypeError("unorderable tpyes: {:s} > {:s}".format(type(self).__name__, type(other).__name__))
        else:
            return self.array.tolist() > other.array.tolist()

    def __eq__(self, other):
        if not isinstance(other, HashArray):
            raise TypeError("unorderable tpyes: {:s} == {:s}".format(type(self).__name__, type(other).__name__))
        else:
            return self.array.tolist() == other.array.tolist()

    def __hash__(self):
        return hash(self.__key)

    def __str__(self):
        return str(self.array)