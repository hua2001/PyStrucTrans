from pystructrans.general_imports import *
import itertools
from functools import reduce
from operator import mul
from .util import divisors


def hnf_from_diag(diag):
    '''
    Find all the hermite normal forms from the given diagonal elements
    
    :param diag: diagonal elments [d1, ..., dN]
    :type diag: list or 1D array
    :return: [n x N x N] matrix, where n is the number of HNFs with the given diagonal.
    '''
    N = len(diag)
    def set_ij(hnfs, modif):
        newH = []
        for H in hnfs:
            if H[modif[0], modif[1]] == 0:
                M = H.copy()
                M[modif[0], modif[1]] = modif[2]
                newH.append(M)
        hnfs.extend(newH)
        return hnfs
    hnfs = [np.diag(diag)]
    modifs = ((i, j, -d) for i in xrange(N) for j in xrange(i) for d in xrange(1, int(diag[i])))
    reduce(set_ij, modifs, hnfs)
    return np.array(hnfs)

def hnf_from_det(n, N=3):
    '''    
    Find all Hermite Normal Forms with determinant n and dimension N.    
    
    :param n: the integer to be factorized
    :type n: integer
    :return: [9 x N] matrix, where N is the number of HNFs with the given diagonal.
    '''
    if n < 0:
        raise ValueError("negative determinant of HNF")
    # get all the divisors of n
    diags = (d for d in itertools.product(divisors(n), repeat=N) if reduce(mul, d) == n)
    hnfs = None
    for diag in diags:
        if hnfs is None:
            hnfs = hnf_from_diag(diag)
        else:
            hnfs = np.append(hnfs, hnf_from_diag(diag), axis=0)
    return hnfs
