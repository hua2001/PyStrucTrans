from pystructrans.general_imports import *
import itertools
from functools import reduce
from operator import mul
from ..util import divisors
import math


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
    :return: [N x 3 x 3] :py:class:`numpy.ndarray`
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


def __swap_cols(M, m, n):
    """swap column m and n of matrix M"""
    T = np.eye(M.shape[1], dtype='int')
    T[m, n] = T[n, m] = 1
    T[m, m] = T[n, n] = 0
    return M.dot(T)

def __add_col(M, c1, c2, t):
    """
    add t * column c2 to column c1
    """
    T = np.eye(M.shape[1], dtype='int')
    T[c2, c1] = t
    return M.dot(T)

def __nonneg_row(Min, r):
    """make row r of matrix M non-negative"""
    T = np.eye(Min.shape[1], dtype='int')
    M = Min.copy()
    for j, e in enumerate(M[r]):
        if e < 0:
            T[j, j] = -1
            M = M.dot(T)
            T[j, j] = 1
    return M

def __gcd_cols(Min, c1, c2):
    M = Min.copy()
    if M[0, c1] < M[0, c2]:
        M = __swap_cols(M, c1, c2)

    while M[0, c2] != 0:
        t = int(math.floor(M[0, c1]/M[0, c2]))
        M = __add_col(M, c1, c2, -t)
        M = __swap_cols(M, c1, c2)
    return M

def hnfd(M, onlyH=False):
    """
    Hermite normal form decomposition

    :param M: the input full-rank integer matrix
    :keyword onlyH: return only the Hermite normal form
    :return H:  the Hermite normal form, :py:class:`numpy.ndarray`
    :return L:  the transformation matrix :py:class:`numpy.ndarray`
    """
    H = np.array(M)
    if H.shape[1] == 1:
        if H[0, 0] > 0:
            return H if onlyH else (H, np.eye(1, dtype='int'))
        else:
            return -H if onlyH else (-H, -np.eye(1, dtype='int'))
    # first row non-negative
    H = __nonneg_row(H, 0)
    # make all but one element of first row zero
    for i in xrange(H.shape[1] - 1):
        H = __gcd_cols(H, - i - 2, - i - 1)
    # process right-bottom submatrix
    Hsub = hnfd(H[1:, 1:], onlyH=True)
    H[1:, 1:] = Hsub[:]
    # make sure H[i, i] is the greatest on row i
    for i in xrange(1, H.shape[0]):
        if H[i, 0] > 0:
            H = __add_col(H, 0, i, -int(math.ceil(H[i, 0]/H[i, i])))
        if H[i, 0] <= - H[i, i]:
            H = __add_col(H, 0, i, int(math.floor(-H[i, 0]/H[i, i])))
    if onlyH:
        return H
    else:
        return H, np.rint(la.inv(H).dot(np.array(M))).astype('int')