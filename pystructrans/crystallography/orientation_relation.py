import numpy as np
import numpy.linalg as la
from numpy.linalg import inv
import itertools

def tuple(start, end):
    x = np.arange(start, end)
    a = [x,x,x]
    Z = np.array(list(itertools.product(*a)))
    return Z

def vec_trans(L_cor, n_list):
    n_list = np.array(n_list)
    oned = False
    if len(n_list.shape) < 2:
        n_list = n_list.reshape(1,len(n_list))
        oned = True
    new_n_list = [np.dot(inv(L_cor),n_list[i]) for i in xrange(len(n_list))]
    if oned:
        return new_n_list[0]
    else:
        return np.array(new_n_list)

def Rvec_trans(L_cor, n_list):
    n_list = np.array(n_list)
    oned = False
    if len(n_list.shape) < 2:
        n_list = n_list.reshape(1,len(n_list))
        oned = True
    new_n_list = [np.dot(L_cor.T, n_list[i]) for i in xrange(len(n_list))]
    if oned:
        return new_n_list[0]
    else:
        return np.array(new_n_list)


