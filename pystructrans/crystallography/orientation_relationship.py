import numpy as np
from numpy.linalg import inv

def direc_trans(L_cor, n_list):
    # convert input list to ndarray
    n_list = np.array(n_list)
    # assume not a 1d-array
    oned = False
    
    if len(n_list.shape) < 2:   # if it is a 1d-array
        # convert to 2d-array
        n_list = n_list.reshape(1,len(n_list))
        # makr the flag. 
        # we want to return an array having the same
        # dimension as the input 
        oned = True
    
    # L_cor * new_n = n ==> new_n  = L_cor^-1 * n    
    new_n_list = [np.dot(inv(L_cor),n_list[i]) for i in xrange(len(n_list))]
    
    if oned:
        return new_n_list[0]
    else:
        return np.array(new_n_list)

def plane_trans(L_cor, n_list):
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


