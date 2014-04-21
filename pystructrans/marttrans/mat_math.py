'''
This file contains functions and data related to 
matrix operation of a large number of matrices
'''
import sys
import os
import logging
import string
import random
import h5py
import numpy as np
from math import sqrt, modf
try:
    from scipy.linalg import eig
except:
    from numpy.linalg import eig
import H5_FILE

def mat_dot (la, lb):
    '''
    multiplying list of dxd matrices in the form of 1x9 vectors
    the dimension of la and lb are both [N x d**2]
    '''
    d = int(sqrt(len(la[0]))) # dimension
    lc = np.zeros_like(la)
    
    
    for i in xrange(d):
        for j in xrange(d):
#             lc[:, i*d+j] = np.zeros(len(la))
            for k in xrange(d):
                # C_ij = A_ik*B_kj                
                lc[:, i*d+j] += la[:,i*d+k] * lb[:,k*d+j]
                    
    return lc
    
    
def mat_trace (l):
    # dimension
    d = int(sqrt(len(l[0])))
    
    tr = np.empty(len(l), dtype=l.dtype)    
        
    for i in xrange(d):
        tr += l[:, i*d+i]
    return tr

def mat_det (l):    
    try:
        d = int(sqrt(len(l[0])))  
    except:
        d = 1
    
    if d == 1:
        return l.flatten()
#         try:
#             det = l
#         except:
#             det = l.T
#         return det
    
    det = np.zeros(len(l), dtype=l.dtype)        
    for n in xrange(d):
        # calculate the minors
        s = (-1)**n # sign
        # remove the 0-th row and the n-th column
        rm_idx = range(d)
        rm_idx.extend([i*d+n for i in xrange(1,d)])
        #rm_idx = set(rm_idx)
        # idx of the minors
        mn_idx = [i for i in xrange(d**2) if not i in rm_idx]
        # the n-th minor    
        minor = mat_det(l[:, mn_idx])
        det[:] += s*l[:, n]*minor  
           
    return det

def mat_T(l):
    d = int(sqrt(len(l[0])))
    
    T = np.empty_like(l) 
        
    for i in xrange(d):
        for j in xrange(d):
            T[:, i*d+j] = l[:, j*d+i]
            
    return T

def mat_eig(l, MAX_MEM=1E9):
    dim = sqrt(l.shape[1])
    
    eigs = np.empty((len(l),dim))
    
    for i in xrange(len(l)):
        M = l[i].reshape(dim, dim)
        eigs[i,:] = np.real(eig(M)[0])
    return eigs
    

def unique_rows(a):
    
    if a.dtype == 'int16':
        N_MAX = H5_FILE.N_MAX/4
    else:
        N_MAX = H5_FILE.N_MAX/16
        
    if len(a) <= N_MAX:
        if isinstance(a, h5py._hl.dataset.Dataset):
            a = a[()]
        b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
        _, idx = np.unique(b, return_index=True)
        return a[idx], idx
    else:
        new_len = 0
        old_len = len(a)
        while new_len < old_len:
            logging.warning('%d matrices are too many, using HDF5 file' % len(a))                
            n_total = len(a)
            n_done = 0
            idx = np.array([])
            while n_done < n_total:
                n_slice = min(N_MAX, n_total-n_done)
                logging.info('checking uniqueness of %s to %s elements out of %d' 
                             % (H5_FILE.num2ord(n_done+1), H5_FILE.num2ord(n_done+n_slice), len(a)))
                a_slice = np.array(a[n_done:n_done+n_slice,:])
                b = np.ascontiguousarray(a_slice).view(np.dtype((np.void, a_slice.dtype.itemsize * a_slice.shape[1])))
                _, c = np.unique(b, return_index=True)
                c = np.sort(c)
                begin = len(idx)
                if begin == 0:
                    idx = H5_FILE.H5_FILE.create_dataset(H5_FILE.rand_str(6), (len(c),), maxshape=(len(a),), dtype='int')
                    d = H5_FILE.H5_FILE.create_dataset(H5_FILE.rand_str(6),(len(c), a.shape[1]),maxshape=a.shape,dtype=a.dtype)
                else:
                    d.resize((len(idx)+len(c), a.shape[1]))
                    idx.resize((len(idx)+len(c),))
                idx[begin:] = c + n_done
                d[begin:] = a_slice[c]
                
                del c
                del a_slice
                del b
                
                n_done = n_done+n_slice
                       
            # update new and old length
            old_len = len(a)
            a = d 
            new_len = len(a)   
            if new_len <= N_MAX:
                return unique_rows(a[()])                 
        return a, idx
    
    
