'''
This file contains functions and data related to 
matrix operation of a large number of matrices
'''
import numpy as np
from math import sqrt
try:
    from scipy.linalg import eig
except:
    from numpy.linalg import eig
import itertools

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
    
    det = np.zeros(len(l), dtype=l.dtype)        
    for n in xrange(d):
        # calculate the minors
        s = (-1)**n # sign
        # remove the 0-th row and the n-th column
        rm_idx = range(d)
        rm_idx.extend([i*d+n for i in xrange(1,d)])
        # idx of the minors
        mn_idx = [i for i in xrange(d**2) if not i in rm_idx]
        # the n-th minor    
        minor = mat_det(l[:, mn_idx])
        det[:] += s*l[:, n]*minor  
           
    return det

def mat_T(l):
    # dimension of the matrices
    dim = int(sqrt(len(l[0])))   
    # initiate an array 
    T = np.empty_like(l)         
    for i in xrange(dim):
        for j in xrange(dim):
            T[:, i*dim+j] = l[:, j*dim+i]            
    return T

def mat_eig(l, MAX_MEM=1E9):
    dim = sqrt(l.shape[1])    
    eigs = np.empty((len(l),dim))    
    for i in xrange(len(l)):
        M = l[i].reshape(dim, dim)
        eigs[i,:] = np.real(eig(M)[0])
    return eigs    

def unique_rows(a):
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    return a[idx], idx
    
def tuple(start, end):
    x = np.arange(start, end)
    a = [x,x,x]
    Z = np.array(list(itertools.product(*a)))
    return Z 
