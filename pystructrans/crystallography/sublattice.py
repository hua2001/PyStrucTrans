import numpy as np

def _divisors(n):
    '''
    Find all divisors of an integer n
    
    :param n: the integer to be factorized
    :type n: integer
    :return: vector - the list of all the divisors of n
    '''

    d = np.array([])
    for i in range(np.int_(np.floor(n/2))):
        if np.remainder(n,i+1) == 0:
            d = np.append(d, [i+1], axis=0)
    if n==1:
        d = np.array([1])
    else:
        d = np.append(d, [n], axis=0)
            
    return np.int_(d)

def hnf_from_diag(diag):
    '''
    Find all the hermite normal forms from the given diagonal elements
    
    :param diag: diagonal elments [d1, ..., dN]
    :type diag: list or 1D array
    :return: [n x N^2] matrix, where n is the number of HNFs with the given diagonal.
    '''        
    N = len(diag)
    H = np.zeros((1, N**2))
    for i in range(N):
        H[0, i+N*i] = diag[i]
    
    for i in range(1, N):
        for j in range(i):
            # modify the (i,j) component
            # copy the H list as a temporary new list
            H_new = np.empty_like(H)
            H_new[:] = H
            for d in range(1, int(max(diag))):
                if d < diag[i]:
                    # modify the (i,j) component of all the new matrices
                    H_new[:, j+N*i] = -d
                    # stack the new list to the old one
                    H = np.vstack((H, H_new))
    return H

def hnf_from_det(n, N=3):
    '''    
    Find all Hermite Normal Forms with determinant n and dimension N.    
    
    :param n: the integer to be factorized
    :type n: integer
    :return: [9 x N] matrix, where N is the number of HNFs with the given diagonal.
    '''
    # get all the divisors of n
    ds = _divisors(n)
    # list all the combinations of the divisors
    diags = np.ones((1,N))  
    if n != 1: # if more than one divisors
        for i in range(N):
            # copy the list of diagonals
            diags_new = diags.copy()
            for d in ds[1:]:
                diags_new[:,i] = d
                diags = np.vstack((diags, diags_new))
    
    # find the diagonals giving the correct determinant
    PI = np.ones(len(diags))
    for i in range(N):
        PI = PI * diags[:,i]
    diags = diags[np.where(PI==n)]
    
    # generate Hermite normal forms from the list of diagonals
    H = np.zeros([1, N**2])
    for row in diags:
        H = np.vstack((H, hnf_from_diag(row)))
    del diags
    return H[1:]
    