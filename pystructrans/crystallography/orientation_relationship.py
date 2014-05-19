import numpy as np
from numpy.linalg import inv, det, eig, norm
from math import sqrt


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

def AM_Solver(A):
    
    kappa = 1
    
    if len(A) - 3 > 1e-9:
        print 'The input should be a 3x3 symmetric matrix!'
        return
      
    if norm((A.T - A).reshape(9)) > 1e-9:
        print 'The input should be a symmetric matrix!'
        return 
    
    e = eig(A)
    e = [np.append(e[0][i],e[1][:,i]).real for i in xrange(3)]
    e = np.array(sorted(e, key=lambda x:x[0]))
    eval = e[:,0]
    evec = np.array([[1.,0,0],[0,1.,0],[0,0,1.]]).dot(e[:,1:4])
    
    
    c = sqrt(eval[2]-eval[0])
    c1 = sqrt(1-eval[0])
    c3 = sqrt(eval[2]-1)
    
    if c < 1e-6:
        print 'solution is b = e, m = - 2 e where |e| = 1.'
        return
    else:
        if abs(eval[1] - 1) <1e-4:
            m1 = ((sqrt(eval[2])-sqrt(eval[0]))/c)*(-c1*evec[0] + kappa*c3*evec[2])
            m2 = ((sqrt(eval[2])-sqrt(eval[0]))/c)*(-c1*evec[0] - kappa*c3*evec[2])
            rho1 = norm(m1)
            rho2 = norm(m2)
            m1 = m1/rho1
            m2 = m2/rho2
            b1 = rho1*((sqrt(eval[2])*c1/c)*evec[0] + kappa*(sqrt(eval[0])*c3/c)*evec[2])
            b2 = rho2*((sqrt(eval[2])*c1/c)*evec[0] - kappa*(sqrt(eval[0])*c3/c)*evec[2])
        
            return np.array([[b1,m1], [b2,m2]])

def TwinSolver(U, e, type):
    
    e = np.array(e)
    
    if len(U) - 3 > 1e-9:
        print 'The input should be a 3x3 symmetric matrix!'
        return
      
    if norm((U.T - U).reshape(9)) > 1e-9:
        print 'The input should be a symmetric matrix!'
        return
    if abs(type - 1) < 1e-9:
        n = e
        denominator = np.dot(inv(U).dot(e), inv(U).dot(e))
#         print denominator
        a = 2*(np.dot(inv(U),e)/denominator - U.dot(e))
    else:
        if abs(type - 2) < 1e-9:
            n = 2*(e - np.dot(U, U.dot(e))/np.dot(U.dot(e), U.dot(e)))
            a = U.dot(e)
        else:
            print 'Please input the type of twin system: 1 for Type I twin; 2 for Type II twin'
            return
    
    return np.array([a, n])
            
        
    
    
        
    
