import numpy as np
from numpy import dot
from numpy.linalg import inv, eig
from pystructrans.mat_math import mat_dot, mat_T, mat_eig

def Eric_dist(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim,dim)
    g1 = dot(E1.T, E1)
    g2 = dot(E2.T, E2)
    M = dot(L.T, dot(g1, L)) - g2
    # normalization factor
    nf = np.trace(dot(g2.T, g2))
    return np.trace(dot(M.T, M))/nf

def eric_dist(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim,dim)
    g1 = dot(E1.T, E1)
    g2 = dot(E2.T, E2)
    M = dot(L.T, dot(g1, L)) - g2
    # normalization factor
    nf = np.trace(dot(g2.T, g2))
    return np.trace(dot(M.T, M))/nf

def eric_dist_mat(l, E1, E2):
    N = len(l)  # number of l's
    dim = len(E1)
    # copy N times G and Gm
    g1 = np.vstack([E1.reshape(dim**2)]*N)
    g2 = np.vstack([E2.reshape(dim**2)]*N)
    
    E = mat_dot(g1, l)
    M = mat_dot(mat_T(E), E) - mat_dot(mat_T(g2), g2)
    # normalization factor
    g2 = dot(E2.T, E2)
    nf = np.trace(dot(g2.T, g2))
    return np.sum(M**2,1)/nf

def strain_dist(x, E1, E2inv):
    dim = len(E1)
    L = x.reshape(dim,dim)
    # Cinv = E1 L E2inv E2invT LT E1T
    Finv = dot(dot(E1, L), E2inv)
    Cinv = dot(Finv, Finv.T)
    lam = np.real(eig(Cinv)[0])
    slam = np.sqrt(lam)
    return np.sum(lam) - 2*np.sum(slam) + 3

def strain_dist_mat(l, E1, E2):
    N = len(l)  # number of l's
    dim = len(E1)
    # copy E1
    E1_list = np.vstack([E1.reshape(dim**2)]*N)
    # copy E2inv
    E2inv = inv(E2)
    E2i_list = np.vstack([E2inv.reshape(dim**2)]*N)
    # Finv
    Finv = mat_dot(E1_list, l)
    Finv = mat_dot(Finv, E2i_list) 
    # Cinv
    Cinv = mat_dot(Finv, mat_T(Finv))
    # eigenvalues of Cinv
    lams = mat_eig(Cinv)
    lams_sum = lams[:,0] + lams[:,1] + lams[:,2]
    # square roots of eigenvalues
    slams = np.sqrt(lams)
    slams_sum = slams[:,0]+slams[:,1]+slams[:,2]
    return lams_sum - 2.0*slams_sum + 3*np.ones_like(lams_sum)

def Cauchy_dist(x, E1, E2inv):
    dim = len(E1)
    L = x.reshape(dim,dim)
    # Cinv = E1 L E2inv E2invT LT E1T
    Finv = dot(dot(E1, L), E2inv)
    Cinv = dot(Finv, Finv.T)
    v = (Cinv - np.eye(dim)).reshape(dim**2)
    return np.dot(v.T, v)

def Cauchy_dist_mat(l, E1, E2):
    N = len(l)  # number of l's
    dim = len(E1)
    # copy E1
    E1_list = np.vstack([E1.reshape(dim**2)]*N)
    # copy E2inv
    E2inv = inv(E2)
    E2i_list = np.vstack([E2inv.reshape(dim**2)]*N)
    # Finv
    Finv = mat_dot(E1_list, l)
    Finv = mat_dot(Finv, E2i_list) 
    # Cinv = Finv*Finv.T
    Cinv = mat_dot(Finv, mat_T(Finv))
    # E = Cinv - I
    Eye = np.vstack([np.eye(dim).reshape(dim**2)]*N)
    E = Cinv - Eye
    return np.sum(E**2, 1)
    