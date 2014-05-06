import numpy as np
from numpy import dot
from numpy.linalg import inv, eig
from pystructrans.mat_math import mat_dot, mat_T, mat_eig

def eric_dist(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim,dim)
    g1 = dot(E1.T, E1)
    g2 = dot(E2.T, E2)
    M = dot(L.T, dot(g1, L)) - g2
    # normalization factor
    nf = np.trace(dot(g2.T, g2))
    return np.trace(dot(M.T, M))/nf

def eric_dist_jac(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim,dim)
    E = dot(E1, L)
    M = dot(E.T, E) - dot(E2.T, E2)
    res = 2.0*dot(dot(E1.T, E),(M+M.T))
    # normalization factor
    g2 = dot(E2.T, E2)
    nf = np.trace(dot(g2.T, g2))
    return res.reshape(dim**2)/nf

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

def eric_dist_isnew(ln, l, LG1, LG2):
    # ln is new in l if none of the element in l
    # equals Q1*ln*Q2, Q1 in LG1, Q2 in LG2
    
    dim = len(LG1[0])
    ln = ln.reshape(dim**2)
    
    # pre and post multiply the group elements
    lg1 = LG1.reshape(len(LG1), dim**2)
    lg2 = LG2.reshape(len(LG2), dim**2)
    l1 = np.vstack([l]*len(lg1))
    g1 = np.tile(lg1,(1,len(l))).reshape(len(lg1)*len(l),dim**2)
    l2 = np.vstack([mat_dot(g1, l1)]*len(lg2))
    g2 = np.tile(lg2,(1,len(l)*len(lg1))).reshape(len(lg1)*len(l)*len(lg2),dim**2)
    c = mat_dot(l2, g2)
    # find the idx of the elements that are the same as the new L
    ex_idx = np.where(np.max(np.abs(c-ln), axis=1) < 1e-10)[0]
    if len(ex_idx) == 0:
        return True,np.array([], 'int16')
    else:
        # mod by the length of l to know which one in l is the same as ln
        ex_idx = np.array(list(set(np.mod(ex_idx, len(l)))), 'int16')
        return False,ex_idx
    
def eric_dist_unique(l, LG1, LG2):
    # unique upto Q1LQ2, Q1 in LG1, Q2 in LG2
    
    if len(l) == 1: # no duplication if l has only one element
        return [0]
    
    l_loc = l.copy()
    idx = np.arange(len(l))
    rem_idx = idx
    uni_idx = []
    
    while len(idx) > 0:
        uni_idx = np.append(uni_idx, [idx[0]])
        if len(idx) == 1:
            idx = []
        else:
            dup = np.array(eric_dist_isnew(l_loc[0], l_loc[1:], LG1, LG2)[1])+1 # duplicate to the first element
            rem_idx = [i for i in xrange(1,len(idx)) if i not in dup] # delete the first and all duplicated elements
            # update l and idx
            l_loc = l_loc[rem_idx]
            idx = idx[rem_idx]
    return np.array(uni_idx, 'i')

def dist_isnew(ln, l, LG1, LG2):
    # ln is new in l if none of the element in l
    # equals Q1*ln*Q2, Q1 in LG1, Q2 in LG2
    
    dim = len(LG1[0])
    ln = ln.reshape(dim**2)
    
    # pre and post multiply the group elements
    lg1 = LG1.reshape(len(LG1), dim**2)
    lg2 = LG2.reshape(len(LG2), dim**2)
    l1 = np.vstack([l]*len(lg1))
    g1 = np.tile(lg1,(1,len(l))).reshape(len(lg1)*len(l),dim**2)
    l2 = np.vstack([mat_dot(g1, l1)]*len(lg2))
    g2 = np.tile(lg2,(1,len(l)*len(lg1))).reshape(len(lg1)*len(l)*len(lg2),dim**2)
    c = mat_dot(l2, g2)
    # find the idx of the elements that are the same as the new L
    ex_idx = np.where(np.max(np.abs(c-ln), axis=1) < 1e-10)[0]
    
    if len(ex_idx) == 0:
        return True,np.array([], 'int16')
    else:
        # mod by the length of l to know which one in l is the same as ln
        ex_idx = np.array(list(set(np.mod(ex_idx, len(l)))), 'int16')
        return False,ex_idx

def dist_unique(l, LG1, LG2):
    # unique upto Q1LQ2, Q1 in LG1, Q2 in LG2
    
    if len(l) == 1: # no duplication if l has only one element
        return [0]
    
    l_loc = l.copy()
    idx = np.arange(len(l))
    rem_idx = idx
    uni_idx = []
    
    while len(idx) > 0:
        uni_idx = np.append(uni_idx, [idx[0]])
        if len(idx) == 1:
            idx = []
        else:
            dup = np.array(dist_isnew(l_loc[0], l_loc[1:], LG1, LG2)[1])+1 # duplicate to the first element
            rem_idx = [i for i in xrange(1,len(idx)) if i not in dup] # delete the first and all duplicated elements
#             logging.info('    found %d duplicated matrices. remain %d.' % (len(dup), len(rem_idx)))
            # update l and idx
            l_loc = l_loc[rem_idx]
            idx = idx[rem_idx]
    return np.array(uni_idx, 'i')

def strain_dist(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim,dim)
    # Cinv = E1 L E2inv E2invT LT E1T
    Finv = dot(dot(E1, L), inv(E2))
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

def Cauchy_dist(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim,dim)
    # Cinv = E1 L E2inv E2invT LT E1T
    Finv = dot(dot(E1, L), inv(E2))
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
    