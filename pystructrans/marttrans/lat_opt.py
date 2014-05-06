import numpy as np
import h5py
import sys
import logging
import os.path

from numpy import dot
from numpy.linalg import inv
from numpy import nanmin

from timeit import default_timer as timer
from dist import eric_dist, eric_dist_mat, eric_dist_unique
from dist import strain_dist, strain_dist_mat
from dist import Cauchy_dist, Cauchy_dist_mat
from pystructrans.mat_math import mat_det, mat_dot
from pystructrans.crystallography import Lattice, LLL

def lprint(s, rank, disp):
    msg = ('Processor %d: ' % rank)+s
    if disp:
        logging.info(msg)

class NullDevice():
    def write(self, s):
        pass

def lat_opt(E1, E2, distance='Ericksen', dim=3, num_sol=1, disp=True, rank=0, nhnf = 1):
    '''
    find the lattice invariant shear L_opt that brings the unit cell E1 closest to E2.
    '''
    if not disp:
        original_stdout = sys.stdout
        sys.stdout = NullDevice()
    log_writer = lambda s: lprint(s, rank, disp)
    
    # set distance functions
    if distance == 'Ericksen':
        dist=eric_dist
#         dist_jac=eric_dist_jac
        dist_mat=eric_dist_mat
        dist_unique = eric_dist_unique
    if distance == 'strain':
        dist=strain_dist
        dist_mat=strain_dist_mat
        dist_unique = eric_dist_unique
    if distance == 'Cauchy':
        dist=Cauchy_dist
        dist_mat=Cauchy_dist_mat
        dist_unique = eric_dist_unique
        
    # start timer
    t_start = timer()
    '''
    find LLL-reduced basis for the sublattice
    '''
    #log_writer(str(E1))
    Ep2 = LLL(E2)
    Er = LLL(E1)
    lo = np.rint(dot(inv(E1), Er).reshape(dim**2))
    chi = np.rint(inv(Ep2).dot(E2).reshape(dim**2))
    
    # search shifts up to 1 lattice point
    #print('Shifting matrices are generated.')    
    h5_filename = 'shift_1_dim_'+str(dim)+'.hdf5'
    data = os.path.join(os.path.dirname(__file__),h5_filename)
    f = h5py.File(data, 'r')
    A = f['shifts'].value
    # add all zero
    A = np.append(A, np.zeros((dim**2,1)), axis=1).T
    f.close()    
    
    # copy lo 3^9 times and add to the list of shift matrices 
    l_shift = A + np.vstack([lo]*len(A))
    # get all the L's with determinant 1    
    l_det = mat_det(l_shift)
    det_idx = np.where(np.abs(l_det-1.0)<1.e-6)[0]
    l_list = l_shift[det_idx]       # get all the L's with determinant 1 
    # release memory
    del A
    del f
    
    # timing integer points search
    t_int = timer()
    
    # calculate all the distance
    if len(l_list)==0:
        log_writer('no determinant 1')
        return np.empty([]), None, False
    d_list = dist_mat(l_list, E1, Ep2)
     
    '''
    look for the required number of solutions
    '''
    # initialization
    dopt = np.array([])
    lopt = np.zeros((1,dim**2))
    # calculate the distance for all the L's
#     d_list = dist_mat(l_list, E1, E2)
    # Instanciate a Lattice object for E1
    Lat1 = Lattice(E1)
    # lattice group of E1
    LG1 = Lat1.getSpecialLatticeGroup().astype('float')
    # Instanciate a Lattice object for E2
    Lat2 = Lattice(Ep2)
    # lattice group of E2
    LG2 = Lat2.getSpecialLatticeGroup().astype('float') 
    SOLG = LG2
#     SOLG = np.array([np.eye(dim).reshape(dim**2)]) 
#     for i,Q in enumerate(LG2):
#         if inPointGroup(Q, np.eye(dim)):
#             # if Q is not in SO(d), add to SOLG
#             SOLG = np.append(SOLG, Q.reshape((1,dim**2)), axis=0)   
    while len(dopt) < num_sol:
        # find minimum distances
        min_dist = nanmin(d_list)
        min_idx = np.where((d_list-min_dist)<=(1E-6)*min_dist)[0]
        # remove duplications caused by symmetry:
        # i.e. the distance between these two lattices is zero
        sol_idx = min_idx[dist_unique(l_list[min_idx], LG1, SOLG)]
        
        # store new solutions and corresponding distances in the return paramters
        lopt = np.append(lopt, l_list[sol_idx], axis=0)
        dopt = np.append(dopt, [min_dist]*len(sol_idx))
#         log_writer('the same distance %s but not symmetry related lopt'%(str(dopt)))
#         log_writer('%s'%(str(lopt)))
        
        
        
        # remove found solutions from l_list
        rem_idx = [i for i in xrange(len(l_list)) if i not in min_idx]
        if len(rem_idx) == 0:
            # no remaining L's
            log_writer('only found %d solutions for this HNF, less than the required number %d!' % (len(dopt), num_sol))
            return lopt[1:], dopt, True
        l_list = l_list[rem_idx]
        d_list = d_list[rem_idx]
        
    # if more than required solutions are found
    if len(dopt) > num_sol:
        log_writer('found %d solutions, which is more than %d!' % (len(dopt), num_sol))
        
    # stop timer
    t_end = timer()
    
    # resume stdout
    if not disp:
        sys.stdout = original_stdout
    
    ltemp = mat_dot(lopt[1:], np.vstack([chi]*(len(lopt)-1)))
    lopt = ltemp
    
    # return solutions
    return np.rint(lopt), dopt, True
