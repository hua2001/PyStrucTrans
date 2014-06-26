import numpy as np
import math
import os
import sys
try:
    import h5py
except:
    pass
import logging

from numpy import dot
from numpy.linalg import inv, det, eig
            
from timeit import default_timer as timer

from lat_opt import lat_opt
from dist import dist_isnew, dist_unique
from mat_math import mat_dot
from crystallography import BravaisLattice, HermiteNormalForms

def _divide_work(W, comm): # divide W into sub-processes
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if rank == 0:
        if size >= len(W): # processes more than HNFs
            w_idx = -np.ones((size,1))
            w_idx[:len(W),0] = range(len(W))
        else:
            quo = np.divide(len(W), size)
            rem = np.mod(len(W), size)
            if rem == 0:
                w_idx = np.arange(len(W)).reshape(size, quo)
            else:
                w_idx = -np.ones((size, quo+1))
                w_idx[:, :quo] = np.arange(len(W))[:len(W)-rem].reshape(size, quo)
                [w_idx[:rem, -1]] = [range(len(W)-rem,len(W))]
    else:
        w_idx = None
    w_idx = comm.scatter(w_idx, root=0)
    w_idx = np.delete(w_idx, np.where(w_idx==-1)).astype('i')
    return w_idx

def _rootprint(s, comm):
    if comm.Get_rank() == 0:
        logging.info(s)
    else:
        pass

class _fake_comm():
    def Get_rank(self):
        return 0
    
class NullDevice():
    def write(self, s):
        pass

def lat_cor(ibrava, pbrava, ibravm, pbravm, 
            dim = 3, # dimension of the lattice
            num_sol=1, # number of solutions
            distance='Cauchy', # distance function
            vol_th = 0.1, # volume change threshold
            disp=True, lat_opt_disp = False, 
            save_results=False, filename='lat_cor.hdf5',
            save_log = False, logfilename='lattcorr.log' # save the logging file            
            ):
    # start the timer
    t_start = timer()
    
    # configure the log file
    if save_log:
        logging.basicConfig(format='%(message)s', filename=logfilename,filemode='w',level=logging.INFO)
        cnsl = logging.StreamHandler()
        cnsl.setLevel(logging.INFO)
        fmt = logging.Formatter('%(message)s')
        cnsl.setFormatter(fmt)
        logging.getLogger('').addHandler(cnsl)
    else:
        logging.basicConfig(format='%(message)s', level=logging.INFO)
    
    # configure display
    if not disp:
        original_stdout = sys.stdout
        sys.stdout = NullDevice()
    
    # configure distance functions
    isnew = dist_isnew
    unique = dist_unique 
    
    # Configure for parallel processing
    USE_MPI = False
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        n_proc = comm.Get_size()
        proc_rank = comm.Get_rank()
        _rootprint("Configure the parallelization environment", comm)
        _rootprint("=========================================", comm)
        _rootprint('mpi4py is imported successifully!', comm)
        if n_proc == 1:
            _rootprint('But only 1 MPI process is launched, still processing in serial.', comm)
        else:
            _rootprint('%d MPI processes are launched.' % n_proc, comm)
            USE_MPI = True
    except:
        n_proc = 1
        proc_rank = 0
        comm = _fake_comm()
        msg = "Configure the parallelization environment"+'\n'+"========================================="+'\n'+ 'Cannot import MPI from mpi4py. Process in serial.'
        logging.info(msg)
    _rootprint('\n', comm)
    
    
    _rootprint("Input data",comm)
    _rootprint("==========", comm)
    
    # construct
    Lat_A = BravaisLattice(ibrava, pbrava, N=dim)
    Lat_M = BravaisLattice(ibravm, pbravm, N=dim)
    
    E_A = Lat_A.getBase()
    E_M = Lat_M.getBase()
    
    C_A = Lat_A.getConventionalTrans()
    C_M = Lat_M.getConventionalTrans()
    
    LG_A = Lat_A.getSpecialLatticeGroup()
    LG_M = Lat_M.getSpecialLatticeGroup()
    
    _rootprint(" - Austenite lattice:\n    %s" % Lat_A, comm)
    _rootprint(" - Martensite lattice:\n    %s" % Lat_M, comm)
    _rootprint('\n', comm)
    
    # Determine the list of Hermite Normal Forms
    r = det(E_M)/det(E_A) # ratio of the volumes of the unit cells
    
    # find all the sizes of sublattices corresponding to 
    # volume changes less than the thershold "vol_th"
    vf = np.arange(#
              math.ceil((1-vol_th)*r), # lower bound
              math.floor((1+vol_th)*r)+1,  # upper bound
              dtype='int'
              )#
    if len(vf)==0:
        vf = (np.round(r),)     
            
    hnfs = np.zeros((1,dim**2))
    for i in vf:
        hnfs = np.vstack((hnfs, HermiteNormalForms(i, N=dim)))
    hnfs = hnfs[1:] # remove the all-zero first row
    
    _rootprint('The ratio between the volume of unit cells is %g.' % r, comm)
    _rootprint('The volume change threshold is {:.2%}.'.format(vol_th), comm)
    _rootprint('So, the possible size(s) of austenite sublattice is (are) %s.' % vf, comm)
    _rootprint('There are %d Hermite Normal Forms in total.' % len(hnfs), comm)
    _rootprint('\n', comm)
    
    _rootprint('Looking for %d best lattice correspondence(s).' % num_sol, comm)
    _rootprint('\n', comm)
            
    # Search for the best unit cell for each Hermite Normal Form
    _rootprint('Search over all the sublattices', comm)
    _rootprint('===============================', comm)
    
    # divide work
    if USE_MPI:
        h_idx = _divide_work(hnfs, comm)
        sub_hnfs = hnfs[h_idx]
    else:
        sub_hnfs = hnfs
        h_idx = np.arange(len(hnfs))
    
    # initialize 
    try: 
        h5_filename = 'shift_1_dim_'+str(dim)+'.hdf5'
        data = os.path.join(os.path.dirname(__file__),h5_filename)
        f = h5py.File(data, 'r')
        A = f['shifts'].value
        # add all zero
        f.close() 
    except:
        txt_filename = 'shift_1_dim_'+str(dim)+'.txt'
        data = os.path.join(os.path.dirname(__file__),txt_filename)
        A = np.loadtxt(data, delimiter=',') 
    min_d = 1e6
    dopt = []
    lopt = []
    hopt = []
    for i in xrange(len(sub_hnfs)):
        msg = 'Processor %d: Optimizing HNF No.%d ...' % (proc_rank, h_idx[i]+1)
        #print(msg)
        logging.info(msg)
        
        hdlr = logging.StreamHandler()
        if lat_opt_disp:
            hdlr.setLevel(logging.INFO)
        else:
            hdlr.setLevel(logging.CRITICAL)
        
        # optimize the unit cell for the given HNF
        H = sub_hnfs[i]
        lopt_H,dopt_H,success = lat_opt(dot(E_A, H.reshape(dim,dim)), E_M, dist = distance, nsol=num_sol, disp=lat_opt_disp, rank=proc_rank,nhnf = h_idx[i]+1, A=A, hdlr=hdlr)
        msg = 'Processor %d: HNF No.%d finished! Checking update ... ' %  (proc_rank, h_idx[i]+1)
        #print(msg)
        logging.info(msg)
        
        # update the final solutions within the process
        if success:
            no_update = True
            for iH in xrange(len(dopt_H)):
                d = dopt_H[iH]
                l = lopt_H[iH]
                if d <= min_d: # if d is smaller than the largest existing one
                    no_update = False
                    if len(dopt) == 0:
                        dopt = np.array([d])
                        lopt = np.array([l])
                        hopt = np.array([H])
                    else:
                        # insert the new solution
                        L_new = dot(H.reshape(dim,dim), l.reshape(dim,dim))
                        if isnew(L_new, mat_dot(hopt, lopt), LG_A, LG_M)[0]:                            
                            if len(np.where(dopt<=d)[0]):
                                ins_idx = np.where(dopt<=d)[0][-1]+1
                            else:
                                ins_idx = 0
                            dopt = np.insert(dopt, ins_idx, d)
                            lopt = np.insert(lopt, ins_idx, l, axis=0)
                            hopt = np.insert(hopt, ins_idx, H, axis=0)
                # cut off redundant solutions
                if len(dopt) > num_sol:
                    i_max = num_sol # suppose every one after num_sol can be deleted
                    while i_max < len(dopt) and dopt[i_max] <= dopt[i_max-1]:
                        # delete only if the distance is larger than maximum of stored solutions
                        i_max += 1
                    dopt = dopt[:i_max]
                    lopt = lopt[:i_max]
                    hopt = hopt[:i_max]
                    min_d = dopt[-1]
            if not no_update:
                msg = 'Processor %d: HNF No.%d has updated the solution(s)!' % (proc_rank, h_idx[i]+1)
                #print(msg)
                logging.info(msg)
        else:
            msg = 'Processor %d: HNF No.%d Failed! No solution is found for this HNF.' % (proc_rank, h_idx[i]+1)
            #print(msg)
            logging.info(msg)
        # blank line after each HNF
        # print('\n')
    if USE_MPI:
        msg = 'Processor %d: All HNFs on this processor have been dealt.' % proc_rank
        #print(msg)
        logging.info(msg)
        
        # collect result from all processors
        _rootprint('Collecting results from all processors.', comm)
        dopt = comm.gather(dopt, root=0)
        lopt = comm.gather(lopt, root=0)
        hopt = comm.gather(hopt, root=0)
        if proc_rank==0:
            
            dopt = np.hstack(dopt)
            lopt = np.vstack(lopt)
            hopt = np.vstack(hopt)
            
            eopt = mat_dot(hopt, lopt)
            uni_idx = unique(eopt, LG_A, LG_M)
            
            # sort combined solutions
            sort_idx = np.argsort(dopt[uni_idx]).astype('i')
            dopt = dopt[uni_idx[sort_idx]]
            lopt = lopt[uni_idx[sort_idx]]
            hopt = hopt[uni_idx[sort_idx]]
            
            # delete redundant ones
            if len(uni_idx) > num_sol:
                dopt = dopt[uni_idx[:num_sol]]
                lopt = lopt[uni_idx[:num_sol]]
                hopt = hopt[uni_idx[:num_sol]]
            else:
                i_max = num_sol # suppose every one after num_sol can be deleted
                while i_max < len(dopt) and dopt[i_max] <= dopt[i_max-1]:
                    # delete only if the distance is larger than maximum of stored solutions
                    i_max += 1
                dopt = dopt[:i_max]
                lopt = lopt[:i_max]
                hopt = hopt[:i_max]
            
    # the following is only for main processor
    if proc_rank == 0:    
        # Print and save the results
        _rootprint('\n', comm)
        _rootprint('Print and save the results', comm)
        _rootprint('==========================', comm)
        
        # list of lattice correspondences
        cor_list = np.ones_like(lopt)
        U_list = np.ones_like(lopt)
        lambda_list = np.ones((len(dopt), dim))
        
        # if more than required solutions are found
        if len(dopt) > num_sol:
            _rootprint('WARNING: found %d solutions, which is more than %d!' % (len(dopt), num_sol), comm)
        
        for i in xrange(len(dopt)):
            _rootprint('Solution %d out of %d:' % (i+1, len(dopt)), comm)
            _rootprint('----------------------', comm)
            
            H = hopt[i].reshape(dim,dim)
            L = lopt[i].reshape(dim,dim)
            
            # lattice correspondence: Cor = CA.H.L.inv(CM)
            cor_list[i] = dot(dot(inv(C_A), dot(H, L)), C_M).reshape(dim**2)
            
            _rootprint(' - Lattice correspondence:', comm)
            for j in xrange(dim):
                msg = '    ['
                for k in xrange(dim):
                    msg += '{:>5.2f} '.format(cor_list[i,dim*k+j])
                msg = msg[:-1] + '] ==> [ '
                for k in xrange(dim):
                    msg += '%1d ' % np.eye(dim)[j,k]
                msg +=']'
                _rootprint(msg, comm)
            
            # Cauchy-Born deformation gradient: F.EA.H.L = EM
            F = dot(E_M, inv(dot(E_A, dot(H, L))))
            C = dot(F.T, F)
            [la, V]=eig(C) # spectrum decomposition
            la = np.sqrt(np.real(la))
            U = dot(dot(V, np.diag(la)), inv(V))
            _rootprint(' - Transformation stretch matrix:', comm)
            for j in xrange(dim):
                if j == 0:
                    msg = '   U = [' 
                else:
                    msg = '       ['
                for k in xrange(dim):
                    msg += '{:>9.6f} '.format(U[j, k])
                msg = msg[:-1] + ']'
                _rootprint(msg, comm)            
            U_list[i] = U.reshape(dim**2)
            
            for j in xrange(dim):
                if j == 0:
                    msg = '   H = [' 
                else:
                    msg = '       ['
                for k in xrange(dim):
                    msg += '{:>9.6f} '.format(H[j, k])
                msg = msg[:-1] + ']'
                _rootprint(msg, comm)  
            
            # ordered eigen strains
            lambda_list[i,:] = np.sort(la)
            _rootprint(' - Sorted eigen strains:', comm)
            msg = '    '
            for j in xrange(dim):
                msg += 'lambda_%d = %g, ' % (j+1, lambda_list[i,j])
            msg = msg[:-2]+'.'
            _rootprint(msg, comm)
            
            # distance
            msg = ' - Assigned distance '
            if distance=='Ericksen':
                msg += '(Ericksen distance):'
            if distance == 'strain':
                msg += '(strain):'
            if distance == 'Cauchy':
                msg += '(Cauchy distance)'
            _rootprint(msg, comm)
            
            _rootprint('    dist = %g' % dopt[i], comm)
            _rootprint('\n', comm) 
           
        # save results in HDF5 format
#         if save_results:
#             _rootprint('Saving results in "%s" ...' % filename, comm)
#             try:
#                 os.remove(filename)
#             except OSError:
#                 pass
#             f = h5py.File(filename, 'w')
#             ds = f.create_dataset('ibrava',(1,),dtype='int')
#             ds[...] = ibrava
#             ds = f.create_dataset('ibravm',(1,),dtype='int')
#             ds[...] = ibravm
#             try:
#                 pbrava = np.array(pbrava)
#                 ds = f.create_dataset('pbrava',pbrava.shape,dtype='float')
#                 ds[...] = pbrava
#             except:
#                 try: 
#                     ds = f.create_dataset('pbrava',(1,),dtype='float')
#                     ds[...] = pbrava
#                 except:
#                     pass
#             try:
#                 pbravm = np.array(pbravm)
#                 ds = f.create_dataset('pbravm',pbravm.shape,dtype='float')
#                 ds[...] = pbravm
#             except:
#                 try: 
#                     ds = f.create_dataset('pbravm',(1,),dtype='float')
#                     ds[...] = pbravm
#                 except:
#                     pass
#             ds = f.create_dataset('Hopt',hopt.shape,dtype='int')
#             ds[...] = hopt
#             ds = f.create_dataset('Lopt',lopt.shape,dtype='int')
#             ds[...] = lopt
#             ds = f.create_dataset('dopt',dopt.shape,dtype='float')
#             ds[...] = dopt
#             ds = f.create_dataset('Corr',cor_list.shape,dtype='float')
#             ds[...] = cor_list
#             ds = f.create_dataset('U',U_list.shape,dtype='float')
#             ds[...] = U_list
#             ds = f.create_dataset('lambda',lambda_list.shape,dtype='float')
#             ds[...] = lambda_list
#             f.close()
        
        # All done! 
        t_end = timer()
        _rootprint('All done in %g secs.' % (t_end - t_start), comm)
        
    # resume stdout
    if not disp:
        sys.stdout = original_stdout
        
    if proc_rank  == 0:
        return cor_list, U_list, dopt, lambda_list, lopt
    else:
        return [], [], [], [], []
