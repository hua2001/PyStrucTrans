from pystructrans.general_imports import *
import os
import logging
import copy
import math

from timeit import default_timer as timer

from lat_opt import lat_opt
from pystructrans.crystallography import BravaisLattice, HermiteNormalForms
    
# create logger
logger = logging.getLogger(__name__)

def lat_cor(ibrava, pbrava, ibravm, pbravm, **kwargs):
    '''
    find the optimal lattice invarient shear move E1 as close as possible to E2
    allowed kwargs:
     - dim: dimension of the Bravais lattice, default is 3, the only other option is 2
     - nsol: number of solutions
     - vol_th: volume change threshold
     - dist: choose from "Cauchy", "Ericksen" and "strain", default is "Cauchy"
     - disp: level of detail of printing, default is 2
             0 => no print
             1 => two lattices and solution
             2 => progress over HNFs
             3 => basic info in the lat_opt
             4 => all info in the lat_opt
     - hdlr: logger handler
     - maxiter: maximum iteration depth, default is 3
     - lat_opt_par: a map function for parallel excution of lat_opt, for example
        define a pool in __main__ first:
             pool = multiprocessing.Pool()
        then define the mapping function
             def lat_opt_par(args, chuncksize):
                return pool.map(lat_opt_unpack, args, chuncksize)
        finally call lat_cor in __main__
     - logfile: save log into the given file
     - loglevel: log level of the log file
    '''
    
    ''' 
    ===========
    Preparation
    =========== 
    '''
    # start the timer
    t_start = timer()
    
    # read kwargs
    def readkw(field, default):
        return kwargs[field] if field in kwargs else default
    
    nsol = readkw('nsol', 1)
    vol_th = readkw('vol_th', 0.1)
    disp = readkw('disp', 1)
    dim = readkw('dim', 3)
    dist = readkw('dist', 'Cauchy')
    maxiter = readkw('maxiter', 3)
    
    def lprint(msg, lev):
        # print with level
        if disp >= lev: print(msg)
        if lev == 1:
            logging.info(msg)
        elif lev >= 2:
            logging.debug(msg)
            
    if 'logfile' in kwargs:
        logfile = kwargs['logfile']
        try:
            os.remove(logfile)
        except OSError:
            pass
        fhdlr = logging.FileHandler(logfile, mode='a')
        fhdlr.setLevel(readkw('loglevel', logging.INFO))
        fhdlr.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(fhdlr)        
    
    ''' 
    ====================
    Preparation - finish
    ====================
    '''
        
    ''' 
    ============
    Core process
    ============
    '''
    
    
    lprint("Input data\n==========", 1)
    
    # construct
    Lat_A = BravaisLattice(ibrava, pbrava, N=dim)
    Lat_M = BravaisLattice(ibravm, pbravm, N=dim)
    
    E_A = Lat_A.getBase()
    E_M = Lat_M.getBase()
    
    C_A = Lat_A.getConventionalTrans()
    C_M = Lat_M.getConventionalTrans()
    
    LG_A = Lat_A.getSpecialLatticeGroup()
    LG_M = Lat_M.getSpecialLatticeGroup()
    
    lprint(" - Austenite lattice:", 1)
    lprint("    {:s}".format(Lat_A), 1)
    lprint(" - Martensite lattice:", 1)
    lprint("    {:s}".format(Lat_M), 1)
    lprint("", 1)
    
    # Determine the list of Hermite Normal Forms
    r = la.det(E_M)/la.det(E_A)  # ratio of the volumes of the unit cells
    
    # find all the sizes of sublattices corresponding to 
    # volume changes less than the threshold "vol_th"
    vf = np.arange(
              math.ceil((1-vol_th)*r), # lower bound
              math.floor((1+vol_th)*r)+1,  # upper bound
              dtype='int')
    if len(vf)==0:
        vf = (np.round(r),)     

    hnfs = None
    for i in vf:
        if hnfs == None:
            hnfs = HermiteNormalForms(i, dim)
        else:
            hnfs = np.append(hnfs, HermiteNormalForms(i, dim), axis=0)
    
    lprint('The ratio between the volume of unit cells is {:g}.'.format(r), 2)
    lprint('The volume change threshold is {:.2%}.'.format(vol_th), 2)
    lprint('So, the possible size(s) of austenite sublattice is (are) {:s}.'.format(vf), 2)
    lprint('There are {:d} Hermite Normal Forms in total.' .format(len(hnfs)), 2)
    lprint('Looking for {:d} best lattice correspondence(s).'.format(nsol), 2)
            
    # Search for the best unit cell for each Hermite Normal Form
    lprint('Search over {:d} sublattices'.format(len(hnfs)), 1)
    lprint('===========================', 1)
    
    # options for lat_opt
    options = {'nsol': nsol, 'dist': dist, 'disp': disp - 1, 'maxiter': maxiter, 'SOLG2': LG_M}
    if 'logfile' in kwargs:
        options['logfile'] = kwargs['logfile']
    if 'loglevel' in kwargs:
        options['loglevel'] = kwargs['loglevel']
    
    EXT_SOLS = {}
    def ext_sols(l):
        """add l to extended solution dictionary"""
        ls = (q1.dot(l).dot(q2) for q1 in LG_A for q2 in LG_M)
        for c in ls:
            EXT_SOLS[c.tostring()] = True

    def add_sol(sols, s):
        """add new solution s to the list of solutions"""
        ltot = np.dot(s['h'], s['l'])
        if not ltot.tostring() in EXT_SOLS:
        # if True:
            pos = len(sols)
            for k in xrange(len(sols)):
                if sols[k]['d'] >= s['d']:
                    pos = k
                    break
            sols.insert(pos, s)
            ext_sols(ltot)
            truncate_sols(sols, nsol)

    def truncate_sols(sols, nsol):
        """truncate sols to number of solutions"""
        if len(sols) > nsol:
            t = 0
            for i in xrange(len(sols)):
                if i < len(sols) - nsol and sols[- 1 - i]['d'] > sols[nsol - 1]['d']:
                    t += 1
                else:
                    break
            for _ in xrange(t):
                sols.pop(-1)

    def merge_sols(sols, new):
        """merge new solutions into old ones"""
        if len(sols) == 0:
            for s in new:
                add_sol(sols, s)
        else:
            for s in xrange(len(new)):
                if new[s]['d'] <= sols[-1]['d']:
                    add_sol(sols, new[s])
                else:
                    return

    # chunck hnfs into groups if there are too many
    def divide_work(W, size): # divide W into sub-groups of at most size
        ngrp = int(np.ceil(len(W)/float(size)))
        return [W[int(g*size) : min(int((g+1)*size), len(W))] for g in xrange(ngrp)]
    # each group has at most 500 solutions
    # grp_sz = 1000.0/nsol
    grp_sz = len(hnfs)
    job_grps = divide_work(range(len(hnfs)), grp_sz)
    sols = []  
    for ig, g in enumerate(job_grps):
        job = [(i, hnfs[i]) for i in g]
        if 'lat_opt_par' in kwargs:
            # parallel excution
            lprint('HNFs are being solved in parallel ...', 1)
            
            lat_opt_par = kwargs['lat_opt_par']
            args = [(np.dot(E_A, h), E_M, options, ih) for ih, h in job]
            async_results = [lat_opt_par([arg]) for arg in args]
            for ir, ares in enumerate(async_results):
                res = ares.get()
                merge_sols(sols, [{'d': d, 'l': l, 'h': job[ir][1]} for l, d in zip(res['lopt'], res['dopt'])])
                
            lprint("\n{:d}/{:d} job groups finished\n".format(ig+1, len(job_grps)), 2)
        else:
            # sequential exuction
            lprint('HNFs are being solved ...', 1)
            for ih, h in job:
                lprint('Processing the HNF No. {:d}'.format(ih+1), 2)
                res = lat_opt(np.dot(E_A, h), E_M, **options)
                merge_sols(sols, [{'d': d, 'l': l, 'h': hnfs[ih]} for l, d in zip(res['lopt'], res['dopt'])])
            lprint('Done.', 1)   
        
    lprint("\nFinally {:d} / {:d} solutions are found.".format(len(sols), nsol), 3) 
    
    ''' 
    =====================
    Core process - finish
    =====================
    '''

    ''' 
    ============
    Print result
    ============
    '''
    # Print and save the results
    lprint('\nPrint and save the results\n==========================', 1)
    
    for i, sol in enumerate(sols):
        lprint('\nSolution {:d} out of {:d}:'.format(i+1, nsol), 1)
        lprint('----------------------', 1)
        H = sol['h']
        L = sol['l']
         
        cor = np.dot(np.dot(la.inv(C_A), np.dot(H, L)), C_M)
        sol['cor'] = copy.copy(cor)
        lprint(' - Lattice correspondence:', 1)
        for j in xrange(dim):
            msg = '    ['
            for k in xrange(dim):
                msg += '{:>5.2f} '.format(cor[k ,j])
            msg = msg[:-1] + '] ==> [ '
            for k in xrange(dim):
                msg += '%1d ' % np.eye(dim)[j, k]
            msg +=']'
            lprint(msg, 1)
         
        # Cauchy-Born deformation gradient: F.EA.H.L = EM
        F = np.dot(E_M, la.inv(np.dot(E_A, np.dot(H, L))))
        C = np.dot(F.T, F)
        [lams, V] = la.eig(C) # spectrum decomposition
        lams = np.sqrt(np.real(lams))
        U = np.dot(np.dot(V, np.diag(lams)), la.inv(V))
        lprint(' - Transformation stretch matrix:', 1)
        for j in xrange(dim):
            if j == 0:
                msg = '   U = [' 
            else:
                msg = '       ['
            for k in xrange(dim):
                msg += '{:>9.6f} '.format(U[j, k])
            msg = msg[:-1] + ']'
            lprint(msg, 1)            
        sol['U'] = copy.copy(U)

        lprint(' - Hermite normal form:', 1)
        for j in xrange(dim):
            if j == 0:
                msg = '   H = [' 
            else:
                msg = '       ['
            for k in xrange(dim):
                msg += '{:>4g} '.format(H[j, k])
            msg = msg[:-1] + ']'
            lprint(msg, 1)  
     
        # ordered eigen strains
        lams = np.sort(lams)
        lprint(' - Sorted eigen strains:', 1)
        msg = '    '
        for j in xrange(dim):
            msg += 'lambda_%d = %g, ' % (j+1, lams[j])
        msg = msg[:-2]+'.'
        lprint(msg, 1)
        sol['lams'] = copy.copy(lams)
            
        # distance
        msg = ' - Assigned distance '
        if dist=='Ericksen':
            msg += '(Ericksen distance):'
        if dist == 'strain':
            msg += '(strain):'
        if dist == 'Cauchy':
            msg += '(Cauchy distance)'
        lprint(msg, 1)
          
        lprint('    dist = {:g}'.format(sol['d']), 1)
        lprint('', 1)
    ''' 
    =====================
    Print result - finish
    =====================
    '''

    # close log file
    if 'logfile' in kwargs: fhdlr.close()
    # timer
    lprint("All done in {:g} secs.".format(timer()-t_start), 1)
    
    return sols
