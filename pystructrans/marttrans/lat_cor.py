from __future__ import print_function

import numpy as np
import numpy.linalg as la
import math
import os
import sys
import logging
import copy
import math
from multiprocessing import cpu_count

from timeit import default_timer as timer

from lat_opt import lat_opt
from lat_opt import logger as lo_logger
from dist import dist_isnew, dist_unique
from mat_math import mat_dot
from crystallography import BravaisLattice, HermiteNormalForms
    
class NullDevice():
    def write(self, s):
        pass

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
             pool = Pool(processes=cpu_count())
        then define the mapping function
             def lat_opt_par(args, chuncksize):
                return pool.map(lat_opt_unpack, args, chuncksize)
        finally call lat_cor in __main__
     - logfile: save log into the given file
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
    
    # display
    if disp > 0:
        logging.basicConfig(level=logging.DEBUG)
        # default logger setting
        cnsl = logging.StreamHandler()
        if disp > 1 :
            cnsl.setLevel(logging.DEBUG)
        else:
            cnsl.setLevel(logging.INFO)
        fmt = logging.Formatter('%(message)s')
        cnsl.setFormatter(fmt)
        logger.addHandler(cnsl)
        logger.propagate = False
    if 'logfile' in kwargs:
        logfile = kwargs['logfile']
        fhdlr = logging.FileHandler(logfile, mode='w')
        fhdlr.setLevel(logging.INFO)
        logger.addHandler(fhdlr)
        
    
    ''' 
    ====================
    Preparation - finish
    ====================
    '''
    
    logger.info("Input data")
    logger.info("==========")
    
    # construct
    Lat_A = BravaisLattice(ibrava, pbrava, N=dim)
    Lat_M = BravaisLattice(ibravm, pbravm, N=dim)
    
    E_A = Lat_A.getBase()
    E_M = Lat_M.getBase()
    
    C_A = Lat_A.getConventionalTrans()
    C_M = Lat_M.getConventionalTrans()
    
    LG_A = Lat_A.getSpecialLatticeGroup()
    LG_M = Lat_M.getSpecialLatticeGroup()
    
    logger.info(" - Austenite lattice:")
    logger.info("    {:s}".format(Lat_A))
    logger.info(" - Martensite lattice:")
    logger.info("    {:s}".format(Lat_M))
    
    logger.info("")
    
    # Determine the list of Hermite Normal Forms
    r = la.det(E_M)/la.det(E_A) # ratio of the volumes of the unit cells
    
    # find all the sizes of sublattices corresponding to 
    # volume changes less than the thershold "vol_th"
    vf = np.arange(
              math.ceil((1-vol_th)*r), # lower bound
              math.floor((1+vol_th)*r)+1,  # upper bound
              dtype='int')
    if len(vf)==0:
        vf = (np.round(r),)     
            
    hnfs = np.zeros((1,dim**2))
    for i in vf:
        hnfs = np.vstack((hnfs, HermiteNormalForms(i, N=dim)))
    hnfs = hnfs[1:] # remove the all-zero first row
    
    logger.debug('The ratio between the volume of unit cells is {:g}.'.format(r))
    logger.debug('The volume change threshold is {:.2%}.'.format(vol_th))
    logger.debug('So, the possible size(s) of austenite sublattice is (are) {:s}.'.format(vf))
    logger.debug('There are {:d} Hermite Normal Forms in total.' .format(len(hnfs)))
    logger.debug('Looking for {:d} best lattice correspondence(s).'.format(nsol))
            
    # Search for the best unit cell for each Hermite Normal Form
    logger.info("")
    logger.info('Search over {:d} sublattices'.format(len(hnfs)))
    logger.info('===========================')
    
    def merge_sols(a, b):
        # combine and sort
        tsols = copy.copy(a)
        tsols.extend(b)
        tsols.sort(key=lambda sol: sol['d'])
        # construct unique sols
        usols = copy.copy(tsols[:1])
        for sol in tsols[1:]:
            ls = [us['l'] for us in usols]
            if dist_isnew(sol['l'], ls, LG_A, LG_M)[0]:
                usols.append(sol)
        # cutoff extra solutions
        if len(usols)>nsol:
            while usols[-1]['d'] > usols[nsol]['d']:
                usols.pop()
        return usols
    
    # options for lat_opt
    options = {'nsol': nsol, 'dist': dist}
    lo_level = logging.CRITICAL
    if disp == 2:
        lo_level = logging.WARNING
    if disp == 3:
        lo_level = logging.INFO
    if disp == 4:
        lo_level = logging.DEBUG
    options['level'] = lo_level
    
    if 'lat_opt_par' in kwargs:
        logger.info('HNFs are being solved in parallel ...')
        lat_opt_par = kwargs['lat_opt_par']
        args = [(np.dot(E_A, h.reshape(dim,dim)), E_M, options, ih)
                for ih,h in enumerate(hnfs)]
        chuncksize = len(hnfs)
        map = lat_opt_par(args, 1)
        logger.info('Parallel excution is done.')
        # chain solutions together
        sols = []    
        for i in xrange(len(hnfs)):
#             sols = merge_sols(sols, list({'d':d, 'l':l, 'h':hnfs[i]} for l,d in zip(map[i][0], map[i][1])))
            sols.extend(list({'d':d, 'l':l, 'h':hnfs[i]} for l,d in zip(map[i][0], map[i][1])))
                
    else:
        sols = []  
        logger.info('HNFs are being solved ...')
        for ih,h in enumerate(hnfs):
            logger.debug('Processing the HNF No. {:d}'.format(ih+1))
            res = lat_opt(np.dot(E_A, h.reshape(dim,dim)), E_M, **options)
#             sols = merge_sols(sols, [{'d':d, 'l':l, 'h':hnfs[i]} for l,d in zip(res[0], res[1])])
            sols.extend([{'d':d, 'l':l, 'h':hnfs[i]} for l,d in zip(res[0], res[1])])
        logger.info('Done.')    
    
    logger.debug("") 
    logger.debug("Got {:d} solutions in total.".format(len(sols)))    
    def unique_sols(sols):
        if len(sols) == 1:
            return copy.copy(sols)
        else:
            # partial solutions
            ps = unique_sols(sols[:-1])
            ls = [s['l'] for s in ps]
            if dist_isnew(sols[-1]['l'], ls, LG_A, LG_M)[0]:
                ps.append(sols[-1])
            return copy.copy(ps)
    sols = unique_sols(sols)
    logger.debug("{:d} solutions are not symmetry related.".format(len(sols))) 
    sols.sort(key = lambda s: s['d'])
    
    # cutoff extra solutions
    if len(sols)>nsol:
        while sols[-1]['d'] > sols[nsol]['d']:
            sols.pop()
    logger.debug("Finally {:d} / {:d} solutions are found.".format(len(sols), nsol)) 
    
    # if disp > 0:
    # Print and save the results
    logger.info('\n')
    logger.info('Print and save the results')
    logger.info('==========================')
    
    for i, sol in enumerate(sols):
        logger.info('')
        logger.info('Solution {:d} out of {:d}:'.format(i+1, nsol))
        logger.info('----------------------')
        H = sol['h'].reshape(dim, dim)
        L = sol['l'].reshape(dim, dim)
         
        cor = np.dot(np.dot(la.inv(C_A), np.dot(H, L)), C_M).reshape(dim*dim)
        sol['cor'] = copy.copy(cor)
        logger.info(' - Lattice correspondence:')
        for j in xrange(dim):
            msg = '    ['
            for k in xrange(dim):
                msg += '{:>5.2f} '.format(cor[dim*k+j])
            msg = msg[:-1] + '] ==> [ '
            for k in xrange(dim):
                msg += '%1d ' % np.eye(dim)[j,k]
            msg +=']'
            logger.info(msg)
         
        # Cauchy-Born deformation gradient: F.EA.H.L = EM
        F = np.dot(E_M, la.inv(np.dot(E_A, np.dot(H, L))))
        C = np.dot(F.T, F)
        [lams, V] = la.eig(C) # spectrum decomposition
        lams = np.sqrt(np.real(lams))
        U = np.dot(np.dot(V, np.diag(lams)), la.inv(V))
        logger.info(' - Transformation stretch matrix:')
        for j in xrange(dim):
            if j == 0:
                msg = '   U = [' 
            else:
                msg = '       ['
            for k in xrange(dim):
                msg += '{:>9.6f} '.format(U[j, k])
            msg = msg[:-1] + ']'
            logger.info(msg)            
        sol['U'] = copy.copy(U.reshape(dim*dim))
         
        for j in xrange(dim):
            if j == 0:
                msg = '   H = [' 
            else:
                msg = '       ['
            for k in xrange(dim):
                msg += '{:>9.6f} '.format(H[j, k])
            msg = msg[:-1] + ']'
            logger.info(msg)  
     
        # ordered eigen strains
        lams = np.sort(lams)
        logger.info(' - Sorted eigen strains:')
        msg = '    '
        for j in xrange(dim):
            msg += 'lambda_%d = %g, ' % (j+1, lams[j])
        msg = msg[:-2]+'.'
        logger.info(msg)
        sol['lams'] = copy.copy(lams)
            
        # distance
        msg = ' - Assigned distance '
        if dist=='Ericksen':
            msg += '(Ericksen distance):'
        if dist == 'strain':
            msg += '(strain):'
        if dist == 'Cauchy':
            msg += '(Cauchy distance)'
        logger.info(msg)
          
        logger.info('    dist = {:g}'.format(sol['d']))
        logger.info('')     
    
    logger.info("All done in {:g} secs.".format(timer()-t_start))
    
    return sols
