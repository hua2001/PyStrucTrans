from pystructrans.general_imports import *

import logging, threading
import numpy as np
import numpy.linalg as la
import copy
from itertools import chain

from dist import Eric_dist, eric_dist_mat
from dist import strain_dist, strain_dist_mat
from dist import Cauchy_dist, Cauchy_dist_mat

from pystructrans.crystallography import Lattice, LLL
from pystructrans.mat_math import mat_dot

from timeit import default_timer as timer

# create logger
logger = logging.getLogger(__name__)

def lat_opt_unpack(args):
#     print("args = {:s}".format(args))
    args[2]['ihnf'] = args[3]
    return lat_opt(args[0], args[1], **(args[2]))

def lat_opt(E1, E2, **kwargs):
    '''
    find the optimal lattice invarient shear move E1 as close as possible to E2
    allowed kwargs:
     - nsol: number of solutions
     - dist: choose from "Cauchy", "Ericksen" and "strain", default is "Cauchy"
     - hdlr: logger handler, default is basicConfigure
     - logfile: log file
     - loglevel: log level for the log file
     - disp: display level, 0 - none, 1 - brief, 2 - verbose 
     - maxiter: maximum iteration depth, default is 3
     - ihnf: number of the HNF
    '''
    
    ''' 
    ===========
    Preparation
    =========== 
    '''
    def print_ary(A):
        return "{:s}".format(A.tolist())
    
    # read kwargs
    nsol = kwargs['nsol'] if 'nsol' in kwargs else 1
    maxiter = kwargs['maxiter'] if 'maxiter' in kwargs else 3
    disp = kwargs['disp'] if 'disp' in kwargs else 1
    
    def lprint(msg, lev):
        # print with level
        if disp >= lev: print(msg)
        if lev == 1:
            logging.info(msg)
        elif lev >= 2:
            logging.debug(msg)
    
    if 'logfile' in kwargs:
        logfile = kwargs['logfile']
        fhdlr = logging.FileHandler(logfile, mode='a')
        fhdlr.setLevel(kwargs['loglevel'] if 'loglevel' in kwargs else logging.INFO)
        fhdlr.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(fhdlr)
    
    if 'ihnf' in kwargs:
        lprint("Processing the HNF No. {:d}".format(kwargs['ihnf']+1), 1)
        
    # get dimension of the problem
    dim = len(E1)
    
    # start timer
    t_start = timer()
        
    # LLL-reduction
    Er2 = LLL(E2)
    Er1 = LLL(E1)
    Er2inv = la.inv(Er2)
    
    # starting point
    lo = np.rint(np.dot(la.inv(E1), Er1))
    chi = np.rint(Er2inv.dot(E2))
    
    # determine distance function
    distance = kwargs['dist'] if 'dist' in kwargs else 'Cauchy'
    
    lprint("distance is \"{:s}\"".format(distance), 2)
    
    # distance functions
    if distance == 'Ericksen':
        dist = lambda x: Eric_dist(x, E1, Er2)
        dist_mat = eric_dist_mat
    if distance == 'strain':
        dist = lambda x: strain_dist(x, E1, Er2inv)
        dist_mat = strain_dist_mat
    if distance == 'Cauchy':
        dist = lambda x: Cauchy_dist(x, E1, Er2inv)
        dist_mat = Cauchy_dist_mat
    
    # lattice groups
    LG1 = Lattice(E1).getSpecialLatticeGroup()
    LG2 = Lattice(Er2).getSpecialLatticeGroup()
    SOLG2 = LG2
    
    ''' 
    ====================
    Preparation - finish
    ==================== 
    '''
    
    ''' 
    ===================
    Travesal of GL(n,Z)
    ===================
    '''
    class GLTree ():
        "Tree structure of GL(n,Z) group"   
        
        # all 12 transvectives in the form of a 6x2 array of matrices
        _EYE = np.eye(dim, dtype="float")
        _T = np.array([
            _EYE + np.tensordot(_EYE[i], _EYE[j], axes=0) 
             for i in xrange(dim) 
             for j in xrange(dim) 
             if i!=j ] + [
            _EYE - np.tensordot(_EYE[i], _EYE[j], axes=0) 
             for i in xrange(dim) 
             for j in xrange(dim) 
             if i!=j ])
        CACHE = {}
        
        def __init__(self, elem):
            "constructed by node element"
            if print_ary(elem) in GLTree.CACHE:
                existing = GLTree.CACHE[tuple(elem.flatten())]
                self.copy(existing)
                self.cached = True
            else:
                self.elem = elem
                self.elem_dist = self.calc_elem_dist()
                self.children = self.generate_children()
                GLTree.CACHE[tuple(elem.flatten())] = self
                self.cached = False
        
        def generate_children(self):
            "generate children nodes"
            return [self.elem.dot(t) for i, t in enumerate(GLTree._T)]
        
        def calc_elem_dist(self):
            "distance value of the node element"
            return dist(self.elem)
         
        def calc_children_dist(self):
            "distance values of childrens as list"
            res = []
            for c in self.children:
                if print_ary(c) in GLTree.CACHE:
                    res.append(GLTree.CACHE[tuple(c.flatten())].elem_dist)
                else:
                    res.append(dist(c))
            return res
         
        def is_min(self):
            if self.children_dist is None:
                children_dist = self.calc_elem_dist()
            return all(children_dist >= self.elem_dist)
        
        def copy(self, that):
            self.elem = that.elem
            self.elem_dist = that.elem_dist
            self.children = that.children
            
        def __str__(self):
            return "GLTree - "+print_ary(self.elem)
        
    ''' 
    ===========================
    Tree datastructure - finish
    =========================== 
    '''
        
    # for each new solution, store all symmetry related l's in memory
    EXT_SOLS = {}
    def ext_sols(l):
        ls = [q1.dot(l).dot(q2) for q1 in LG1 for q2 in SOLG2]
        for c in ls:
            EXT_SOLS[tuple(c.flatten())] = True
        
    # update strategy
    def update_solutions(ds, ls, tree):
        # otherwise, update existing solutions
        if len(ds) >= nsol and tuple(tree.elem.flatten()) not in EXT_SOLS:
            for i, d in enumerate(ds):
                if tree.elem_dist <= d:
                    ds.insert(i, tree.elem_dist)
                    ls.insert(i, tree.elem.flatten())
                    ext_sols(tree.elem)
                    # TODO: truncation strategy
                    # delete last one as long as it is not the same as the supposed tail
                    while len(ds) > nsol and abs(ds[-1]-ds[nsol-1])>1E-10:
                        ds = ds[:-1]
                        ls = ls[:-1]
                    return ds[:], ls[:]
        # directly append if not full
        elif tuple(tree.elem.flatten()) not in EXT_SOLS:
            ds.append(tree.elem_dist)
            ls.append(tree.elem.flatten())
            ext_sols(tree.elem)
            return ds[:], ls[:]
        return None, None
           
    ''' 
    ==============
    Core procedure
    ============== 
    '''
    roots = [GLTree(lo)]    
    dopt = []
    lopt = []
     
    depth = 0
    dopt,lopt = update_solutions([],[],roots[0])
    lprint("loop starts with the first trial {:s} => {:g}".format(print_ary(lopt[0]), dopt[0]), 2)
    updated = True;
    while updated and depth < maxiter:
        # update roots, generator first
        t0 = timer()
        # clear update flag
        updated = False;
        # change depth tracker
        depth += 1
        # going down on the tree by iteration    
        new_roots = []
        uncached = 0
        for root in roots:
            for elem in root.children:
                if not tuple(elem.flatten()) in GLTree.CACHE:
                    t = GLTree(elem)    
                    new_roots.append(t)
                    # try to update the solution if the node element is closer than the last solution
                    if (len(dopt) < nsol or t.elem_dist <= dopt[-1]):
                        ds, ls = update_solutions(dopt, lopt, t)
                        if ds:
                            dopt, lopt = ds, ls
                            updated = True
        roots = new_roots
        # debug messages
        update_msg = "found update" if updated else "no update" 
        lprint("number of roots at depth {:d} is {:d}, construction time: {:g}, {:s}.".format(depth, len(roots), timer()-t0, update_msg), 2)
    
    if depth == maxiter and updated:
        lprint("WARNING: maximum depth {:d} reached before solutions guaranteed".format(maxiter), 2)
     
    ''' 
    =======================
    Core procedure - finish
    ======================= 
    '''
     
    # finish timer
    t_elapsed = timer() - t_start
    
    # change back to the original E2 
    lopt = mat_dot(np.array(lopt), np.vstack([chi.flatten()]*len(lopt)))
    lprint("{:d} / {:d} solution(s) found by {:d} iterations and in {:g} sec.".format(len(dopt), nsol, depth, t_elapsed), 2)
    return lopt, dopt