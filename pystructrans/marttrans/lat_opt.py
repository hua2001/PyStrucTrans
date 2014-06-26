from __future__ import print_function

import logging, threading
import numpy as np
import numpy.linalg as la
import copy
from itertools import chain

from dist import Eric_dist, eric_dist_mat, eric_dist_unique
from dist import strain_dist, strain_dist_mat
from dist import Cauchy_dist, Cauchy_dist_mat
from dist import dist_isnew

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
     - level: logging level 
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
    
    if 'hdlr' in kwargs:
        hdlr = kwargs['hdlr']
        logger.propagate = False
        logger.addHandler(hdlr)    
    else:
        lev = kwargs['level'] if 'level' in kwargs else logging.INFO 
        logger.setLevel(lev)   
    
    if 'logfile' in kwargs:
        logfile = kwargs['logfile']
        fhdlr = logging.FileHandler(logfile, mode='a')
        fhdlr.setLevel(kwargs['loglevel'] if 'loglevel' in kwargs else logging.INFO)
        fhdlr.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(fhdlr)
    
    if 'ihnf' in kwargs:
        logger.warning("Processing the HNF No. {:d}".format(kwargs['ihnf']+1))
        
    # get dimension of the problem
    dim = len(E1) 
    logger.debug("E1 = "+print_ary(E1))
    logger.debug("E2 = "+print_ary(E2))
    logger.debug("dimemsion = {:d}".format(dim))
    
    # start timer
    t_start = timer()
    logger.debug("timer starts")
        
    # LLL-reduction
    Er2 = LLL(E2)
    Er1 = LLL(E1)
    logger.debug("LLL-reduced E1 = "+print_ary(Er1))
    logger.debug("LLL-reduced E2 = "+print_ary(Er2))
    
    # starting point
    lo = np.rint(np.dot(la.inv(E1), Er1))
    chi = np.rint(la.inv(Er2).dot(E2))
    logger.debug("starting point = "+print_ary(lo))
    logger.debug("E2 recovery matrix = "+print_ary(chi))
    
    # determine distance function
    distance = kwargs['dist'] if 'dist' in kwargs else 'Cauchy'
    logger.debug("distance is \"{:s}\"".format(distance))
    if distance == 'Ericksen':
        dist = eric_dist
        dist_mat = eric_dist_mat
        dist_unique = eric_dist_unique
    if distance == 'strain':
        dist = strain_dist
        dist_mat = strain_dist_mat
        dist_unique = eric_dist_unique
    if distance == 'Cauchy':
        dist = Cauchy_dist
        dist_mat = Cauchy_dist_mat
        dist_unique = eric_dist_unique
    
    # lattice groups
    LG1 = Lattice(Er1).getSpecialLatticeGroup().astype('float')
    LG2 = Lattice(Er1).getSpecialLatticeGroup().astype('float') 
    SOLG2 = LG2
    
    ''' 
    ====================
    Preparation - finish
    ==================== 
    '''
    
    ''' 
    ==================
    Tree datastructure
    ================== 
    '''
    class GLTree:
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
        
        def __init__(self, elem, generator=-2, parent=None):
            "constructed by node element and the generator of getting this node from its parent"
            if print_ary(elem) in GLTree.CACHE:
                existing = GLTree.CACHE[print_ary(elem)]
                self.copy(existing)
                self.cached = True
            else:
                self.elem = elem
                self.elem_dist = self.calc_elem_dist()
                self.parent = parent
                self.generator = generator
                self.children = self.generate_children()
                self.children_dist = self.calc_children_dist
                GLTree.CACHE[print_ary(elem)] = self
                self.cached = False
        
        def generate_children(self):
            "generate children nodes"
            # reverse generator, by default is -1   
            rgen = self.generator + 1 if (self.generator % 2 == 0) else self.generator - 1 
            return [((self.elem).dot(t), i) for i,t in enumerate(GLTree._T) if not i == rgen]
        
        def calc_elem_dist(self):
            "distance value of the node element"
            return dist(self.elem, E1, Er2)
         
        def calc_children_dist(self):
            "distance values of childrens as list"
            return [dist(c, E1, Er2) for c,_ in self.children]
         
        def is_min(self):
            if self.generator == -2:
                return all(self.chil_dist()>=self.elem_dist())
            else:
                # restore parent
                parent = self.elem.dot(GLTree._T[self.r_generator])
                return all(self.chil_dist()>=self.elem_dist()) and dist(parent, E1, Er2)
        
        def copy(self, that):
            self.elem = that.elem
            self.elem_dist = that.elem_dist
            self.parent = that.parent
            self.generator = that.generator
            self.children = that.children
            self.children_dist = that.children_dist
            
        def __str__(self):
            return "GLTree - "+print_ary(self.elem)
        
        
    ''' 
    ===========================
    Tree datastructure - finish
    =========================== 
    '''
           
    ''' 
    ==============
    Core procedure
    ============== 
    '''
    lo, _ = GLTree(lo).children[3]
    GLTree.CACHE = {}
    roots = [GLTree(lo)]    
    dopt = []
    lopt = []
     
    # update strategy
    def update_solutions(ds, ls, tree):
        # otherwise, update existing solutions
        for i, d in enumerate(ds):
            if tree.elem_dist <= d and dist_isnew(tree.elem.flatten(), ls, LG1, SOLG2)[0]:
                ds.insert(i, tree.elem_dist)
                ls.insert(i, tree.elem.flatten())
                # TODO: truncation strategy
                # delete last one as long as it is not the same as the supposed tail
                while len(ds) > nsol and abs(ds[-1]-ds[nsol-1])>1E-10:
                    ds = ds[:-1]
                    ls = ls[:-1]
                return copy.copy(ds), copy.copy(ls)
        # directly append if not full
        if len(ds) < nsol:
            ds.append(tree.elem_dist)
            ls.append(tree.elem.flatten())
            return copy.copy(ds), copy.copy(ls)
        return None, None
     
    depth = 0
    dopt,lopt = update_solutions([],[],roots[0])
    logger.debug("loop starts with the first trial {:s} => {:g}".format(print_ary(lopt[0]), dopt[0]))
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
        for root in roots:
            for elem, gen in root.children:
                t = GLTree(elem, generator=gen, parent=root)
                # try to update the solution if the node element is colser than the last solution
                if not t.cached and (len(dopt) < nsol or t.elem_dist <= dopt[-1]):
                    ds, ls = update_solutions(dopt, lopt, t)
                    if ds:
                        logger.debug("solutions updated by {:s} => {:g}".format(print_ary(t.elem), t.elem_dist))
                        dopt, lopt = ds, ls
                        updated = True
                new_roots.append(t)
        roots = new_roots
        logger.debug("number of roots at depth {:d} is {:d}, construction time: {:g}".format(depth, len(roots), timer()-t0))
    
    if depth == maxiter and updated:
        logger.debug("WARNING: maximum depth {:d} reached before solutions gauranteed".format(maxiter))
     
    ''' 
    =======================
    Core procedure - finish
    ======================= 
    '''
     
    # finish timer
    t_elapsed = timer() - t_start
    
    # change back to the original E2 
    lopt = mat_dot(np.array(lopt), np.vstack([chi.flatten()]*len(lopt)))
    logger.info("{:d} / {:d} solution(s) found by {:d} iterations and in {:g} sec.".format(len(dopt), nsol, depth, t_elapsed))
    return lopt, dopt
    
