from pystructrans.general_imports import *

import logging
from timeit import default_timer as timer

from .dist import Eric_dist, strain_dist, Cauchy_dist
from ..crystallography import Lattice, LLL
# create logger
logger = logging.getLogger(__name__)
# logger = multiprocessing.get_logger()

def lat_opt_unpack(args):
    args[2]['ihnf'] = args[3]
    return lat_opt(args[0], args[1], **(args[2]))

def lat_opt(E1, E2, **kwargs):
    '''
    find the optimal lattice invarient shear move E1 as close as possible to E2
    allowed kwargs:
     - LG2: special lattice group of E2
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
    if 'logfile' in kwargs:
        logging.basicConfig(level=logging.INFO)
        logger.propagate = False
        logfile = kwargs['logfile']
        fhdlr = logging.FileHandler(logfile, mode='a')
        fhdlr.setLevel(kwargs['loglevel'] if 'loglevel' in kwargs else logging.INFO)
        fhdlr.setFormatter(logging.Formatter('%(message)s'))
        logger.addHandler(fhdlr)

    def lprint(msg, lev):
        # print with level
        if disp >= lev:
            print(msg)
        if lev == 1:
            logger.info(msg)
        elif lev >= 2:
            logger.debug(msg)

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
    lo = np.array(np.rint(np.dot(la.inv(E1), Er1)), dtype='int')
    chi = np.array(np.rint(Er2inv.dot(E2)), dtype='int')
    
    # determine distance function
    distance = kwargs['dist'] if 'dist' in kwargs else 'Cauchy'
    
    lprint("distance is \"{:s}\"".format(distance), 2)
    
    # distance functions
    if distance == 'Ericksen':
        dist = lambda x: Eric_dist(x, E1, Er2)
    if distance == 'strain':
        dist = lambda x: strain_dist(x, E1, Er2inv)
    if distance == 'Cauchy':
        dist = lambda x: Cauchy_dist(x, E1, Er2inv)
    
    # lattice groups
    LG1 = Lattice(E1).getSpecialLatticeGroup()
    SOLG2 = kwargs['SOLG2'] if 'SOLG2' in kwargs else Lattice(Er2).getSpecialLatticeGroup()
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
        T1 = np.array([np.eye(dim, dtype="int")
                       + np.tensordot(np.eye(dim, dtype="int")[i], np.eye(dim, dtype="int")[j], axes=0)
                       for i in xrange(dim) for j in xrange(dim) if i != j])
        T2 = np.array([np.eye(dim, dtype="int")
                       - np.tensordot(np.eye(dim, dtype="int")[i], np.eye(dim, dtype="int")[j], axes=0)
                       for i in xrange(dim) for j in xrange(dim) if i != j])
        _T = np.append(T1, T2, axis=0)
        CACHE = {}
        
        def __init__(self, elem):
            "constructed by node element"
            if elem.tostring() in GLTree.CACHE:
                existing = GLTree.CACHE[elem.tostring()]
                self.copy(existing)
                self.cached = True
            else:
                self.elem = elem
                self.hashcode = elem.tostring()
                self.elem_dist = dist(self.elem)
                self.children = self.generate_children()
                self.children_dist = None
                GLTree.CACHE[self.hashcode] = self
                self.cached = False
        
        def generate_children(self):
            "generate children nodes"
            return (self.elem.dot(t) for t in GLTree._T)
         
        def calc_children_dist(self):
            "distance values of childrens as list"
            return (dist(c) for c in self.generate_children())
         
        def is_min(self):
            if self.children_dist is None:
                children_dist = self.calc_elem_dist()
            return all(children_dist >= self.elem_dist)
        
        def copy(self, that):
            self.elem = that.elem
            self.elem_dist = that.elem_dist
            self.children =that.children
            # self.children = self.generate_children()
            self.hashcode = that.hashcode
            
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
            EXT_SOLS[c.tostring()] = True
        
    # update strategy
    def update_solutions(ds, ls, tree):
        # otherwise, update existing solutions
        if len(ds) >= nsol and tree.hashcode not in EXT_SOLS:
            for i, d in enumerate(ds):
                if tree.elem_dist <= d:
                    ds.insert(i, tree.elem_dist)
                    ls.insert(i, tree.elem)
                    ext_sols(tree.elem)
                    # delete last one as long as it is not the same as the supposed tail
                    while len(ds) > nsol and abs(ds[-1]-ds[nsol-1]) > 1.0E-12:
                        ds = ds[:-1]
                        ls = ls[:-1]
                    return ds[:], ls[:]
        # directly append if not full
        elif tree.hashcode not in EXT_SOLS:
            ds.append(tree.elem_dist)
            ls.append(tree.elem)
            ext_sols(tree.elem)
            return ds[:], ls[:]
        return None, None
           
    ''' 
    ==============
    Core procedure
    ============== 
    '''
    roots = [GLTree(lo)]
    depth = 0
    dopt, lopt = update_solutions([], [], roots[0])
    lprint("loop starts with the first trial {:s} => {:g}".format(str(lopt[0]), dopt[0]), 2)
    updated = True
    while updated and depth < maxiter:
        # update roots, generator first
        t0 = timer()
        # clear update flag
        updated = False
        # change depth tracker
        depth += 1
        # going down on the tree by iteration    
        new_roots = []
        for root in roots:
            for elem in root.children:
                if not elem.tostring() in GLTree.CACHE:
                    t = GLTree(elem)
                    new_roots.append(t)
                    # try to update the solution if the node element is closer than the last solution
                    if len(dopt) < nsol or t.elem_dist <= dopt[-1]:
                        ds, ls = update_solutions(dopt, lopt, t)
                        if ds is not None:
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
    lopt = [l.dot(chi) for l in lopt]
    lprint("{:d} / {:d} solution(s) found by {:d} iterations and in {:g} sec.".format(len(dopt), nsol, depth, t_elapsed), 2)
    return {'lopt': lopt, 'dopt': dopt}