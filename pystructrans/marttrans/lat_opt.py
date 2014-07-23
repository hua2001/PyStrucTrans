from pystructrans.general_imports import *
from timeit import default_timer as timer
import json

from .dist import Eric_dist, strain_dist, Cauchy_dist
from ..crystallography import Lattice, LLL
# create logger
logger = logging.getLogger(__name__)
# logger = multiprocessing.get_logger()

def lat_opt_unpack(args):
    args[2]['ihnf'] = args[3]
    return args[3], lat_opt(args[0], args[1], **(args[2]))

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

    def lprint(msg, lev=2):
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
    Er1 = LLL(E1)

    # starting point
    lo = np.array(np.rint(np.dot(la.inv(E1), Er1)), dtype='int')
    # chi = np.array(np.rint(Er2inv.dot(E2)), dtype='int')

    # determine distance function
    distance = kwargs['dist'] if 'dist' in kwargs else 'Cauchy'

    lprint("distance is \"{:s}\"".format(distance), 2)

    # distance functions
    if distance == 'Ericksen':
        dist = lambda x: Eric_dist(x, E1, E2)
    if distance == 'strain':
        dist = lambda x: strain_dist(x, E1, la.inv(E2))
    if distance == 'Cauchy':
        dist = lambda x: Cauchy_dist(x, E1, la.inv(E2))

    # lattice groups
    LG1 = Lattice(E1).getspeciallatticegroup().matrices()
    SOLG2 = kwargs['SOLG2'] if 'SOLG2' in kwargs else Lattice(Er2).getspeciallatticegroup().matrices()
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

        onetotwo = [(i, j) for i in xrange(dim) for j in xrange(dim) if i != j]

        def __init__(self, elem, parent=None, grandpa=None):
            """
            constructed by node element
            """
            self.elem = elem
            self.parent = parent
            self.grandpa = grandpa
            self.elem_dist = dist(self.elem)
            self.children = self.generate_children()
            self.children_dist = None

        def generate_children(self):
            """
            generate children nodes
            """
            nobacktrans = []
            for i in xrange(len(GLTree._T)):
                if dim != 3 or (
                    not GLTree.gobacktoparent(i, self.parent) and
                    not GLTree.gobacktouncle(i, self.parent, self.grandpa)
                ):
                    nobacktrans.append(i)
            return ((j, self.elem.dot(GLTree._T[j])) for j in nobacktrans)

        @classmethod
        def gobacktoparent(cls, child, parent):
            N = dim * (dim - 1)
            return False if parent is None else child == (parent + N) % N

        @classmethod
        def gobacktouncle(cls, child, parent, grandpa):
            if parent is None or grandpa is None:
                return False
            else:
                N = dim * (dim - 1)
                ic, kc = GLTree.onetotwo[child % N]
                jp, kp = GLTree.onetotwo[parent % N]
                ig, jg = GLTree.onetotwo[grandpa % N]
                return ic == ig and jg == jp and kp == kc

        def calc_neighbor_dist(self):
            """distance values of neighbors"""
            nb = (self.elem.dot(t) for t in GLTree._T)
            return (dist(c) for c in nb)

        def steepest_des(self):
            """direction of the steepest descendant"""
            dmin = 1E10
            nbmin = None
            for i, nbdist in enumerate(self.calc_neighbor_dist()):
                if nbdist - self.elem_dist < dmin:
                    nbmin = i
                    dmin = nbdist - self.elem_dist
            return GLTree(self.elem.dot(GLTree._T[nbmin])) if dmin < 0 else None

        def is_min(self):
            return all(list(self.calc_neighbor_dist()) >= self.elem_dist)

        def copy(self, that):
            self.elem = that.elem
            self.elem_dist = that.elem_dist
            self.children = that.generate_children()

        def __str__(self):
            return "GLTree - "+print_ary(self.elem)

    '''
    ============================
    Tree data structure - finish
    ============================
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
        if len(ds) >= nsol and tree.elem.tostring() not in EXT_SOLS:
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
        elif len(ds) < nsol:
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

    # find local minimum
    lprint("Finding the first local minimum ... ")
    minloc = lo
    minnode = GLTree(minloc)
    nextmin = minnode.steepest_des()
    miniter = 0
    while nextmin is not None:
        # print(minnode.elem_dist)
        # print(nextmin.elem_dist)
        minnode = nextmin
        nextmin = minnode.steepest_des()
        miniter += 1
    lprint("Found local min at {:g} after {:d} iterations".format(minnode.elem_dist, miniter))

    # searching for other solutions
    EXPLORED = {}
    # node = GLTree(lo)
    node = minnode
    EXPLORED[node.elem.tostring()] = True
    roots = [node]
    depth = 0
    dopt, lopt = update_solutions([], [], roots[0])
    lprint("loop starts with the first trial {:s} => {:g}".format(str(lopt[0]), dopt[0]), 2)
    updated = True

    # for debug store variants of all appeared elem
    # EXT_ELEMS = {}
    # ls = [q1.dot(lo).dot(q2) for q1 in LG1 for q2 in SOLG2]
    # for c in ls:
    #     EXT_ELEMS[c.tostring()] = True

    while updated and depth < maxiter:
    # while depth < maxiter:
        # update roots, generator first
        t0 = timer()
        # clear update flag
        updated = False
        # change depth tracker
        depth += 1
        # going down on the tree by iteration
        new_roots = []

        # dmin1 = dmin2 = 1E10

        for root in roots:
            for gen, elem in root.children:
                hashcode = elem.tostring()
                if not hashcode in EXPLORED:
                    EXPLORED[hashcode] = True
                    t = GLTree(elem, parent=gen, grandpa=root.parent)
                    new_roots.append(t)

                    # if t.is_min():
                    #     print("found another local min at {:g}".format(t.elem_dist))

                    # if elem.tostring() not in EXT_ELEMS:
                    #     oldmap = EXT_SOLS.copy()
                    #     ls = [q1.dot(elem).dot(q2) for q1 in LG1 for q2 in SOLG2]
                    #     for c in ls:
                    #         EXT_ELEMS[c.tostring()] = True
                    #     dmin1 = min(t.elem_dist, dmin1)
                    #
                    # dmin2 = min(t.elem_dist, dmin2)

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

    del EXT_SOLS
    del EXPLORED

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
    # lopt = [l.dot(chi) for l in lopt]
    lprint("{:d} / {:d} solution(s) found by {:d} iterations and in {:g} sec.".format(len(dopt), nsol, depth, t_elapsed), 2)
    return {'lopt': lopt, 'dopt': dopt}