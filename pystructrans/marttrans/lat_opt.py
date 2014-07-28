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
    global miniter, maxiter
    maxiter = kwargs['maxiter'] if 'maxiter' in kwargs else 2
    miniter = kwargs['miniter'] if 'miniter' in kwargs else 2
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
    SOLG2 = kwargs['SOLG2'] if 'SOLG2' in kwargs else Lattice(E2).getspeciallatticegroup().matrices()
    '''
    ====================
    Preparation - finish
    ====================
    '''

    '''
    ===================
    Travesal of SL(n,Z)
    ===================
    '''
    class SLNode ():
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
            for i in xrange(len(SLNode._T)):
                if dim != 3 or (
                    not SLNode.gobacktoparent(i, self.parent) and
                    not SLNode.gobacktouncle(i, self.parent, self.grandpa) and
                    not SLNode.gobacktosibling(i, self.parent)
                ):
                    nobacktrans.append(i)
            return ((j, self.elem.dot(SLNode._T[j])) for j in nobacktrans)

        @classmethod
        def gobacktoparent(cls, child, parent):
            N = dim * (dim - 1)
            return False if parent is None else child == (parent + N) % N

        @classmethod
        def gobacktouncle(cls, child, parent, grandpa):
            """
            [Tij, Tjk] = Tik => Tij Tjk = Tik Tjk Tij
            """
            if parent is None or grandpa is None:
                return False
            else:
                N = dim * (dim - 1)
                ic, kc = SLNode.onetotwo[child % N]
                jp, kp = SLNode.onetotwo[parent % N]
                ig, jg = SLNode.onetotwo[grandpa % N]
                return ic == ig and jg == jp and kp == kc

        @classmethod
        def gobacktosibling(cls, child, parent):
            """
            [Tij Tkl] = 1 =>
            """
            if parent is None:
                return False
            else:
                N = dim * (dim - 1)
                i, j = SLNode.onetotwo[child % N]
                k, l = SLNode.onetotwo[parent % N]
                return (i, j) > (k, l)

        def calc_neighbor_dist(self):
            """distance values of neighbors"""
            nb = (self.elem.dot(t) for t in SLNode._T)
            return (dist(c) for c in nb)

        def steepest_des(self):
            """direction of the steepest descendant"""
            dmin = 1E10
            nbmin = None
            for i, nbdist in enumerate(self.calc_neighbor_dist()):
                if nbdist - self.elem_dist < dmin:
                    nbmin = i
                    dmin = nbdist - self.elem_dist
            return SLNode(self.elem.dot(SLNode._T[nbmin])) if dmin < 0 else None

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

    '''
    ==============
    Core procedure
    ==============
    '''
    # initialization
    global dopt, lopt, depth
    dopt = []
    lopt = []
    depth = 0

    # DEBUG
    global nlocmin
    nlocmin = 0

    # # look for new minnode
    # minnode_changed = True

    # for each new solution, store all symmetry related l's in memory
    EXT_SOLS = {}

    def ext_sols(l):
        for c in (q1.dot(l).dot(q2) for q1 in LG1 for q2 in SOLG2):
            EXT_SOLS[c.tostring()] = True

    def truncate_sols():
        global dopt, lopt
        # delete last one as long as it is not the same as the supposed tail
        while len(dopt) > nsol and dopt[-1] > dopt[nsol - 1]:
            dopt.pop(-1)
            lopt.pop(-1)

    # update strategy
    def update_solutions(tree):
        global dopt, lopt
        # otherwise, update existing solutions
        if len(dopt) > 0 and tree.elem.tostring() not in EXT_SOLS:
            for i, d in enumerate(dopt):
                if tree.elem_dist <= d:
                    dopt.insert(i, tree.elem_dist)
                    lopt.insert(i, tree.elem)
                    ext_sols(tree.elem)
                    truncate_sols()
                    return True
            if len(dopt) < nsol:
                dopt.append(tree.elem_dist)
                lopt.append(tree.elem)
                ext_sols(tree.elem)
                return True
        # directly append if it is the first one
        elif len(dopt) == 0:
            dopt.append(tree.elem_dist)
            lopt.append(tree.elem)
            ext_sols(tree.elem)
            return True
        return False

    def findlocmin(start):
        # find local minimum
        lprint("Finding the first local minimum ... ")
        minnode = start
        nextmin = minnode.steepest_des()
        minlociter = 0
        while nextmin is not None:
            minnode = nextmin
            nextmin = minnode.steepest_des()
            minlociter += 1
        lprint("Found local min at {:g} after {:d} iterations".format(minnode.elem_dist, minlociter))

        # DEBUG
        global nlocmin
        nlocmin += 1

        return minnode


    def bfs(source):
        """
        breath first search of solutions
        :return: new starting point or None
        """
        global dopt, lopt, depth, miniter, maxiter

        EXPLORED = {}
        # bfs search from minnode
        EXPLORED[source.elem.tostring()] = True
        roots = [source]
        update_solutions(roots[0])
        lprint("loop starts with the first trial {:s} => {:g}".format(str(source.elem.flatten()), source.elem_dist), 2)
        updated = True

        depth = 0
        while depth < miniter or (updated and depth < maxiter):
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
                        t = SLNode(elem, parent=gen, grandpa=root.parent)
                        new_roots.append(t)

                        if t.elem_dist < source.elem_dist:
                            # a new starting point
                            del EXPLORED
                            return t

                        # try to update the solution if the node element is closer than the last solution
                        if update_solutions(t):
                            updated = True

            roots = new_roots
            # debug messages
            update_msg = "found update" if updated else "no update"
            lprint("number of roots at depth {:d} is {:d}, construction time: {:g}, {:s}.".format(depth, len(roots), timer()-t0, update_msg), 2)

        # all done
        del EXPLORED
        return None

    new_start = SLNode(lo)
    while new_start is not None:
        # local minimization
        minnode = findlocmin(new_start)
        # BFS
        new_start = bfs(minnode)

#     del EXT_SOLS

    # if depth == maxiter and updated:
    #     lprint("WARNING: maximum depth {:d} reached before solutions guaranteed".format(maxiter), 2)

    '''
    =======================
    Core procedure - finish
    =======================
    '''
    # finish timer
    t_elapsed = timer() - t_start

    # change back to the original E2
    # lopt = [l.dot(chi) for l in lopt]
    lprint("{:d} / {:d} solution(s) found by {:d} iterations and in {:g} sec.".format(len(dopt), nsol, depth, t_elapsed))
    return {'lopt': lopt, 'dopt': dopt}
