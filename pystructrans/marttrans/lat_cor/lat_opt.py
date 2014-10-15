from __future__ import division, absolute_import
from timeit import default_timer as timer

from pystructrans.general_imports import *
from pystructrans.marttrans.lat_cor.dist import Eric_dist, strain_dist, Cauchy_dist
from pystructrans.crystallography import Lattice, LLL
from .dfs import dfs
from .slnode import SLNode

# global SLNode

# create logger
logger = logging.getLogger(__name__)

def lat_opt_unpack(args):
    import multiprocessing
    logger = multiprocessing.get_logger()
    args[2]['ihnf'] = args[3]
    return args[3], lat_opt(args[0], args[1], **(args[2]))

def lat_opt(E1, E2r, **kwargs):
    '''
    find the optimal lattice invarient shear move E1 as close as possible to E2
    allowed kwargs:
     - LG2: special lattice group of E2r
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

    def maxiter_by_nsol(nsol):
        """
        automatically choose maxiter by nsol
        """
        if nsol < 3:
            return 3
        if nsol < 10:
            return 4
        return 5

    # read kwargs
    nsol = kwargs['nsol'] if 'nsol' in kwargs else 1
    global miniter, maxiter
    maxiter = kwargs['maxiter'] if 'maxiter' in kwargs else maxiter_by_nsol(nsol)
    miniter = kwargs['miniter'] if 'miniter' in kwargs else 3
    disp = kwargs['disp'] if 'disp' in kwargs else 1
    if 'logfile' in kwargs:
        logging.basicConfig(level=logging.INFO)
        logger.propagate = False
        logfile = kwargs['logfile']
        fhdlr = logging.FileHandler(logfile, mode='a')
        fhdlr.setLevel(kwargs['loglevel'] if 'loglevel' in kwargs else logging.INFO)
        fhdlr.setFormatter(logging.Formatter(' %(message)s'))
        logger.addHandler(fhdlr)

    if 'ihnf' in kwargs:
        logger.info("Processing the HNF No. {:d}".format(kwargs['ihnf'] + 1))
        if disp > 0:
            print("Processing the HNF No. {:d}".format(kwargs['ihnf'] + 1))

    # get dimension of the problem
    dim = len(E1)
    SLNode.dim = dim

    # start timer
    t_start = timer()

    # LLL-reduction
    Er1 = LLL(E1)
    # starting point
    lo = np.array(np.rint(np.dot(la.inv(E1), Er1)), dtype='int')

    # determine distance function
    distance = kwargs['dist'] if 'dist' in kwargs else 'Cauchy'

    logger.debug("distance is \"{:s}\"".format(distance))

    # distance functions
    if distance == 'Ericksen':
        dist = lambda x: Eric_dist(x, E1, E2r)
    if distance == 'strain':
        E2inv = la.inv(E2r)
        dist = lambda x: strain_dist(x, E1, E2inv)
    if distance == 'Cauchy':
        E2inv = la.inv(E2r)
        dist = lambda x: Cauchy_dist(x, E1, E2inv)

    # lattice groups
    LG1 = Lattice(E1).getspeciallatticegroup().matrices()
    SOLG2 = kwargs['SOLG2'] if 'SOLG2' in kwargs else Lattice(E2r).getspeciallatticegroup().matrices()

    SLNODE_CACHE = {'dist': dist, 'lg1': LG1, 'lg2': SOLG2}
    DFS_CACHE = {'counter': 0}

    '''
    ====================
    Preparation - finish
    ====================
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
        update_solutions(source)
        logger.debug("loop starts with the first trial {:s} => {:g}".format(str(source.elem.flatten()), source.elem_dist))
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

            for root in roots:
                for gen, elem in root.children:
                    hashcode = elem.tostring()
                    if not hashcode in EXPLORED:
                        EXPLORED[hashcode] = True
                        t = SLNode(elem, SLNODE_CACHE, parent=gen, grandpa=root.parent)
                        new_roots.append(t)

                        # if in min iter
                        if depth <= miniter:
                            potential_newstart = dfs(t, DFS_CACHE)
                            if potential_newstart.elem_dist < source.elem_dist:
                                # a new starting point
                                del EXPLORED
                                return potential_newstart

                        # try to update the solution if the node element is closer than the last solution
                        if update_solutions(t):
                            updated = True

            roots = new_roots
            # debug messages
            update_msg = "found update" if updated else "no update"
            logger.debug("number of roots at depth {:d} is {:d}, construction time: {:g}, {:s}.".format(depth, len(roots), timer()-t0, update_msg))

        # all done
        del EXPLORED
        return None

    logger.debug("starting point: {:s}".format(str(lo.flatten())))
    new_start = SLNode(lo, SLNODE_CACHE)
    nrestart = 0

    while new_start is not None:
        # local minimization
        minnode = dfs(new_start, DFS_CACHE)
        # DEBUG
        nrestart += 1
        # BFS
        new_start = bfs(minnode)

    DFS_CACHE = {'counter': 0}
    SLNODE_CACHE = {'dist': dist, 'lg1': LG1, 'lg2': SOLG2}

    logger.debug("restarted {:d} times in total".format(nrestart))

    # if depth == maxiter and updated:
    #     lprint("DEBUG: maximum depth {:d} reached before solutions guaranteed".format(maxiter), 2)

    '''
    =======================
    Core procedure - finish
    =======================
    '''
    # finish timer
    t_elapsed = timer() - t_start

    if 'ihnf' in kwargs:
        logger.debug("HNF No. {:d} finished.".format(kwargs['ihnf']+1))
    for d, l in zip(dopt, lopt):
        logger.debug("{:.4f}, {:s}".format(d, str(l.flatten())))

    final_msg = "{:d} / {:d} solution(s) found by {:d} iterations and in {:g} sec.".format(len(dopt), nsol, depth, t_elapsed)
    if 'ihnf' in kwargs:
        final_msg = "HNF No. {:d}: ".format(kwargs['ihnf']+1) + final_msg
    logger.info(final_msg)
    if disp > 1:
        print(final_msg)
    return {'lopt': lopt, 'dopt': dopt}
