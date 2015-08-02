import logging
from timeit import default_timer as timer
from main_direct import vec_generator as vgen
import numpy as np
import numpy.linalg as la
import multiprocessing as mp

from structrans import BravaisLattice
from structrans.marttrans.lat_cor.quartic_opt import quart_min
from structrans.marttrans.lat_cor.dist import Cauchy_dist, Cauchy_dist_vec


class Solver:
    logger = logging.getLogger(__name__)

    def __init__(self, brava, pa, bravm, pm, **kwargs):
        """
        find the optimal lattice invarient shear move E1 as close as possible to E2
        allowed kwargs:
         - dim: dimension of the Bravais lattice, default is 3, the only other option is 2
         - nsol: number of solutions
         - slice_sz: size of each slice of L's for vectorized calculation, default is 1000
         - disp: level of detail of printing, default is 2
                 0 => no print
                 1 => two lattices and solution
                 2 => progress over HNFs
                 3 => more info in the lat_opt
                 4 => all info in the lat_opt (not implemented)
         - logfile: name of the logfile
         - loglevel: level of the logging
         - mpmgr: multiprocessing manager
        """
        self.brava = brava
        self.pa = pa
        self.bravm = bravm
        self.pm = pm

        # read kwargs
        def readkw(field, default):
            return kwargs[field] if field in kwargs else default

        self.nsol = readkw('nsol', 1)
        self.vol_th = readkw('vol_th', 0.1)
        self.disp = readkw('disp', 1)
        self.dim = readkw('dim', 3)
        self.slice_sz = readkw('slice_sz', 100)
        self.distName = readkw('dist', 'Cauchy')

        def lprint(msg, lev=1):
            # print with level
            if self.disp >= lev:
                print(msg)
            if lev == 1:
                Solver.logger.info(msg)
            elif lev >= 2:
                Solver.logger.debug(msg)

        self.lprint = lprint

        def log_progress(rmax, dmax, checkpoint):
            self.lprint('searching radius = {:>9.6f}, max distance = {:>9.6f}, starting at {:s}'
                        .format(rmax, dmax, str(checkpoint)), 2)

        self.log_progress = log_progress

        # setup logging
        self.loglevel = readkw('loglevel', logging.INFO)
        if 'logfile' in kwargs:
            FORMAT = '[%(levelname)s] %(asctime)-15s %(name)s: %(message)s'
            self.logfile = readkw('logfile', 'untitled.log')
            logging.basicConfig(filename=self.logfile, filemode='w', level=self.loglevel, format=FORMAT)

        self.alpha = None
        self.gamma = None
        self.rmax = None
        self.add_sol = None
        self.update_sol = None

    def load_input(self):
        self.lprint("Input data")
        self.lprint("==========")

        # construct
        Lat_A = BravaisLattice(self.brava, self.pa, dim=self.dim)
        Lat_M = BravaisLattice(self.bravm, self.pm, dim=self.dim)

        E_A = Lat_A.base()
        E_M = Lat_M.base()
        self.E_A = E_A
        self.E_M = E_M

        E_Minv = la.inv(E_M)
        self.dist = lambda x: Cauchy_dist(x, E_A, E_Minv)
        self.dist_vec = lambda xs: Cauchy_dist_vec(xs, E_A, E_Minv)

        C_A = Lat_A.conventional_trans()
        C_M = Lat_M.conventional_trans()

        LG_A = Lat_A.lattice_group().matrices()
        LG_M = Lat_M.lattice_group().matrices()

        def eqlatcor(l):
            L = l.reshape(self.dim, self.dim)
            return set(tuple(M1.dot(L).dot(M2).flatten()) for M1 in LG_A for M2 in LG_M)

        self.lprint(" - Austenite lattice:")
        self.lprint("    {:s}".format(str(Lat_A)))
        self.lprint(" - Martensite lattice:")
        self.lprint("    {:s}".format(str(Lat_M)))
        self.lprint("")

        # minimum gamma
        self.gamma = lambda x: 3 + 2 * np.sqrt(3 * x) + x

        C1 = np.tensordot(E_A, E_Minv, axes=([], [])).transpose(0, 3, 1, 2)
        C2 = np.tensordot(C1, C1, axes=([1], [1])) \
            .transpose(0, 3, 1, 2, 4, 5) \
            .reshape(self.dim ** 2, self.dim ** 2, self.dim ** 2)
        D = np.tensordot(C2, C2, axes=([0], [0]))
        self.lprint('finding alpha by quartic minimization ...', 2)
        self.alpha = quart_min(D)
        self.lprint('alpha = {:>9.6f}'.format(self.alpha), 2)

        # start searching
        def add_sol(sol, sols, ext_sols):
            # if can add to the solution list
            if not tuple(sol['l']) in ext_sols:
                # insert the solution and update EXT_SOLS
                idx = np.searchsorted([s['d'] for s in sols], sol['d'])
                if la.det(sol['l'].reshape(self.dim, self.dim)) != 0:
                    sols.insert(idx, sol)
                ext_sols = ext_sols.union(eqlatcor(sol['l']))
                if len(sols) > self.nsol:
                    while sols[-1]['d'] > sols[self.nsol - 1]['d']:
                        delSol = sols.pop()
                        ext_sols = ext_sols.difference(eqlatcor(delSol['l']))
            return sols, ext_sols, sols[-1]['d']

        self.add_sol = add_sol

        def update_sol(ls, sols, ext_sols, dmax):
            dmax_ref = dmax
            ds = self.dist_vec(ls)
            while len(ds) > 0:
                minidx = np.argmin(ds)
                dmin = ds[minidx]
                ds = np.delete(ds, minidx)
                lmin = ls[minidx].reshape(self.dim ** 2)
                ls = np.delete(ls, minidx, axis=0)
                if dmin <= dmax:
                    sols, ext_sols, dmax = add_sol({'l': lmin, 'd': dmin}, sols, ext_sols)
                else:
                    break
            updated = (dmax < dmax_ref)
            return updated, sols, ext_sols, dmax

        self.update_sol = update_sol

    # noinspection PyTypeChecker,PyUnresolvedReferences
    def initialize(self):
        sols = []
        ext_sols = set()
        dmax = float('inf')

        L0 = np.floor(la.inv(self.E_A).dot(self.E_M)).astype(np.int)

        for _ in xrange(max(20, 5 * self.nsol)):
            shift = np.random.random_integers(-1, 1, (self.dim, self.dim))
            L = L0 + shift
            l = L.reshape(self.dim ** 2)
            while tuple(l) in ext_sols or not (1 - self.vol_th < la.det(L0 + shift) < 1 + self.vol_th):
                shift = np.random.random_integers(-1, 1, (self.dim, self.dim))
                L = L0 + shift
                l = L.reshape(self.dim ** 2)
            sol = {'l': l, 'd': self.dist(l)}
            sols, ext_sols, dmax = self.add_sol(sol, sols, ext_sols)

        self.rmax = lambda d: (self.gamma(d) / self.alpha) ** 0.25

        maxRadius = self.rmax(dmax)
        updated = True
        checkpoint = np.zeros(self.dim ** 2, np.int)

        return sols, ext_sols, dmax
