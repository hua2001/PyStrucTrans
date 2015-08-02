import logging
import itertools
from timeit import default_timer as timer
from main_direct import vec_generator as vgen
from main_direct import _sols_tostr
import numpy as np
import numpy.linalg as la
import multiprocessing as mp

from structrans import BravaisLattice
from structrans.marttrans.lat_cor.quartic_opt import quart_min
from structrans.marttrans.lat_cor.dist import Cauchy_dist, Cauchy_dist_vec

# create logger
logger = logging.getLogger(__name__)

if __name__ == "__main__":
    # brava = 2
    # pa = 6.1606
    # bravm = 12
    # pm = [4.4580, 5.7684, 40.6980, 86.80]

    brava = 2
    pa = 2
    bravm = 6
    pm = [1.414, 2]

    nsol = 3
    vol_th = 0.1
    disp = 2

    dim = 3
    distName = 'Cauchy'

    ncpu = 4
    batch = 20
    slice_sz = 100

    '''
    ===========
    Preparation
    ===========
    '''
    # start the timer
    t_start = timer()

    def lprint(msg, lev=1):
        # print with level
        if disp >= lev:
            print(msg)
        if lev == 1:
            logger.info(msg)
        elif lev >= 2:
            logger.debug(msg)

    def log_progress(rmax, dmax, checkpoint):
            lprint('searching radius = {:>9.6f}, max distance = {:>9.6f}, starting at {:s}'
                   .format(rmax, dmax, str(checkpoint)), 2)

    # setup logging
    # loglevel = logging.INFO
    # FORMAT = '[%(levelname)s] %(asctime)-15s %(name)s: %(message)s'
    # logfile = readkw('logfile', 'untitled.log')
    # logging.basicConfig(filename=logfile, filemode='w', level=loglevel, format=FORMAT)

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

    lprint("Input data")
    lprint("==========")

    # construct
    Lat_A = BravaisLattice(brava, pa, dim=dim)
    Lat_M = BravaisLattice(bravm, pm, dim=dim)

    E_A = Lat_A.base()
    E_M = Lat_M.base()

    E_Minv = la.inv(E_M)
    dist = lambda x: Cauchy_dist(x, E_A, E_Minv)
    dist_vec = lambda xs: Cauchy_dist_vec(xs, E_A, E_Minv)

    C_A = Lat_A.conventional_trans()
    C_M = Lat_M.conventional_trans()

    LG_A = Lat_A.lattice_group().matrices()
    LG_M = Lat_M.lattice_group().matrices()

    def eqlatcor(l):
        L = l.reshape(dim, dim)
        return set(tuple(M1.dot(L).dot(M2).flatten()) for M1 in LG_A for M2 in LG_M)

    lprint(" - Austenite lattice:")
    lprint("    {:s}".format(str(Lat_A)))
    lprint(" - Martensite lattice:")
    lprint("    {:s}".format(str(Lat_M)))
    lprint("")

    # minimum gamma
    gamma = lambda x: 3 + 2 * np.sqrt(3 * x) + x

    lprint("Searching")
    lprint("=========")

    C1 = np.tensordot(E_A, E_Minv, axes=([], [])).transpose(0, 3, 1, 2)
    C2 = np.tensordot(C1, C1, axes=([1], [1]))\
        .transpose(0, 3, 1, 2, 4, 5)\
        .reshape(dim**2, dim**2, dim**2)
    D = np.tensordot(C2, C2, axes=([0], [0]))
    lprint('finding alpha by quartic minimization ...', 2)
    alpha = quart_min(D)
    lprint('alpha = {:>9.6f}'.format(alpha), 2)

    # start searching
    def add_sol(sol, sols, ext_sols):
        # if can add to the solution list
        if not tuple(sol['l']) in ext_sols:
            # insert the solution and update EXT_SOLS
            idx = np.searchsorted([s['d'] for s in sols], sol['d'])
            if la.det(sol['l'].reshape(dim, dim)) != 0:
                sols.insert(idx, sol)
            ext_sols = ext_sols.union(eqlatcor(sol['l']))
            if len(sols) > nsol:
                while sols[-1]['d'] > sols[nsol-1]['d']:
                    # sols.pop()
                    delSol = sols.pop()
                    ext_sols = ext_sols.difference(eqlatcor(delSol['l']))
        return sols, ext_sols, sols[-1]['d']

    def update_sol(ls, sols, ext_sols, dmax):
        dmax_ref = dmax
        ds = dist_vec(ls)
        while len(ds) > 0:
            minidx = np.argmin(ds)
            dmin = ds[minidx]
            ds = np.delete(ds, minidx)
            lmin = ls[minidx].reshape(dim**2)
            ls = np.delete(ls, minidx, axis=0)
            if dmin <= dmax:
                sols, ext_sols, dmax = add_sol({'l': lmin, 'd': dmin}, sols, ext_sols)
            else:
                break
        updated = (dmax < dmax_ref)
        return updated, sols, ext_sols, dmax

    # initialization
    L0 = np.floor(la.inv(E_A).dot(E_M)).astype(np.int)
    SOLS = []
    EXT_SOLS = set()
    DMAX = float('inf')

    for _ in xrange(max(50, 20 * nsol)):
        shift = np.random.random_integers(-1, 1, (dim, dim))
        L = L0 + shift
        l = L.reshape(dim**2)
        while tuple(l) in EXT_SOLS or not (1 - vol_th < la.det(L0 + shift) < 1 + vol_th):
            shift = np.random.random_integers(-1, 1, (dim, dim))
            L = L0 + shift
            l = L.reshape(dim**2)
        sol = {'l': l, 'd': dist(l)}
        SOLS, EXT_SOLS, DMAX = add_sol(sol, SOLS, EXT_SOLS)

    rmax = lambda d: (gamma(d) / alpha) ** 0.25

    maxRadius = rmax(DMAX)
    lprint('start searching ...')
    checkpoint = np.zeros(dim**2, np.int)
    log_progress(maxRadius, DMAX, checkpoint)

    def calc_dist(ary, dmax, n_visited, n_total):
        ds = dist_vec(ary)
        candidates = None
        dsmin = np.min(ds)
        if dsmin <= dmax:
            dmax = dsmin
            candidates = ary[np.where(ds <= dmax)]
        workder_id = int(mp.current_process()._identity[0])
        msg = '[%s] Processed %d/%g points within rmax=%g using %g secs' % \
              (mp.current_process().name, n_visited + len(ary) * workder_id, n_total, rmax(dmax), timer() - t_start)
        if candidates is not None:
            msg += ' Found %d candidates' % len(candidates)
        print msg
        return candidates

    pool = mp.Pool(ncpu)
    n_visited = 0
    while True:
        updated = False
        g = vgen(maxRadius, dim**2, checkpoint)
        n_total = (32. * (np.pi ** 4) * (maxRadius ** 9)) / 945.
        results = []
        for _ in xrange(ncpu * batch):
            ary = np.array(list(itertools.islice(g, slice_sz)))
            results.append(pool.apply_async(calc_dist, args=(ary, DMAX, n_visited, n_total)))
            n_visited += len(ary)
        try:
            checkpoint = np.array(next(g))
        except StopIteration:
            break

        candidates = None
        for res in results:
            new_cands = res.get()
            if new_cands is not None:
                if candidates is None:
                    candidates = new_cands
                else:
                    candidates = np.append(candidates, new_cands, axis=0)

        if candidates is not None:
            updated, SOLS, EXT_SOLS, DMAX = update_sol(candidates, SOLS, EXT_SOLS, DMAX)

        maxRadius = rmax(DMAX)
        log_progress(maxRadius, DMAX, checkpoint)

    for sol in SOLS:
        eqls = eqlatcor(sol['l'])
        sol['l'] = sol['l'].reshape(dim, dim)
        for l in eqls:
            L = np.array(l).reshape(dim, dim)
            if la.det(L) > 0:
                sol['l'] = L
                break

        cor = np.dot(np.dot(la.inv(C_A), sol['l']), C_M)
        sol['cor'] = cor[:]
        # Cauchy-Born deformation gradient: F.EA.H.L = EM
        F = np.dot(E_M, la.inv(np.dot(E_A, sol['l'])))
        C = np.dot(F.T, F)
        [lams, V] = la.eig(C)  # spectrum decomposition
        lams = np.sqrt(np.real(lams))
        U = np.dot(np.dot(V, np.diag(lams)), la.inv(V))
        sol['u'] = U[:]
        lams = np.sort(lams)
        sol['lams'] = lams[:]

    lprint("Finally {:d} / {:d} solutions are found.".format(len(SOLS), nsol), 3)
    lprint("")

    '''
    ============
    Print result
    ============
    '''
    # Print and save the results
    lprint('Print and save the results')
    lprint('==========================')

    output = _sols_tostr(SOLS, nsol, dim, distName)
    if disp > 0:
        print(output)
    for line in output.split("\n"):
        logger.info(line)
    '''
    =====================
    Print result - finish
    =====================
    '''

    # timer
    lprint("All done in {:g} secs.".format(timer()-t_start), 1)


