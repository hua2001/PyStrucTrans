from __future__ import absolute_import
from structrans.general_imports import *
from timeit import default_timer as timer
import itertools

from structrans import BravaisLattice
from structrans.marttrans.lat_cor.quartic_opt import _quart_min
from structrans.marttrans.lat_cor.dist import Cauchy_dist, Cauchy_dist_vec

# create logger
logger = logging.getLogger(__name__)

def _vec_generator(maxlen, dim, x0=None):
    """
    generate dim-dimensional vectors with maximum length maxlen
    """
    if x0 is None:
        x0 = [0] * dim
    if dim == 1:
        n = x0[0]
        while n <= maxlen:
            yield [n]
            if n > 0:
                yield [-n]
            n = abs(n) + 1
    else:
        firstDigit = x0[0]
        xs = x0[1:]
        while firstDigit <= maxlen:
            g = _vec_generator(np.sqrt(max(0, maxlen**2 - firstDigit**2)), dim - 1, xs)
            for subvec in g:
                v = [firstDigit] * dim
                v[1:] = subvec[:]
                yield v
            if firstDigit > 0:
                g = _vec_generator(np.sqrt(max(0, maxlen**2 - firstDigit**2)), dim - 1, xs)
                for subvec in g:
                    v = [-firstDigit] * dim
                    v[1:] = subvec[:]
                    yield v
            firstDigit = abs(firstDigit) + 1
            xs = None

def lat_cor(brava, pa, bravm, pm, **kwargs):
    '''
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
    slice_sz = readkw('slice_sz', 100)
    distName = readkw('dist', 'Cauchy')


    def lprint(msg, lev=1):
        # print with level
        if disp >= lev:
            print(msg)
        if lev == 1:
            logger.info(msg)
        elif lev >= 2:
            logger.debug(msg)

    # setup logging
    loglevel =readkw('loglevel', logging.INFO)
    if 'logfile' in kwargs:
        # FORMAT = '[%(levelname)s] %(asctime)-15s %(name)s: %(message)s'
        FORMAT = '%(message)s'
        logging.basicConfig(filename=kwargs['logfile'], filemode='w', level=loglevel, format=FORMAT)

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

    C_A = Lat_A.getConventionalTrans()
    C_M = Lat_M.getConventionalTrans()

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
    alpha = _quart_min(D)
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

    for _ in xrange(20 * nsol):
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
    updated = True
    checkpoint = np.zeros(dim**2)

    # iteration
    while updated:
        updated = False
        g = _vec_generator(maxRadius, dim**2, checkpoint)
        lprint('searching radius = {:>9.6f}, max distance = {:>9.6f}, starting at {:s}'.format(
            maxRadius, DMAX, str(checkpoint)
        ), 2)
        ary = np.array(list(itertools.islice(g, slice_sz)))
        while len(ary) > 0:
            updated, SOLS, EXT_SOLS, DMAX = update_sol(ary, SOLS, EXT_SOLS, DMAX)
            if updated:
                maxRadius = rmax(DMAX)
                checkpoint = np.array(next(g))
                break
            else:
                ary = np.array(list(itertools.islice(g, slice_sz)))

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

    return SOLS

def _sols_tostr(sols, nsol, dim, dist):
    """convert solutions to formatted string so that it can be logged or printed to the screen"""
    output = ""
    for i, sol in enumerate(sols):

        output += 'Solution {:d} out of {:d}:\n'.format(i+1, nsol)
        output += '----------------------\n'

        cor = sol['cor']
        output += ' - Lattice correspondence:\n'
        for j in xrange(dim):
            msg = '    ['
            for k in xrange(dim):
                msg += '{:>5.2f} '.format(cor[k, j])
            msg = msg[:-1] + '] ==> [ '
            for k in xrange(dim):
                msg += '%1d ' % np.eye(dim)[j, k]
            msg += ']'
            output += msg + "\n"

        output += ' - Transformation stretch matrix:\n'
        for j in xrange(dim):
            if j == 0:
                msg = '   U = ['
            else:
                msg = '       ['
            for k in xrange(dim):
                msg += '{:>9.6f} '.format(sol['u'][j, k])
            msg = msg[:-1] + ']'
            output += msg + "\n"


        # ordered eigen strains
        lams = sol['lams']
        output += ' - Sorted eigen strains:\n'
        msg = '    '
        for j in xrange(dim):
            msg += 'lambda_%d = %g, ' % (j+1, lams[j])
        msg = msg[:-2]+'.'
        output += msg + "\n"

        # distance
        msg = ' - Assigned distance '
        if dist == 'Ericksen':
            msg += '(Ericksen distance):'
        if dist == 'strain':
            msg += '(strain):'
        if dist == 'Cauchy':
            msg += '(Cauchy distance)'
        output += msg + "\n"
        output += '    dist = {:g}\n \n'.format(sol['d'])
    return output