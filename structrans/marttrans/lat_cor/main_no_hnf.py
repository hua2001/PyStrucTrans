from __future__ import absolute_import
from structrans.general_imports import *
import copy
import math
from timeit import default_timer as timer

from structrans import BravaisLattice
from structrans.crystallography import HermiteNormalForms, HNFDecomposition, LLL, Lattice
from structrans.marttrans.lat_cor.quartic_opt import _quart_min
from structrans.marttrans.lat_cor.dist import Cauchy_dist, Eric_dist, strain_dist

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
            n = n + 1
    else:
        firstDigit = x0[0]
        xs = x0[1:]
        while firstDigit <= maxlen:
            g = _vec_generator(np.sqrt(maxlen**2 - firstDigit**2), dim - 1, xs)
            for subvec in g:
                v = subvec[:]
                v.append(firstDigit)
                yield v
                if firstDigit != 0:
                    v = subvec[:]
                    v.append(-firstDigit)
                    yield v
            firstDigit += 1
            xs = None

def lat_cor(brava, pa, bravm, pm, **kwargs):
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
             3 => more info in the lat_opt
             4 => all info in the lat_opt (not implemented)
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
    distName = readkw('dist', 'Cauchy')


    def lprint(msg, lev=1):
        # print with level
        if disp >= lev:
            print(msg)
        if lev == 1:
            logger.info(msg)
        elif lev >= 2:
            logger.debug(msg)

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
    #E_Mr = LLL(E_M)
    #chi = np.array(np.rint(la.inv(E_Mr).dot(E_M)), dtype='int')

    E_Minv = la.inv(E_M)
    distFuncs = {
        'Ericksen': lambda x: Eric_dist(x, E_A, E_M),
        'Cauchy': lambda x: Cauchy_dist(x, E_A, E_Minv),
        'strain': lambda x: strain_dist(x, E_A, E_Minv)
    }
    dist = distFuncs[distName]

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
    
    # Determine the list of Hermite Normal Forms
    r = la.det(E_M)/la.det(E_A)  # ratio of the volumes of the unit cells

    # minimum gamma
    gamma = lambda x: 3 + 2 * np.sqrt(3 * x) + x

    # 9 x 9 matrix
    C1 = np.tensordot(E_A, E_Minv, axes=([], [])).transpose(0, 3, 1, 2)
    C2 = np.tensordot(C1, C1, axes=([1], [1]))\
        .transpose(0, 3, 1, 2, 4, 5)\
        .reshape(dim**2, dim**2, dim**2)
    D = np.tensordot(C2, C2, axes=([0], [0]))
    alpha = _quart_min(D)

    lprint("Searching")
    lprint("=========")

    lprint('alpha = {:>9.6f}'.format(alpha), 2)

    # start searching
    def add_sol(sol, sols, ext_sols):
        # if can add to the solution list
        if not tuple(sol['l']) in ext_sols:
            # insert the solution and update EXT_SOLS
            idx = np.searchsorted([s['d'] for s in sols], sol['d'])
            sols.insert(idx, sol)
            ext_sols = ext_sols.union(eqlatcor(sol['l']))
            if len(sols) > nsol:
                while sols[-1]['d'] > sols[nsol-1]['d']:
                    # sols.pop()
                    delSol = sols.pop()
                    ext_sols = ext_sols.difference(eqlatcor(delSol['l']))
        return sols, ext_sols, sols[-1]['d']

    # initialization
    L0 = np.floor(la.inv(E_A).dot(E_M)).astype(np.int)
    SOLS = []
    EXT_SOLS = set()
    DMAX = float('inf')
    for _ in xrange(nsol):
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
    lprint('searching radius = {:>9.6f}, max distance = {:>9.6f}'.format(maxRadius, DMAX), 2)
    updated = True

    # iteration
    while updated:
        updated = False
        g = _vec_generator(maxRadius, dim**2)
        for l in g:
            l = np.array(l)
            d = dist(l)
            if d <= DMAX:
                oldDMAX = DMAX
                SOLS, EXT_SOLS, DMAX = add_sol({'l': l, 'd': d}, SOLS, EXT_SOLS)
                if oldDMAX > DMAX:
                    maxRadius = rmax(DMAX)
                    lprint('searching radius = {:>9.6f}, max distance = {:>9.6f}'.format(maxRadius, DMAX), 2)
                    updated = True
                    break

    for sol in SOLS:
        sol['l'] = sol['l'].reshape(dim, dim)
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