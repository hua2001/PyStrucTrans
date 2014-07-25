from pystructrans.general_imports import *
import FuncDesigner as FD
import openopt as oo

def M3toM9(M):
    """
    convert a 3 x 3 matrix M to 9 x 9 matrix M' so that
    MN = M'n
    where N is 3 x 3 and n is 9 x 1
    """
    M9 = np.zeros((9, 9))
    M9[0:9:3, 0:9:3] = M[:]
    M9[1:9:3, 1:9:3] = M[:]
    M9[2:9:3, 2:9:3] = M[:]
    return M9

def lat_glp(M):
    """
    global minimization in a box
    """
    dim = len(M)
    l = FD.oovars(dim)  # create 4 variables in a single vector
    # d = np.random.random(9)
    # M = M3toM9(M)

    F = FD.sum(l * FD.dot(M, l))

    l0 = np.zeros(dim)
    l0[0] = 1
    startPoint = {l: l0}
    #
    # # set box-bound domain:
    constraints = [l >= -1, l <= 1]

    # # set some general constraints:
    constraints += [
        (FD.sum(l*l) >= 1.0)(tol=1.e-6),
        (FD.sum(l*l) <= 1.0)(tol=1.e-6)
    ]
    #
    # # choose required objective function tolerance:
    # # |f-f*| < fTol, where f* is objective function value in optimal point
    fTol = 1.e-3
    #
    solver = 'interalg'
    # # another global solver to compare (it cannot handle required tolerance fTol)
    # # solver=oosolver('de', iprint=10, maxFunEvals = 10000, maxIter = 1500)
    # # or this solver with some non-default parameters:
    # #solver=oosolver('interalg', fStart = 5.56, maxIter = 1000,maxNodes = 1000000, maxActiveNodes = 15)
    #
    p = oo.GLP(F, startPoint, fTol=fTol, constraints=constraints, dataHandling="raw")
    r = p.minimize(solver, iprint=10)
    return r(l)

if __name__ == '__main__':
    M = np.random.rand(4, 4)
    while la.det(M) <= 0:
        M = np.random.rand(4, 4)
    M = M.T.dot(M)
    v = lat_glp(M)
    print("mineig = {:g}".format(min(la.eig(M)[0])))
    print("sol = {:g}".format(np.dot(v, M.dot(v))))
    print("univector: {:g}".format(v.dot(v)))

