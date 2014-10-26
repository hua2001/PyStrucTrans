from __future__ import absolute_import
from structrans.general_imports import *

def Eric_dist(x, E1, E2):
    dim = len(E1)
    L = x.reshape(dim, dim)
    g1 = np.dot(E1.T, E1)
    g2 = np.dot(E2.T, E2)
    M = np.dot(L.T, np.dot(g1, L)) - g2
    # normalization factor
    nf = np.trace(np.dot(g2.T, g2))
    return np.trace(np.dot(M.T, M))/nf

def strain_dist(x, E1, E2inv):
    dim = len(E1)
    L = x.reshape(dim,dim)
    # Cinv = E1 L E2inv E2invT LT E1T
    Finv = np.dot(np.dot(E1, L), E2inv)
    Cinv = np.dot(Finv, Finv.T)
    lam = np.real(la.eig(Cinv)[0])
    slam = np.sqrt(lam)
    return np.sum(lam) - 2*np.sum(slam) + 3

def Cauchy_dist(x, E1, E2inv):
    dim = len(E1)
    L = x.reshape(dim, dim)
    # Cinv = E1 L E2inv E2invT LT E1T
    Finv = np.dot(np.dot(E1, L), E2inv)
    Cinv = np.dot(Finv, Finv.T)
    v = (Cinv - np.eye(dim)).reshape(dim**2)
    return np.dot(v.T, v)