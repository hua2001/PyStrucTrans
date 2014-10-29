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
    L = x.reshape(dim, dim)
    # F E1 L = E2 => F = E2 (E1 L)^-1 => F^-1 = E1 L E2^-1
    Finv = np.dot(np.dot(E1, L), E2inv)
    Cinv = np.dot(Finv, Finv.T)
    lam = np.real(la.eig(Cinv)[0])
    slam = np.sqrt(lam)
    return np.sum(lam) - 2*np.sum(slam) + 3


def Cauchy_dist(x, E1, E2inv):
    dim = len(E1)
    L = x.reshape(dim, dim)
    # F E1 L = E2 => F = E2 (E1 L)^-1 => F^-1 = E1 L E2^-1
    Finv = np.dot(np.dot(E1, L), E2inv)
    Cinv = np.dot(Finv, Finv.T)
    v = (Cinv - np.eye(dim)).reshape(dim**2)
    return np.dot(v.T, v)


def Cauchy_dist_vec(xs, E1, E2inv):
    if len(xs) < 100:
        return np.array([Cauchy_dist(x, E1, E2inv) for x in xs])
    else:
        dim = len(E1)
        dimsquare = dim**2
        xs = xs.reshape(len(xs), dim, dim)
        F1 = np.tensordot(xs, E2inv, axes=([2], [0]))
        Fts = np.tensordot(F1, E1, axes=([1], [1]))
        Cs = np.diagonal(np.tensordot(Fts, Fts, axes=([1], [1])), axis1=0, axis2=2).transpose(2, 0, 1)
        I = np.eye(dim)
        Es = (Cs - I).reshape(len(xs), dimsquare)
        return np.array([e.dot(e) for e in Es])



