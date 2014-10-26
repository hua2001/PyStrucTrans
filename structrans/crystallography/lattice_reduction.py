from structrans.general_imports import *
from math import copysign


def Gram_Schmidt(E):
    d = len(E)  # dimension
    Eo = np.empty_like(E)
    M = np.eye(d, dtype='float')
    B = np.zeros(d)
    for i in xrange(d):
        b = E[:, i].copy()
        for j in xrange(i):
            M[i, j] = np.dot(E[:, i], Eo[:, j]) / B[j]
            b -= M[i, j] * Eo[:, j]
        Eo[:, i] = b.copy()
        B[i] = np.dot(b, b)
    return Eo, M


def Gauss_reduction(E):
    d = len(E)
    Eo, M = Gram_Schmidt(E)
    B = E.copy()
    for i in xrange(d):
        for j in xrange(i):
            m = np.round(M[i, j])
            B[:, i] -= m * B[:, j]
    return B


def projection(i, x, E):
    d = len(x)
    p = np.zeros_like(x)
    Eo, _ = Gram_Schmidt(E)
    for j in xrange(i, d):
        p += Eo[:, j] * np.dot(x, Eo[:, j]) / (la.norm(Eo[:, j]) ** 2)

    return p


def swapping(E):
    d = len(E)
    for i in xrange(d-1):
        p1 = projection(i, E[:, i], E)
        p2 = projection(i, E[:, i+1], E)
        # if la.norm(projection(i, E[:, i], E))**2 > (4.0/3.0)*la.norm(projection(i, E[:, i+1], E))**2:
        if p1.dot(p1) > (4. / 3) * p2.dot(p2):
            return True, i    
    return False, -1


def LLL(E):
    swap = True

    Niter = 0

    B = E.copy()
    while Niter < 1000 and swap:
        B = Gauss_reduction(B).copy()
        swap, i = swapping(B)
        if swap:
            b = B[:, i].copy()
            B[:, i] = B[:, i + 1].copy()
            B[:, i + 1] = b
        Niter += 1

    if Niter == 100:
        return B, False

    if not copysign(1, la.det(B)) == copysign(1, la.det(E)):
        B[:, 0] = -B[:, 0]

    return B