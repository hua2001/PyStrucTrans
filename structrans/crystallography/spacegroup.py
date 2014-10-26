r'''
This files provides methods of calculating the space group symmetry operations.
i.e. skew and glide symmetry.
Meanwhile, for a definite symmetry operation, 
it generates all symmetry related atomic positions for a given generalized equivalent position (gep).
'''
from structrans.general_imports import *
from numpy.linalg import norm
from structrans.util import rotation

def skew(n, p, lattvec):
    '''
    It generates a rotation and translation vector
    with order of rotation n 
    along with lattice vector lattvec in terms of miller index
    '''
    T = lattvec if isinstance(lattvec, np.ndarray) else np.array(lattvec)
    t = T/norm(T)
    theta = 2*np.pi/n
    if p>n:
        p = p % n
    Q = rotation(theta, t, unit="rad")
    c = (p/n)*T
    return Q, c

