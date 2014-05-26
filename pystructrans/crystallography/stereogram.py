import numpy as np
from numpy.linalg import inv, norm
from math import sqrt, acos, cos, tan

def stereo_proj(z, p):
    r'''
    For a given zone axis z in Z^3 of a crystal,
    it calculates the stereoprojection of vector p.
    Noticed that if p denotes a plane of the crystal,
    it should be written in its reciprocal base.
    It always computes the projection from the south pole.
    '''
    R = norm(z)
    if p.dot(z)>1e-6:
        phat = R*p/norm(p)
    else:
        phat = -R*p/norm(p)
    c_theta = phat.dot(z)/(R**2)
    theta = acos(c_theta)
    R_ps = R*sqrt(2*(1+c_theta))
    
    return np.array(phat - (R_ps - R/cos(theta/2))*(z+phat)/R_ps)








