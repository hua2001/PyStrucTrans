import numpy as np
from numpy.linalg import inv, norm
from math import sqrt, acos, cos, tan

def stereo_proj(p, z):
    r'''
    For a given zone axis z in Z^3 of a crystal,
    it calculates the stereoprojection of vector p.
    Noticed that if p denotes a plane of the crystal,
    it should be written in its reciprocal base.
    It always computes the projection from the south pole.
    '''
    R = norm(z)
    if norm(p)>1e-6:
        if p.dot(z)>1e-6:
            phat = R*p/norm(p)
        else:
            phat = -R*p/norm(p)
        c_theta = phat.dot(z)/(R**2)
        theta = acos(c_theta)
        R_ps = R*sqrt(2*(1+c_theta))
    
        return np.array(phat - (R_ps - R/cos(theta/2))*(z+phat)/R_ps)
    else:
        return np.array([0, 0, 0])    

def pole_figure(plist,z):
    
    if z[0]**2+z[1]**2>1e-6:
        z_1 = np.array([-z[1],z[0],0])
        z_1 = z_1/norm(z_1)
        z_2 = np.cross(z, z_1)
        z_2 = z_2/norm(z_2)
    else:
        z_1 = np.array([1,0,0])
        z_2 = np.array([0,1,0])
        
    plist = np.array(plist)
    if len(plist.shape)>2:
        print('The input should be triplets of numbers.')
        return None;
    plist = np.array([plist]) if len(plist.shape)==1 else plist
    if plist.shape[1]!=3:
        print('The input should be triplets of numbers.')
        return None;
    p_pole = []
    for p in plist:
        p_proj = stereo_proj(p, z)
        x = p_proj.dot(z_1)
        y = p_proj.dot(z_2)
        p_pole.append([x, y])
    return np.array(p_pole)
    
    
    
    
    
    
    
    
    






