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
    zhat = z/norm(z)
    Q = -np.eye(3) + 2*np.outer(zhat,zhat)
    if z[0]**2+z[1]**2>1e-6:
        z_1 = np.array([z[1],-z[0],0])
        z_1 = z_1/norm(z_1)
        z_2 = np.cross(z, z_1)
        z_2 = z_2/norm(z_2)
    else:
        z_1 = np.array([1,0,0])
        z_2 = np.array([0,1,0])
        
    
    if norm(p)>1e-6:
        if p.dot(z)>1e-6:
            phat = p/norm(p)
        else:
            phat = -(Q.dot(p))/norm(p)
        c_theta = phat.dot(z)/R
#         print c_theta*180/np.pi
#         R_ps = R*sqrt(2*(1+c_theta))
#         semi_theta = sqrt((c_theta + 1)/2)
        sq_temp = round(1.0 - c_theta**2, 6)
        sq_theta = sqrt(sq_temp)
        if sq_theta<1e-6:
            return np.array([0,0])
        else:
            return np.array([R*sqrt((1-c_theta)/(1+c_theta))*(phat.dot(z_1))/sq_theta,
                             R*sqrt((1-c_theta)/(1+c_theta))*(phat.dot(z_2))/sq_theta ])
#     
#         return np.array(phat - (R_ps - R/semi_theta)*(z+phat)/R_ps)
    else:
        return np.array([0, 0])    

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
#         x = p_proj.dot(z_1)
#         y = p_proj.dot(z_2)
        p_pole.append(p_proj)
    return np.array(p_pole)
    
    
    
    
    
    
    
    
    






