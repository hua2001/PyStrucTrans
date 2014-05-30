import numpy as np
from numpy.linalg import inv, det, eig, norm
from math import sqrt, cos, sin


def direc_trans(L_cor, n_list):
    # convert input list to ndarray
    n_list = np.array(n_list)
    # assume not a 1d-array
    oned = False
    
    if len(n_list.shape) < 2:   # if it is a 1d-array
        # convert to 2d-array
        n_list = n_list.reshape(1,len(n_list))
        # makr the flag. 
        # we want to return an array having the same
        # dimension as the input 
        oned = True
    
    # L_cor * new_n = n ==> new_n  = L_cor^-1 * n    
    new_n_list = [np.dot(inv(L_cor),n_list[i]) for i in xrange(len(n_list))]
    
    if oned:
        return new_n_list[0]
    else:
        return np.array(new_n_list)

def plane_trans(L_cor, n_list):
    n_list = np.array(n_list)
    oned = False
    if len(n_list.shape) < 2:
        n_list = n_list.reshape(1,len(n_list))
        oned = True
    new_n_list = [np.dot(L_cor.T, n_list[i]) for i in xrange(len(n_list))]
    if oned:
        return new_n_list[0]
    else:
        return np.array(new_n_list)
    
def rad(t):
    return t*np.pi/180
def rZ(t):
    return np.array([[cos(t), -sin(t), 0], 
                     [sin(t), cos(t), 0],
                     [0, 0, 1]])
def rX(t):
    return np.array([[1, 0, 0], 
                     [0, cos(t), -sin(t)],
                     [0, sin(t), cos(t)]])
def rY(t):
    return np.array([[cos(t), 0, sin(t)],
                     [0, 1, 0], 
                     [-sin(t), 0, cos(t)]])

def euler(phi_1, Psi, phi_2):
    r'''
    calculate rotation matrix with
    given three Euler angles by
    order Z-X-Z 
    '''
    phi_1 = rad(phi_1)
    Psi = rad(Psi)
    phi_2 = rad(phi_2)
    return np.dot(rZ(phi_2), rY(Psi).dot(rX(phi_1)))
def Rot(t, z):
    r'''
    calculate the rotation matrix about axis z
    '''
    z = np.array(z)
    z = z/norm(z)
    t = rad(t)
    if z[0]**2+z[1]**2<1e-6:
        z2 = np.array([0,1,0])
        z1 = np.array([1,0,0])
    else:
        z2 = np.array([-z[1], z[0], 0])
        z2 = z2/norm(z2)
        z1 = np.cross(z, z2)
    return cos(t)*(np.outer(z1, z1)+np.outer(z2,z2))+sin(t)*(np.outer(z2,z1)-np.outer(z1,z2))+np.outer(z, z)
    



            
        
    
    
        
    
