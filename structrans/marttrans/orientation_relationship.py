from ..general_imports import *

def direc_trans(L_cor, n_list):
    """
    transform a list of directions in initial phase
    to final phase based on lattice correspondence

    :param L_cor: lattice correspondence
    :param n_list: list of directions
    """
    # convert input list to ndarray
    n_list = np.array(n_list)
    if len(n_list.shape) < 2:   # if it is a 1d-array
        # convert to 2d-array
        n_list = n_list.reshape(1, len(n_list))

    # L_cor * new_n = n ==> new_n  = L_cor^-1 * n
    Linv = la.inv(L_cor)
    new_n_list = Linv.dot(n_list.T).T.tolist()

    return new_n_list[0] if len(new_n_list) == 1 else new_n_list

def plane_trans(L_cor, n_list):
    """
    transform a list of planes in initial phase
    to final phase based on lattice correspondence

    :param L_cor: lattice correspondence
    :param n_list: list of planes
    """
    n_list = np.array(n_list)
    if len(n_list.shape) < 2:
        n_list = n_list.reshape(1, len(n_list))
    # new_n  = L_cor.T * n
    LT = L_cor.T
    new_n_list = LT.dot(n_list.T).T.tolist()
    return new_n_list[0] if len(new_n_list) == 1 else new_n_list
    



            
        
    
    
        
    
