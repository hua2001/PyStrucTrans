'''
This file contains functions and data related to lattice invariant shear matrices
'''
import numpy as np
import h5py

def get_shift(N, dim=3): 
    # generate the list of shift matrices for shifting from -N to N
    bd = np.arange(-N,N+1)
    M = [bd]
    for i in np.arange(1,dim**2):
        # tile the existing matrix
        M_tile = np.tile(M,(1, 2*N+1))
        # repeat the basic array
        M_rep = bd.repeat((2*N+1)**(i))
        # assemble the new M
        M = np.append(M_tile, [M_rep], axis = 0)    
        #print(M.shape)
    # delete small shifts (all components < N)
#     del_idx = np.array([])
#     for i in range(np.alen(M[0,:])):
#         if np.max(np.absolute(M[:,i])) < N:        
#             del_idx = np.append(del_idx, [i])
    del_idx = np.array([0])
    M = np.delete(M, del_idx, 1)
    #print(M.shape)
    
    f = h5py.File('shift_{:d}_dim_{:d}.hdf5'.format(int(N), int(dim)),'w')
    ds = f.create_dataset('shifts',M.shape,dtype='i')
    ds[...] = M
    f.close()
    #np.savetxt('shift_matrices_{:d}.txt'.format(int(N)), M, delimiter=',')
    
if __name__ == "__main__":
    get_shift(2, dim=2)    
