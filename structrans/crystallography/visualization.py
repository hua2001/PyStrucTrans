from structrans.general_imports import *
from .bravais_lattice import BravaisLattice
from numpy.linalg import inv, norm

class visualError(Exception):
    pass

def vertex(o, E, n=3):
    dim = np.array(E).shape
    if n==3:
        if dim[0]==3 and dim[1]==3 and len(o)==3:
            return np.array([o,
                              o+E[0],o+E[1],o+E[2],
                              o+E[0]+E[1],o+E[1]+E[2],o+E[0]+E[2],
                              o+E[0]+E[1]+E[2]])
        else:
            print('E has to be a 3x3 matrix.')
            return None
    elif n==2:
        if dim[0]==2 and dim[1]==2 and len(o)==2:
            return np.array([o,
                             o+E[0], o+E[1],
                             o+E[0]+E[1]])
        else:
            print('E has to be a 2x2 matrix.')
            return None
    else:
        raise visualError('The dimension n has to be 2 or 3!')
    
def bonds(verts):
        '''
        Verts are a list of vertices, 
        which are 8x3 for 3D, 4x2 for 2D
        '''
        verts = np.array(verts)
                
        if len(verts[0])==3:
            idx_1 = [0,0,0,1,1,4,4,7,7,6,3,5]
            idx_2 = [1,2,3,4,6,2,7,6,5,3,5,2]
            idx = np.vstack((idx_1, idx_2)).T
            bds = [[verts[i[0]], verts[i[1]]] for i in idx]
            return np.array(bds)

        elif len(verts[0])==2:
            idx = [[0,1],[0,2],[1,3],[2,3]]
            return np.array([[verts[i[0]], verts[i[1]]] for i in len(idx)])
        else:
            raise visualError('The number of vertices is not valid.')
    
def box(v_list, E, n=3):
    dim = np.array(E).shape
    v_list = np.array(v_list)
    E = np.array(E)
    inside = []
    if n==3:
        if dim[0]==3 and dim[1]==3 and len(v_list[0])==3:
#             cr = norm(E[:,0])+norm(E[:,1])+norm(E[:,2])
            
            for v in v_list:
                v1 = v.dot(E[:,0]/norm(E[:,0]))
                v2 = v.dot(E[:,1]/norm(E[:,1]))
                v3 = v.dot(E[:,2]/norm(E[:,2]))
                if v1<1.05*norm(E[:,0]) and v2<1.05*norm(E[:,1]) and v3<1.05*norm(E[:,2]):
                    if v1>=0 and v2>=0 and v3>=0:
                        inside.append(v)
        else:
            print("E has to be a 3x3 matrix!")
            return None
    elif n==2:
        if dim[0]==2 and dim[1]==2 and len(v_list[0])==2:
            cr = norm(E[:,0])+norm(E[:,1])
            for v in v_list:
                if norm(v.dot(E))<=1.05*cr:
                    inside.append(v)
        else:
            print("E has to be a 2x2 matrix!")
    else:
        raise visualError("The dimension n has to be 2 or 3!")
    if np.array(inside).shape[0]>0:
        return np.array(inside)
    else:
        return np.array([0, 0, 0])
    

# def bonds(o, E, n=3):
#     verts = vertex(o, E, n)
#     if len(verts)==8:
#         idx_1 = [0,0,0,1,1,4,4,7,7,6,3,5]
#         idx_2 = [1,2,3,4,6,2,7,6,5,3,5,2]
# #         idx = np.array([0,1,0,2,0,3,1,4,1,5,2,4,3,5,3,6,4,7,5,7,7,6,2,6])
#         idx = np.vstack((idx_1, idx_2)).T
#         bds = [[verts[i[0]], verts[i[1]]] for i in idx]
#         return np.array(bds)
#     elif len(verts)==4:
#         idx = [[0,1],[0,2],[1,3],[2,3]]
#         return np.array([[verts[i[0]], verts[i[1]]] for i in len(idx)])
#     else:
#         raise visualError('The number of vertices is not valid.')
        
            
class UnitCell():
    r'''
    unit_cell class gives the vertices
    it can be constructed by a lattice (2D or 3D).
    
    '''
    def __init__(self, origin, lat_ID, lat_Params, N=3):
        r'''
        Constructed by a Bravais lattice
        '''
        self._lattice = BravaisLattice(lat_ID, lat_Params, N=3)
        self._o = np.array(origin)
        self._dim = N
    def setOrigin(self, new_origin):
        self._o = new_origin
        
    def bonds(self, verts):
        '''
        Verts are a list of vertices, 
        which are 8x3 for 3D, 4x2 for 2D
        '''
        verts = np.array(verts)
                
        if len(verts[0])==3:
            idx_1 = [0,0,0,1,1,4,4,7,7,6,3,5]
            idx_2 = [1,2,3,4,6,2,7,6,5,3,5,2]
            idx = np.vstack((idx_1, idx_2)).T
            bds = [[verts[i[0]], verts[i[1]]] for i in idx]
            return np.array(bds)

        elif len(verts[0])==2:
            idx = [[0,1],[0,2],[1,3],[2,3]]
            return np.array([[verts[i[0]], verts[i[1]]] for i in len(idx)])
        else:
            raise visualError('The number of vertices is not valid.')
        
    def getPrimitive(self):
        E = self._lattice.getbase().T
        return vertex(self._o, E, self._dim)
    def getConventionalVex(self):
        C = self._lattice.getConventionalBase().T
        return vertex(self._o, C, self._dim)
    def getConventional(self):
        C = self._lattice.getConventionalBase().T
        inside = box(self.getPrimitive(), C, self._dim)
#         print np.array(inside).shape
        return np.vstack((vertex(self._o, C, self._dim), inside))
    def getConventionalbonds(self):
        verts = vertex(self._o, self._lattice.getConventionalBase().T, self._dim)
        return self.bonds(verts)
    def getPrimitivebonds(self):
        verts = vertex(self._o, self._lattice.getBase().T, self._dim)
        return self.bonds(verts)
    def copyTo(self, new_origin, v_list):
#         verts = self.getConventional()
        return v_list+np.vstack([new_origin]*len(v_list))

            
        
        
        
    
    