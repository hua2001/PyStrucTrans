import numpy as np
import numpy.linalg as la
from math import sqrt
from lattice import Lattice, inPointGroup

class BravaisLatticeError(Exception):
    pass

class BravaisLattice(Lattice):
    '''
    This is the class of Bravais lattices according to conventional classification. 
     
    Attributes:
     
    .. py:attribute:: __latID (integer) - ID of the bravais lattices.
    .. py:attribute:: __latParam (1D numpy array) - the list of lattice parameters.
    .. py:attribute:: __C (2D numpy array) - conversion between primitive and conventional unit cells
     
    Methods:
    '''  
 
    def __init__(self, latID, latParam, N=3):
        '''
        Constructor
        '''
        if isinstance(latParam, int) or isinstance(latParam, float):
            latParam = np.array([latParam])
        if len(np.where(latParam<=0)[0])>0:
            raise BravaisLatticeError('All lattice parameters musth be positive.')
        # 3D Bravais lattices        
        if N == 3:
            if not latID in range(1,15):
                raise BravaisLatticeError('no such Bravais lattice. ID must be between 1 and 14 for a 3D Bravais lattice.')
            pramlen = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 6] 
            if not np.alen(latParam) == pramlen[latID-1]:
                raise BravaisLatticeError('the number of lattice parameters does not match the type of Bravais lattice.')
            E, C = self.__primitiveBase3D__(latID, latParam)
        if N == 2:
            if not latID in range(1,6):
                raise BravaisLatticeError('no such Bravais lattice. ID must be between 1 and 5 for a 2D Bravais lattice.')
            pramlen = [1, 2, 2, 1, 3] 
            if not np.alen(latParam) == pramlen[latID-1]:
                raise BravaisLatticeError('the number of lattice parameters does not match the type of Bravais lattice.')
            E, C = self.__primitiveBase2D__(latID, latParam)
        Lattice.__init__(self, E)
        self.__latID = latID        
        self.__latParam = latParam
        self.__C = C
        
    def __primitiveBase2D__(self, n, p):
        '''
        generate the base matrix for primitive cell of 2D Bravais lattices
        '''
        C = np.eye(2)
        if n == 1: 
            # square
            e1 = p[0]*np.array([1,0])
            e2 = p[0]*np.array([0,1])
        elif n == 2: 
            # rectangular
            e1 = p[0]*np.array([1,0])
            e2 = p[1]*np.array([0,1])
        elif n == 3:
            # centered rectangular
            e1 = 0.5*np.array([p[0], p[1]])
            e2 = 0.5*np.array([-p[0], p[1]])
            C = 0.5*np.array([[1, -1],[1, 1]])
        elif n == 4:
            # hexagonal            
            e1 = p[0]*np.array([1, 0])
            e2 = p[0]*np.array([0.5, np.sqrt(3)])
        elif n == 5:
            # oblique
            gamma = p[2]*np.pi/180.
            e1 = p[0]*np.array([1,0])
            e2 = p[1]*np.array([np.cos(gamma), np.sin(gamma)])        
        return np.array([e1, e2]).T, la.inv(C)
    
    def __primitiveBase3D__(self, n, p):
        '''
        generate the base matrix for primitive cell of 3D Bravais lattices
        '''
        C = np.eye(3)
        if n == 1: 
            # simple cubic
            e1 = p[0]*np.array([1,0,0])
            e2 = p[0]*np.array([0,1,0])
            e3 = p[0]*np.array([0,0,1])
        elif n == 2: 
            # face centered cubic
            e1 = 0.5*p[0]*np.array([1,1,0])
            e2 = 0.5*p[0]*np.array([0,1,1])
            e3 = 0.5*p[0]*np.array([1,0,1])
            C = 0.5*np.array([[ 1,  0,  1],
                              [ 1,  1,  0],
                              [ 0,  1,  1]])        
        elif n==3:
            #bcc
            e1 = 0.5*p[0]*np.array([ 1, 1, 1])
            e2 = 0.5*p[0]*np.array([-1, 1, 1])
            e3 = 0.5*p[0]*np.array([-1,-1, 1])
            C = 0.5*np.array([[ 1, -1, -1],
                              [ 1,  1, -1],
                              [ 1,  1,  1]])
        elif n==4:
            #hexagonal
            e1 = p[0]*np.array([1,0,0])
            e2 = p[0]*np.array([0.5,np.sqrt(3)/2, 0])
            e3 = np.array([0,0,p[1]])
        elif n==5:
            #trigonal
            #<111> is the 3fold axis
            c = np.cos(p[1]*np.pi/180)
            a = p[0]
            ty = np.sqrt((1-c)/6)
            tz = np.sqrt((1+2*c)/3)
            u = tz - 2*np.sqrt(2)*ty; 
            v = tz + np.sqrt(2)*ty;
            e1 = a/np.sqrt(3)*np.array([u,v,v])
            e2 = a/np.sqrt(3)*np.array([v,u,v])
            e3 = a/np.sqrt(3)*np.array([v,v,u])
        elif n==6:
            #simple tetragonal
            a = p[0]; c = p[1];
            e1 = a*np.array([1,0,0])
            e2 = a*np.array([0,1,0])
            e3 = c*np.array([0,0,1])
        elif n==7:
            #body centered tetragonal
            a = p[0]; c = p[1];
            e1 = (a/2)*np.array([ 1, 1, c/a])
            e2 = (a/2)*np.array([-1, 1, c/a])
            e3 = (a/2)*np.array([-1,-1, c/a])
            C = 0.5*np.array([[ 1, -1, -1],
                              [ 1,  1, -1],
                              [ 1,  1,  1]])
        elif n==8:
            #simple orthorhombic
            a = p[0]; b = p[1]; c = p[2];
            e1 = np.array([a,0,0])
            e2 = np.array([0,b,0])
            e3 = np.array([0,0,c])
        elif n==9:
            #base centered orthorhombic
            a = p[0]; b = p[1]; c = p[2];
            e1 = np.array([a/2, b/2,0])
            e2 = np.array([-a/2,b/2,0])
            e3 = np.array([0,0,c])
            C = np.array([[0.5, -0.5, 0],
                          [0.5, 0.5, 0],
                          [  0,   0, 1]])
        elif n==10:
            #face centered orthorhombic
            a = p[0]; b = p[1]; c = p[2];
            e1 = np.array([a/2, b/2,   0])
            e2 = np.array([  0, b/2, c/2])
            e3 = np.array([a/2,   0, c/2])
            C = 0.5*np.array([[ 1,  0,  1],
                              [ 1,  1,  0],
                              [ 0,  1,  1]])  
        elif n==11:
            #body centered orthorhombic
            a = p[0]; b = p[1]; c = p[2];
            e1 = np.array([a/2,b/2,c/2])
            e2 = np.array([-a/2,b/2,c/2])
            e3 = np.array([-a/2,-b/2,c/2])
            C = 0.5*np.array([[ 1, -1, -1],
                              [ 1,  1, -1],
                              [ 1,  1,  1]])
        elif n==12:
            #monoclinic unique axis b
            a = p[0]; b = p[1]; c = p[2]; beta = p[3]*np.pi/180.;
            e1 = np.array([a,0,0])
            e2 = np.array([0,b,0])
            e3 = np.array([c*np.cos(beta),0,c*np.sin(beta)])
        elif n== 13:
            #base centered monoclinic
            a = p[0]; b = p[1]; c = p[2]; beta = p[3]*np.pi/180;
            e1 = np.array([a/2, b/2, 0])
            e2 = np.array([-a/2, b/2, 0])
            e3 = np.array([c*np.cos(beta),0,c*np.sin(beta)])
            C = np.array([[0.5, -0.5, 0],
                          [0.5, 0.5, 0],
                          [  0,    0, 1]])
        elif n==14:
            #triclinic
            a = p[0]; b = p[1]; c = p[2]; alpha = p[3]*np.pi/180; beta = p[4]*np.pi/180; gamma = p[5]*np.pi/180;
            e1 = np.array([a,0,0])
            e2 = np.array([b*np.cos(gamma), b*np.sin(gamma),0])
            e3 = np.array([c*np.cos(beta),c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma),
                        c*np.sqrt(1+2*np.cos(alpha)*np.cos(beta)*np.cos(gamma) - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2)/np.sin(gamma)])
    
        return np.array([e1, e2, e3]).T, la.inv(C)
    
    def getLatID(self):
        '''
        get the lattice ID
        '''
        return self.__latID
    
    def getLatParam(self):
        '''
        get the lattice parameters
        '''
        return self.__latParam
    
    def __str__(self):
        '''
        Get the description of the Bravais lattice
        
        :return: (*string*) lattice description
        '''
        if self.getDimension() == 3:
            des_str = '3D Bravais Lattice - '
            des_list = ['Cubic P', 
                        'Cubic F (face centered cubic)',
                        'Cubic I (body centered cubic)',
                        'Hexagonal P',
                        'Rhombohedral P',
                        'Tetragonal P',
                        'Tetragonal I (body centered tetragonal)',
                        'Orthorhombic P',
                        'Orthorhombic C (base centered orthorhombic)',
                        'Orthorhombic F (face centered orthorhombic)',
                        'Orthorhombic I (body centered orthorhombic)',
                        'Monoclinic P',
                        'Monoclinic C (base centered monoclinic)',
                        'Triclinic P']
            # print the type of Bravais lattice
            des_str = des_str + des_list[self.__latID-1]
            # print the first lattice parameter
            des_str = des_str+':    a = {:g}'.format(self.__latParam[0])
            # print b if it exists
            if self.__latID in range(8,15):
                des_str = des_str+',    b = {:g}'.format(self.__latParam[1])
            #print c if it exists
            if self.__latID in np.array([4,6,7]):
                des_str = des_str+',    c = {:g}'.format(self.__latParam[1])
            if self.__latID in range(8,15):
                des_str = des_str+',    c = {:g}'.format(self.__latParam[2])
            # print alpha if it exists
            if self.__latID == 5.:
                des_str = des_str+u',    alpha = {:g}'.format(self.__latParam[1])
            if self.__latID == 14.:
                des_str = des_str+u',    alpha = {:g}'.format(self.__latParam[3])
            # print beta if it exists
            if self.__latID in np.array([12., 13.]):
                des_str = des_str+u',    beta = {:g}'.format(self.__latParam[3])
            if self.__latID == 14.:
                des_str = des_str+u',    beta = {:g}'.format(self.__latParam[4])
            # print gamma if it exists
            if self.__latID == 14.:
                des_str = des_str+u',    gamma = {:g}'.format(self.__latParam[5])
        elif self.getDimension() == 2:
            des_str = '2D Bravais Lattice - '
            des_list = ['Square P', 
                        'Rectangular P',
                        'Rectangular I (centered rectangular)',
                        'Hexagonal P',
                        'Oblique P']
            # print the type of Bravais lattice
            des_str = des_str + des_list[self.__latID-1]
            # print the first lattice parameter
            des_str = des_str+':    a = {:g}'.format(self.__latParam[0])
            # print b if it exists
            if self.__latID in np.array([2, 3, 5]):
                des_str = des_str+',    b = {:g}'.format(self.__latParam[1])
            # print gamma if it exists
            if self.__latID == 5.:
                des_str = des_str+u',    gamma = {:g}'.format(self.__latParam[2])            
        else:
            raise BravaisLatticeError('only 2D and 3D Bravais lattices are acceptable.')  
              
        return des_str
    
    def getConventionalTrans(self):
        '''
        get the conversion matrix from the primitive base to the conventional base
        
        prim_base = conv_base * trans_mat
        
        or 
        
        conv_idx = trans_mat * prim_idx
        
        :return: (*2D numpy array*) conversion matrix
        '''
        return self.__C
    
    def getConventionalBase(self):
        '''
        get the base matrix for conventional unit cell
        
        :return: (*2D numpy array*) conventional base matrix
        '''
        return np.dot(self.getBase(), self.getConventionalTrans())

    
    