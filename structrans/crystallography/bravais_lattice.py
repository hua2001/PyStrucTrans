from structrans.general_imports import *
from .lattice import Lattice
from math import radians

class BravaisLattice(Lattice):
    '''
    This is the class of Bravais lattices according to conventional classification.

    :param latID: id of the Bravais lattice type
    :param latParam: lattice parameters
    :param dim: dimension of the Bravais lattice, 2 or 3.
    '''  
 
    def __init__(self, latID, latParam, dim=3):
        '''
        Constructor
        '''
        try:
            if not isinstance(latParam, np.ndarray):
                if isinstance(latParam, list):
                    latParam = np.array(latParam)
                else:
                    latParam = np.array([latParam])
            if len(latParam.shape) > 1:
                raise ValueError('lattice parameters is not a 1D list nor single value')
            elif any(latParam <= 0):
                raise ValueError('non-positive lattice parameters')
        except Exception:
            raise ValueError('invalid lattice parameters')

        # 3D Bravais lattices
        if dim == 3:
            if not latID in range(1, 15):
                raise ValueError('3D Bravais lattice ID must between 1 and 14')
            pramlen = [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 6] 
            if not np.alen(latParam) == pramlen[latID - 1]:
                raise ValueError('the number of lattice parameters does not match the type of Bravais lattice.')
            E, C = BravaisLattice.__primitiveBase3D__(latID, latParam)
        if dim == 2:
            if not latID in range(1, 6):
                raise ValueError('2D Bravais lattice ID must be between 1 and 5.')
            pramlen = [1, 2, 2, 1, 3] 
            if not np.alen(latParam) == pramlen[latID-1]:
                raise ValueError('the number of lattice parameters does not match the type of Bravais lattice.')
            E, C = BravaisLattice.__primitiveBase2D__(latID, latParam)

        Lattice.__init__(self, E)
        self.__latID = latID        
        self.__latParam = latParam
        self.__C = C

    @staticmethod
    def __primitiveBase2D__(n, p):
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
            C = 0.5*np.array([[1, -1], [1, 1]])
        elif n == 4:
            # hexagonal            
            e1 = p[0] * np.array([1, 0])
            e2 = p[0] * 0.5 *np.array([1, np.sqrt(3)])
        elif n == 5:
            # oblique
            gamma = p[2]*np.pi/180.
            e1 = p[0]*np.array([1,0])
            e2 = p[1]*np.array([np.cos(gamma), np.sin(gamma)])        
        return np.array([e1, e2]).T, la.inv(C)

    @staticmethod
    def __primitiveBase3D__(n, p):
        '''
        generate the base matrix for primitive cell of 3D Bravais lattices
        '''
        C = np.eye(3)
        if n == 1:
            # simple cubic
            e1 = p[0] * np.array([1, 0, 0])
            e2 = p[0] * np.array([0, 1, 0])
            e3 = p[0] * np.array([0, 0, 1])
        elif n == 2:
            # face centered cubic
            e1 = 0.5 * p[0] * np.array([1, 1, 0])
            e2 = 0.5 * p[0] * np.array([0, 1, 1])
            e3 = 0.5 * p[0] * np.array([1, 0, 1])
            C = 0.5 * np.array([[1, 0, 1],
                                [1, 1, 0],
                                [0, 1, 1]])
        elif n == 3:
            # bcc
            e1 = 0.5 * p[0] * np.array([1, 1, 1])
            e2 = 0.5 * p[0] * np.array([-1, 1, 1])
            e3 = 0.5 * p[0] * np.array([-1, -1, 1])
            C = 0.5 * np.array([[1, -1, -1],
                                [1, 1, -1],
                                [1, 1, 1]])
        elif n == 4:
            # hexagonal
            e1 = p[0] * np.array([1, 0, 0])
            e2 = p[0] * np.array([0.5, np.sqrt(3) / 2, 0])
            e3 = np.array([0, 0, p[1]])
        elif n == 5:
            # trigonal
            #<111> is the 3fold axis
            c = np.cos(p[1] * np.pi / 180)
            a = p[0]
            ty = np.sqrt((1 - c) / 6)
            tz = np.sqrt((1 + 2 * c) / 3)
            u = tz - 2 * np.sqrt(2) * ty
            v = tz + np.sqrt(2) * ty
            e1 = a / np.sqrt(3) * np.array([u, v, v])
            e2 = a / np.sqrt(3) * np.array([v, u, v])
            e3 = a / np.sqrt(3) * np.array([v, v, u])
        elif n == 6:
            # simple tetragonal
            a = p[0]
            c = p[1]
            e1 = a * np.array([1, 0, 0])
            e2 = a * np.array([0, 1, 0])
            e3 = c * np.array([0, 0, 1])
        elif n == 7:
            # body centered tetragonal
            a = p[0]
            c = p[1]
            e1 = (a / 2) * np.array([1, 1, c / a])
            e2 = (a / 2) * np.array([-1, 1, c / a])
            e3 = (a / 2) * np.array([-1, -1, c / a])
            C = 0.5 * np.array([[1, -1, -1],
                                [1, 1, -1],
                                [1, 1, 1]])
        elif n == 8:
            # simple orthorhombic
            a = p[0]
            b = p[1]
            c = p[2]
            e1 = np.array([a, 0, 0])
            e2 = np.array([0, b, 0])
            e3 = np.array([0, 0, c])
        elif n == 9:
            # base centered orthorhombic
            a = p[0]
            b = p[1]
            c = p[2]
            e1 = np.array([a / 2, b / 2, 0])
            e2 = np.array([-a / 2, b / 2, 0])
            e3 = np.array([0, 0, c])
            C = np.array([[0.5, -0.5, 0],
                          [0.5, 0.5, 0],
                          [0, 0, 1]])
        elif n == 10:
            # face centered orthorhombic
            a = p[0]
            b = p[1]
            c = p[2]
            e1 = np.array([a / 2, b / 2, 0])
            e2 = np.array([0, b / 2, c / 2])
            e3 = np.array([a / 2, 0, c / 2])
            C = 0.5 * np.array([[1, 0, 1],
                                [1, 1, 0],
                                [0, 1, 1]])
        elif n == 11:
            # body centered orthorhombic
            a = p[0]
            b = p[1]
            c = p[2]
            e1 = np.array([a / 2, b / 2, c / 2])
            e2 = np.array([-a / 2, b / 2, c / 2])
            e3 = np.array([-a / 2, -b / 2, c / 2])
            C = 0.5 * np.array([[1, -1, -1],
                                [1, 1, -1],
                                [1, 1, 1]])
        elif n == 12:
            # monoclinic unique axis b
            a = p[0]
            b = p[1]
            c = p[2]
            beta = radians(p[3])
            e1 = np.array([a, 0, 0])
            e2 = np.array([0, b, 0])
            e3 = np.array([c * np.cos(beta), 0, c * np.sin(beta)])
        elif n == 13:
            # base centered monoclinic
            a = p[0]
            b = p[1]
            c = p[2]
            beta = radians(p[3])
            e1 = np.array([a / 2, b / 2, 0])
            e2 = np.array([-a / 2, b / 2, 0])
            e3 = np.array([c * np.cos(beta), 0, c * np.sin(beta)])
            C = np.array([[0.5, -0.5, 0],
                          [0.5, 0.5, 0],
                          [0, 0, 1]])
        elif n == 14:
            # triclinic
            a = p[0]
            b = p[1]
            c = p[2]
            alpha = radians(p[3])
            beta = radians(p[4])
            gamma = radians(p[5])
            e1 = np.array([a, 0, 0])
            e2 = np.array([b * np.cos(gamma), b * np.sin(gamma), 0])
            e3 = np.array([c * np.cos(beta), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                           c * np.sqrt(1 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.cos(alpha) ** 2
                                       - np.cos(beta) ** 2 - np.cos(gamma) ** 2) / np.sin(gamma)])

        return np.array([e1, e2, e3]).T, la.inv(C)

    def __str__(self):
        '''
        Get the description of the Bravais lattice
        
        :return: (*string*) lattice description
        '''
        if self.dimension() == 3:
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
            des_str += des_list[self.__latID - 1]
            # print the first lattice parameter
            des_str = des_str+':    a = {:g}'.format(self.__latParam[0])
            # print b if it exists
            if self.__latID in range(8, 15):
                des_str = des_str+',    b = {:g}'.format(self.__latParam[1])
            #print c if it exists
            if self.__latID in np.array([4,6,7]):
                des_str = des_str+',    c = {:g}'.format(self.__latParam[1])
            if self.__latID in range(8, 15):
                des_str = des_str+',    c = {:g}'.format(self.__latParam[2])
            # print alpha if it exists
            if self.__latID in [5, 14]:
                des_str = des_str+u',    alpha = {:g}'.format(self.__latParam[1])
            # print beta if it exists
            if self.__latID in np.array([12., 13.]):
                des_str = des_str+u',    beta = {:g}'.format(self.__latParam[3])
            if self.__latID == 14:
                des_str = des_str+u',    beta = {:g}'.format(self.__latParam[4])
                des_str = des_str+u',    gamma = {:g}'.format(self.__latParam[5])
        else:
            des_str = '2D Bravais Lattice - '
            des_list = ['Square P', 
                        'Rectangular P',
                        'Rectangular I (centered rectangular)',
                        'Hexagonal P',
                        'Oblique P']
            # print the type of Bravais lattice
            des_str += des_list[self.__latID - 1]
            # print the first lattice parameter
            des_str = des_str+':    a = {:g}'.format(self.__latParam[0])
            # print b if it exists
            if self.__latID in np.array([2, 3, 5]):
                des_str = des_str+',    b = {:g}'.format(self.__latParam[1])
            # print gamma if it exists
            if self.__latID == 5.:
                des_str = des_str+u',    gamma = {:g}'.format(self.__latParam[2])
              
        return des_str

    def Bravais_id(self):
        '''
        get the Bravaise lattice ID
        '''
        return self.__latID

    def getLatID(self):
        warn("getLatID() is deprecated, please use Bravais_id() instead.", DeprecationWarning)
        return self.Bravais_id()

    def lat_param(self):
        '''
        get the lattice parameters
        '''
        return self.__latParam

    def getLatParam(self):
        warn("getLatParam() is deprecated, please use lat_param() instead.", DeprecationWarning)
        return self.lat_param()

    def conventional_trans(self):
        '''
        get the conversion matrix from the primitive base to the conventional base

        - conv_base = prim_base * **cnvr_mat**
        - prim_idx = **cnvr_mat** * conv_idx
        '''
        return self.__C

    def getConventionalTrans(self):
        warn("getConventionalTrans() is deprecated, please use conventional_trans() instead.", DeprecationWarning)
        return self.conventional_trans()

    def conventional_base(self):
        '''
        get the base matrix for conventional unit cell
        '''
        return np.dot(self.base(), self.conventional_trans())

    def getConventionalBase(self):
        warn("getConventionalBase() is deprecated, please use conventional_base() instead.", DeprecationWarning)
        return self.conventional_base()

    def to_primitive_Miller(self, idx):
        """
        convert a Miller index or a list of Miller indices
        from conventional base to primitive base

        :param idx: a Miller index or a list of Miller indices
        :return: the Miller index or the list of Miller indices in primitive base
        """
        try:
            idx = np.array(idx)
            if idx.ndim == 1:
                idx = np.array([idx])
            elif idx.ndim > 2:
                raise ValueError("invalid input")
            if idx.shape[1] is not self.dimension():
                raise ValueError("dimensions of lattice and indices not match")
        except:
            raise ValueError("invalid input of index (indices)")

        res = (self.__C.dot(idx.T)).T
        return res[0] if len(idx) == 1 else res

    def toPrimitive(self, idx):
        warn("toPrimitive() is deprecated, please use to_primitive() instead.", DeprecationWarning)
        return self.to_primitive()

    def to_conventional_Miller(self, idx):
        """
        convert a Miller index or a list of Miller indices
        from primitive to conventional

        :param idx: a Miller index or a list of Miller indices
        :return: the Miller index or the Miller list of indices in conventional base
        """
        try:
            if not isinstance(idx, np.ndarray):
                idx = np.array(idx)
            if idx.ndim == 1:
                idx = np.array([idx])
            elif idx.ndim > 2:
                raise ValueError("invalid input")
            if idx.shape[1] is not self.dimension():
                raise ValueError("dimensions of lattice and indices not match")
        except:
            raise ValueError("invalid input of index (indices)")

        res = (la.inv(self.__C).dot(idx.T)).T
        return res[0] if len(idx) == 1 else res

    def toConventional(self, idx):
        warn("toConventional() is deprecated, please use to_conventional_Miller() instead.", DeprecationWarning)
        return self.to_conventional_Miller()