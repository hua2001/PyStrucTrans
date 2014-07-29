from ..general_imports import *
from math import radians

def is_3_by_3(U):
    """
    U is a [3 x 3] matrix
    """
    return isinstance(U, np.ndarray) and U.shape is not (3, 3)

def pos_def_sym(U):
    """
    check if U is [3 x 3] positive definite and symmetric
    """
    return is_3_by_3(U) and la.det(U) > 0 and np.array_equal(U, U.T)

def sort_eig(A):
    """
    sorted eigenvalues and eigenvectors
    :param A:
    :return:
    """
    e = la.eig(A)
    e = [np.append(e[0][i], e[1][:, i]).real for i in xrange(3)]
    e = np.array(sorted(e, key=lambda x: x[0]))
    eval = e[:, 0]
    evec = np.array([[1., 0, 0], [0, 1., 0], [0, 0, 1.]]).dot(e[:, 1:4])
    return eval, evec

def divisors(n):
    """
    Find all divisors of an integer n

    :param n: the integer to be factorized
    :type n: integer
    :return: vector - the tuple of all the divisors of n
    """

    dvrs = [i + 1
            for i in xrange(int(np.floor(n / 2)))
            if n % (i + 1) == 0]
    dvrs.append(n)
    return tuple(dvrs)

def Euler(angles, order="ZXZ", unit="deg"):
    """
    the rotation matrix of applying Euler angles in the given order.
    """
    if len(angles) != 3 or len(order) != 3:
        raise ValueError("number of rotations is not three")
    canonical_rotation = {
        'X': lambda x: np.array([[1, 0, 0], [0, np.cos(x), -np.sin(x)], [0, np.sin(x), np.cos(x)]]),
        'Y': lambda x: np.array([[np.cos(x), 0, np.sin(x)], [0, 1, 0], [-np.sin(x), 0, np.cos(x)]]),
        'Z': lambda x: np.array([[np.cos(x), -np.sin(x), 0], [np.sin(x), np.cos(x), 0], [0, 0, 1]])
    }
    R = np.eye(3)
    for i, a in enumerate(order):
        if a not in "XYZxyz":
            raise ValueError("unknown rotation axis {:s}", a)
        else:
            angle = radians(angles[i]) if unit is "deg" else angles[i]
            R = canonical_rotation[a.upper()](angle).dot(R)
    return R


def rotation(angle, axis, unit="deg"):
    """
    calculate the rotation matrix about axis z (in 3D only)
    """

    a = axis if isinstance(axis, np.ndarray) else np.array(axis)
    z = a/la.norm(a)
    t = radians(angle) if unit is "deg" else angle
    if z.dot(np.array([1, 0, 0])) < 0.9:
        # z is not parallel to e1
        z1 = np.cross(z, np.array([1, 0, 0]))
        z2 = np.cross(z, z1)
    else:
        # z is parallel to e1
        z1 = np.cross(z, np.array([0, 1, 0]))
        z2 = np.cross(z, z1)
    z1 /= la.norm(z1)
    z2 /= la.norm(z2)
    R1 = np.cos(t)*(np.outer(z1, z1) + np.outer(z2, z2))
    R2 = np.sin(t)*(np.outer(z2, z1) - np.outer(z1, z2))
    return R1 + R2 + np.outer(z, z)
