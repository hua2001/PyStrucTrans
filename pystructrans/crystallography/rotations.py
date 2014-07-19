from ..general_imports import *
from math import radians

def Euler(angles, order="ZXZ", unit="deg"):
    """
    the rotation matrix of applying `Euler angles`_ in the given order,
    using the convention of `intrinsic rotations`_

    :param angles: three Euler angles, :class:`list`, :class:`tuple`, or :class:`numpy.ndarray`
    :keyword order: axis sequence, default is "ZXZ", :class:`str` of length 3
    :keyword unit: **"deg"** or **"rad"**, default is **"deg"**

    .. _Euler angles: http://en.wikipedia.org/wiki/Euler_angles
    .. _intrinsic rotations: http://en.wikipedia.org/wiki/Euler_angles#Intrinsic_rotations
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
    '''
    calculate the rotation matrix about axis z
    (in 3D only)

    :param angle: angle
    :param axis: axis as a 3D vector, :class:`list`, :class:`tuple` or :class:`numpy.ndarray`
    :keyword unit: "deg" or "rad", default is "deg"

    :return: 3 x 3 rotation matrix, :class:`numpy.ndarray`
    '''

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