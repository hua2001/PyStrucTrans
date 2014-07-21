from ..general_imports import *
from .. import rotation

def solvetwin(U, e, twintype, unsafe=False):
    """
    find the solutions to the twinning equation

    :param U: positive symmetric matrix, convertible to [3 x 3] :py:class:`numpy.ndarray`
    :param e: axis of a 180 degree rotation, convertible to [3 x 1] :py:class:`numpy.ndarray`
    :param str twintype: **"I"** (Type-I), **"II"** (Type-II) or **"C"** (compound)
    :keyword boolean unsafe: skip input check for a better performance, default is ``False``.
    :return: a tuple of two :py:class:`numpy.ndarray`, the first is `a` and the second is `n`
    :raises ValueError:

        * U cannot be converted to :py:class:`numpy.ndarray`
        * U is not 3 x 3
        * U is not positive definite
        * U is not symmetric

    """

    try:
        U = U if isinstance(U, np.ndarray) else np.array(U)
    except:
        raise ValueError("unrecognizable input for U")

    e = np.array(e)
    Q = rotation(180, e)

    if len(U) - 3 > 1e-9:
        print('The input should be a 3x3 symmetric matrix!')
        return

    if norm((U.T - U).reshape(9)) > 1e-9:
        print('The input should be a symmetric matrix!')
        return
    if abs(type - 1) < 1e-9:
        n = e
        denominator = np.dot(inv(U).dot(e), inv(U).dot(e))
#         print denominator
        a = 2*(np.dot(inv(U),e)/denominator - U.dot(e))
    else:
        if abs(type - 2) < 1e-9:
            n = 2*(e - np.dot(U, U.dot(e))/np.dot(U.dot(e), U.dot(e)))
            a = U.dot(e)
        else:
            print('Please input the type of twin system: 1 for Type I twin; 2 for Type II twin')
            return
    rho = norm(n)
    a = np.round(rho*a, 6)
    n = np.round(n/rho, 6)
    return np.array([a, n])
