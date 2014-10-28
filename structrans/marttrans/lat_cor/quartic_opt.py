"""
This file is about finding the minimum value of a quartic form
"""
from structrans.general_imports import *
from scipy.optimize import minimize

def _quart_min(T):
    """
    :param T: a forth order tensor
    :return: the minimum quartic form of unit vectors
    """

    dim = int(np.sqrt(len(T)))
    # l0 = np.eye(dim).reshape(dim**2)
    l0 = - np.ones(dim**2) / dim
    constrains = (
        {'type': 'ineq', 'fun': lambda x: x.dot(x) - 1},
        {'type': 'ineq', 'fun': lambda x: 1 - x.dot(x)}
    )
    res = minimize(lambda x: T.dot(x).dot(x).dot(x).dot(x), l0, method="SLSQP", constraints=constrains)
    return res.fun




