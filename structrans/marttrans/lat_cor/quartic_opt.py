"""
This file is about finding the minimum value of a quartic form
"""
from structrans.general_imports import *
try:
    from scipy.optimize import minimize

    def quart_min(T):
        """
        :param T: a forth order tensor
        :return: the minimum quartic form of unit vectors
        """

        dim = int(np.sqrt(len(T)))
        constrains = (
            {'type': 'eq', 'fun': lambda x: 1 - x.dot(x)}
        )
        def jac(x):
            D1 = T.transpose(0, 1, 2, 3).dot(x).dot(x).dot(x)
            D2 = T.transpose(1, 2, 3, 0).dot(x).dot(x).dot(x)
            D3 = T.transpose(2, 3, 0, 1).dot(x).dot(x).dot(x)
            D4 = T.transpose(3, 0, 1, 2).dot(x).dot(x).dot(x)
            return D1 + D2 + D3 + D4
        fmin = float('inf')
        for _ in xrange(100):
            l0 = np.random.rand(dim**2)
            l0 /= la.norm(l0)
            res = minimize(lambda x: T.dot(x).dot(x).dot(x).dot(x), l0, jac=jac, method="SLSQP", constraints=constrains)
            fmin = min(fmin, res.fun)
        return fmin

except ImportError:
    def _quartic_roots(T, x0, maxiter=1000, ftol=1e-4, rtol=1e-3):
        """

        :param T: fourth order tensor
        :param x0: initial guess
        :return: closest extremal
        """
        dim = len(x0)

        def func(xs):
            x = xs[:dim]
            lam = xs[-1]
            f = np.zeros(dim + 1)
            D1 = T.transpose(0, 1, 2, 3).dot(x).dot(x).dot(x)
            D2 = T.transpose(1, 2, 3, 0).dot(x).dot(x).dot(x)
            D3 = T.transpose(2, 3, 0, 1).dot(x).dot(x).dot(x)
            D4 = T.transpose(3, 0, 1, 2).dot(x).dot(x).dot(x)
            f[:dim] = D1 + D2 + D3 + D4 - 2 * lam * x
            f[-1] = 1 - x.dot(x)
            return f

        def jac(xs):
            x = xs[:dim]
            grad = np.zeros((dim + 1, dim + 1))

            D11 = T.transpose(0, 1, 2, 3).dot(x).dot(x)
            D12 = T.transpose(0, 2, 3, 1).dot(x).dot(x)
            D13 = T.transpose(0, 3, 1, 2).dot(x).dot(x)

            D21 = T.transpose(1, 2, 3, 0).dot(x).dot(x)
            D22 = T.transpose(1, 3, 0, 2).dot(x).dot(x)
            D23 = T.transpose(1, 0, 2, 3).dot(x).dot(x)

            D31 = T.transpose(2, 3, 0, 1).dot(x).dot(x)
            D32 = T.transpose(2, 0, 1, 3).dot(x).dot(x)
            D33 = T.transpose(2, 1, 3, 0).dot(x).dot(x)

            D41 = T.transpose(3, 0, 1, 2).dot(x).dot(x)
            D42 = T.transpose(3, 1, 2, 0).dot(x).dot(x)
            D43 = T.transpose(3, 2, 0, 1).dot(x).dot(x)

            grad[:dim, :dim] = (D11 + D12 + D13 + D21 + D22 + D23
                                + D31 + D32 + D33 + D41 + D42 + D43)
            grad[-1, :dim] = - 2 * x
            grad[:dim, -1] = - 2 * x

            return grad

        x = np.append(x0, [0.1], axis=0)
        f = func(x)
        df = float('inf')
        itr = 0
        while itr < maxiter and la.norm(f) > ftol and df > rtol and max(x) < 1e+20:
            g = jac(x)
            x1 = x - 0.2 * la.inv(g).dot(f)
            f1 = func(x1)
            df = la.norm((f1 - f))/la.norm(f)
            x, f = x1, f1
            itr += 1
        return x[:dim], x[-1]

    def quart_min(T):
        """
        :param T: a forth order tensor
        :return: the minimum quartic form of unit vectors
        """
        dim = int(np.sqrt(len(T)))
        fmin = float('inf')
        for _ in xrange(20):
            l0 = np.random.rand(dim**2)
            l0 /= la.norm(l0)
            x, mu = _quartic_roots(T, l0)
            f = T.dot(x).dot(x).dot(x).dot(x)
            fmin = min(fmin, f)
        return fmin