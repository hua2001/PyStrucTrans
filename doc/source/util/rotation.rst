Rotation matrices
=================

Generic rotation
----------------
.. autofunction:: pystructrans.rotation

    :param angle: angle
    :param axis: axis as a 3D vector, :class:`list`, :class:`tuple` or :class:`numpy.ndarray`
    :keyword unit: "deg" or "rad", default is "deg"

    :return: 3 x 3 rotation matrix, :class:`numpy.ndarray`

Giving an ``angle`` (in the unit of |deg| by default)
and the rotation ``axis`` (not necessarily normalized),
we can get a matrix ``R``

.. doctest::

    >>> import numpy as np
    >>> import numpy.linalg as la
    >>> np.set_printoptions(suppress=True)
    >>> from pystructrans import rotation
    >>> angle = 120
    >>> axis = [1, 1, 1]
    >>> R = rotation(angle, axis)
    >>> print(R)
    [[ 0. -0.  1.]
     [ 1.  0. -0.]
     [-0.  1.  0.]]

We can check that it is in :math:`SO(3)`:

.. doctest::

  >>> print(la.det(R))
  1.0
  >>> print(R.dot(R.T) - np.eye(3))
  [[ 0.  0.  0.]
   [ 0.  0.  0.]
   [ 0.  0.  0.]]

We know that this ``R`` rotates the :math:`x` axis to :math:`y` axis, :math:`y` to :math:`z`,
and :math:`z` back to :math:`x`:

.. doctest::

    >>> x = [1, 0, 0]
    >>> y = R.dot(x)
    >>> print(y)
    [ 0.  1. -0.]

    >>> z = R.dot(y)
    >>> print(z)
    [-0.  0.  1.]

    >>> x2 = R.dot(z)
    >>> print(x2)
    [ 1. -0.  0.]

Euler angles
------------
.. autofunction:: pystructrans.Euler

    :param angles: three `Euler angles`_, :class:`list`, :class:`tuple`, or :class:`numpy.ndarray`
    :keyword order: sequence of `intrinsic rotations`_, default is "ZXZ", :class:`str` of length 3
    :keyword unit: **"deg"** or **"rad"**, default is **"deg"**

    .. _Euler angles: http://en.wikipedia.org/wiki/Euler_angles
    .. _intrinsic rotations: http://en.wikipedia.org/wiki/Euler_angles#Intrinsic_rotations

Define three ``angles`` and the ``order`` of three intrinsic rotation axis

.. doctest::

    >>> import numpy as np
    >>> import numpy.linalg as la
    >>> np.set_printoptions(suppress=True)

    >>> from pystructrans import Euler

    >>> angles = [90, 90, 90]
    >>> order = 'XYZ'
    >>> R = Euler(angles, order)
    >>> print(R)
    [[ 0.  0.  1.]
     [ 0.  1.  0.]
     [-1.  0.  0.]]

Again, it is in :math:`SO(3)`:

.. doctest::

  >>> print(la.det(R))
  1.0

  >>> print(R.dot(R.T) - np.eye(3))
  [[ 0.  0.  0.]
   [ 0.  0.  0.]
   [ 0.  0.  0.]]



.. |deg| unicode:: 0xB0 .. degree sign