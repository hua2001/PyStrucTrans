Rotation matrices
===================

There are two 3D rotation matrix generating functions:
:py:func:`pystructrans.crystallography.rotation`
and
:py:func:`pystructrans.crystallography.Euler`
.

Generic rotation
----------------

.. autofunction:: pystructrans.crystallography.rotation

.. testsetup:: *

    from pystructrans.crystallography import rotation, Euler

Giving an ``angle`` (in the unit of |deg| by default)
and the rotation ``axis`` (not necessarily normalized),
we can get a matrix ``R``

.. testcode::

    import numpy as np
    import numpy.linalg as la
    np.set_printoptions(suppress=True)

    angle = 120
    axis = [1, 1, 1]
    R = rotation(angle, axis)
    print(R)

.. testoutput::

    [[ 0. -0.  1.]
     [ 1.  0. -0.]
     [-0.  1.  0.]]

We can check that it is in :math:`SO(3)`:

.. testcode::

  print(la.det(R))
  print(R.dot(R.T) - np.eye(3))

.. testoutput::

  1.0
  [[ 0.  0.  0.]
   [ 0.  0.  0.]
   [ 0.  0.  0.]]

We know that this ``R`` rotates the :math:`x` axis to :math:`y` axis, :math:`y` to :math:`z`,
and :math:`z` back to :math:`x`:

.. testcode::

    x = [1, 0, 0]
    y = R.dot(x)
    print(y)

    z = R.dot(y)
    print(z)

    x2 = R.dot(z)
    print(x2)

.. testoutput::

    [ 0.  1. -0.]
    [-0.  0.  1.]
    [ 1. -0.  0.]


Euler angles
------------
.. autofunction:: pystructrans.crystallography.Euler

Define three ``angles`` and the ``order`` of three intrinsic rotation axis

.. testcode::

    angles = [90, 90, 90]
    order = 'XYZ'
    R = Euler(angles, order)
    print(R)

.. testoutput::

    [[ 0.  0.  1.]
     [ 0.  1.  0.]
     [-1.  0.  0.]]

Again, it is in :math:`SO(3)`:

.. testcode::

  print(la.det(R))
  print(R.dot(R.T) - np.eye(3))

.. testoutput::

  1.0
  [[ 0.  0.  0.]
   [ 0.  0.  0.]
   [ 0.  0.  0.]]



.. |deg| unicode:: 0xB0 .. degree sign