'''
Here is the description of the StructPhaseTrans package
'''
from __future__ import absolute_import

__verison__ = "1.0.0b"

def version():
    return __verison__

from .crystallography import Lattice, BravaisLattice
from .marttrans import *
from .util import Euler, rotation
