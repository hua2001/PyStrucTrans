#!/usr/bin/env python

import sys
from optparse import OptionParser
from pystructrans import lat_cor
from pystructrans.marttrans.lat_opt import lat_opt_unpack
from pystructrans.marttrans.lat_cor import reduce_hnfs
from multiprocessing import Pool
from numpy import fromstring
import logging

usage = '''
    usage: %prog [options] ibrava pbrava ibravm pbravm
        brava - Bravais lattice ID for austenite
        pa    - lattice parameters for austenite
        bravm - Bravais lattice ID for martensite
        pm    - lattice parameters for martensite

    For example, the following command searching for
    the best 2 lattice correspondences
    from f.c.c. (id=2) austenite with a = 2
    to simple tetragonal (id=6) martensite
    with a = 1.414 and b = 2.

        lattcorr -n 2 2 2 6 '1.414 2'

    '''
parser = OptionParser(usage=usage)
# define options
parser.add_option("-n", "--num-sols", dest="nsol", default=1, type="int", help="number of solutions. default = 1")
parser.add_option("-d", "--display", dest="disp", default=2, help="verbosty of the calculation process and results. default = 2")
parser.add_option("-l", "--log-level", dest="loglevel", default="NONE", help="log level: DEBUG, INFO, WARNING, CRITICAL, NONE")
parser.add_option("-f", "--log-file", dest="logfile", help="log file")

# parse the options
(opts, args) = parser.parse_args()

# run the program
kwargs = {'nsol': opts.nsol, 'disp': opts.disp}


loglevel_map = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'CRITICAL': logging.CRITICAL
}
if opts.loglevel is not "NONE":
    kwargs['loglevel'] = loglevel_map[opts.loglevel]
if opts.logfile:
    kwargs['logfile'] = opts.logfile

def reduce_hnfs_par(args):
    return pool.map(reduce_hnfs, args)

def lat_opt_par(args):
    return pool.apply_async(lat_opt_unpack, args)

kwargs['lat_opt_par'] = lat_opt_par
kwargs['reduce_hnfs_par'] = reduce_hnfs_par

if __name__ == '__main__':
    pool = Pool()
    lat_cor(int(args[0]), fromstring(args[1], sep=" "),
            int(args[2]), fromstring(args[3], sep=" "),
            **kwargs)