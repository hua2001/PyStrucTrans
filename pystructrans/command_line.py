from marttrans import lat_cor
from marttrans.lat_cor import logger as lc_logger
from marttrans.lat_opt import logger as lo_logger
import numpy as np
import sys
from optparse import OptionParser


def main():
    print('Welcome!')
    print sys.argv

def run_lat_cor():
    usage = '''
    usage: %prog [options] ibrava pbrava ibravm pbravm
        ibrava - Bravais lattice ID for austenite
        pbrava - lattice parameters for austenite
        ibravm - Bravais lattice ID for martensite
        pbravm - lattice parameters for martensite
        
    For example, the following command searching for 
    the best 2 lattice correspondences 
    from f.c.c. (id=2) austenite with a = 2 
    to simple tetragonal (id=6) martensite 
    with a = 1.414 and b = 2.
    
        lattcorr -n 2 2 2 6 '1.414 2'
        
    '''
    parser = OptionParser(usage=usage)
    # define options
    parser.add_option("-n","--num-sols", dest="nsol", default=1, type="int", help="number of solutions. default = 1")
    parser.add_option("-d","--display", dest="disp", default=2, help="verbosty of the calculation process and results. default = 2")
    parser.add_option("-l","--log-level", dest="loglevel", help="log level: DEBUG, INFO, WARNING, CRITICAL, NONE")
    parser.add_option("-f","--log-file", dest="logfile", help="log file")
    
    # parse the options
    (opts, args) = parser.parse_args()
    
    # run the program
    kwargs = {'nsol':opts.nsol, 'disp':opts.disp}
    
    if opts.loglevel:
        kwargs['loglevel'] = opts.loglevel
    if opts.logfile:
        kwargs['logfile'] = opts.logfile
        
    lat_cor(int(args[0]), np.fromstring(args[1],sep=" "),
            int(args[2]), np.fromstring(args[3],sep=" "),
            **kwargs)
    
    
    
    
    