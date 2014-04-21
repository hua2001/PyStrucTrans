from marttrans import lat_cor
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
#    parser.add_option("-f", "--file", dest="filename",
#                      help="write report to FILE", metavar="FILE")
#    parser.add_option("-q", "--quiet",
#                      action="store_false", dest="verbose", default=True,
#                      help="don't print status messages to stdout")
    # define options
    parser.add_option("-n","--number-of-solutions", dest="N", default=1, type="int", help="Number of solutions. default = 1")
    parser.add_option("-s","--save-results", dest="save", action='store_true', default=False, help="Save the results in an HDF5 file. default = False")
    parser.add_option("-f","--filename",dest="filename",metavar="FILE",default="lat_cor.hdf5",
                      help="The name of the HDF5 file that stores the results. default=lat_cor.hdf5")
    parser.add_option("-d","--display",dest="display",action='store_true', default=False,help="Display the calculation process and results. default = False")
    parser.add_option("-l","--logging",dest="logging",action='store_true', default=False,help="Save the log in \"lattcorr.log\". default = False")
    
    # parse the options
    (opts, args) = parser.parse_args()
    
    # read number of solutions
    num_sol = opts.N   
        
    # run the program
    lat_cor(int(args[0]), np.fromstring(args[1],sep=" "),
            int(args[2]), np.fromstring(args[3],sep=" "), 
            num_sol=num_sol, 
            disp=opts.display,
            save_results=opts.save,
            filename=opts.filename,
            save_log=opts.logging)
    
    
    
    
    