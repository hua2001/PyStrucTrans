# PyStrucTrans #

PyStrucTrans is a python package for crystallography and structural phase transformaiton analysis.

Download and upzip the repository. Install the package by running the following in the unzipped folder:

    python setup.py install

Then import the library in your python script

    import structrans

There is one shell command, `lattcorr`, that can be directly used after the installation. This command will invoke the python interpreter and search for the optimal lattice correspondences based on the provided lattice parameters. For example

    lattcorr -n 3 2 2 6 '1.414 2'

will search for the 3 best lattice correspondences from f.c.c. (id = 2) with *a*=2 to tetragonal (id = 6) with *a*=1.414 and *c*=2.

For more detail, refer to the documentation.

## Dependencies

Numpy >= 1.6.0

## Licence

[The MIT Lincense](http://opensource.org/licenses/MIT)
