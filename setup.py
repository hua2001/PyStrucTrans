from setuptools import setup
import pystructrans

setup(name='pystructrans',
      version=pystructrans.__verison__,
      description='A python package for structural phase transformation',
      url='http://github.com/storborg/funniest',
      author=['Xian Chen', 'Yintao Song'],
      author_email=['x.ch.msti@gmail.com', 'yintaosong@gmail.com'],
      license='GPL',
      classifiers = [
          'Development Status :: 2 - Pre-Alpha',
          'Environment :: Console',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: Freeware',
          'Natural Language :: English',
          'Programming Language :: Python :: 2.7'
      ],
      packages=[
          'pystructrans', 
          'pystructrans.crystallography',
          'pystructrans.marttrans',
          'pystructrans.tests'
      ],
      install_requires=[
          'numpy'
      ],
      entry_points = {
        'console_scripts': [
                            'pst = pystructrans.command_line:main',
                            'lattcorr = pystructrans.command_line:run_lat_cor'
                            ],
      },
      data_files=[('pystructrans/marttrans', ['pystructrans/marttrans/shift_1_dim_3.txt','pystructrans/marttrans/shift_1_dim_3.hdf5', 'pystructrans/marttrans/shift_1_dim_2.hdf5'])],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
