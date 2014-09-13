from setuptools import setup
import pystructrans

setup(name='pystructrans',
      version=pystructrans.__verison__,
      description='A python package for structural phase transformation',
      url='http://github.com/yintaosong/pystructrans',
      author=('Yintao Song', 'Xian Chen'),
      author_email=('yintaosong@gmail.com', 'x.ch.msti@gmail.com'),
      license='LICENSE',
      classifiers = (
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Education',
          'Intended Audience :: Science/Research',
          'License :: Freeware',
          'Natural Language :: English',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: Implementation :: CPython'
      ),
      packages=(
          'pystructrans', 
          'pystructrans.crystallography',
          'pystructrans.marttrans',
          'pystructrans.util',
          'pystructrans.tests'
      ),
      install_requires=(
          'numpy >= 1.6.0'
      ),
      scripts=('bin/lattcorr',),
      entry_points={
        # 'console_scripts': [
        #                     'pst = pystructrans.command_line:main',
        #                     'lattcorr.py = pystructrans.command_line:run_lat_cor',
        #                     'main = pystructrans.command_line'
        #                     ],
      },
      data_files=[],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
