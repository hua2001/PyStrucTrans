from setuptools import setup
import structrans

setup(name='PyStrucTrans',
      version=structrans.__verison__,
      description='A python package for structural phase transformation',
      url='http://github.com/structrans/pystructrans',
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
          'structrans',
          'structrans.crystallography',
          'structrans.marttrans',
          'structrans.util',
          'structrans.tests'
      ),
      install_requires=(
          'numpy >= 1.6.0'
      ),
      scripts=('bin/lattcorr',),
      entry_points={
        # 'console_scripts': []
      },
      data_files=[],
      test_suite='nose.collector',
      tests_require=['nose'],
      zip_safe=False)
