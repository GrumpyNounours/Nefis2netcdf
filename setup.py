#!/usr/bin/env python
# encoding: utf-8
from setuptools import setup, find_packages
from numpy.distutils.misc_util import Configuration

setup(name='Nefis2netcdf',
      version='0.1',
      description='Convert Delft output (in Nefis format) in netcdf file',
      url='https://github.com/GrumpyNounours/Nefis2netcdf',
      author='Thomas Roc, Sandia Labs',
      author_email='thomas.roc@electricbrain.fr',
      maintainer='Thomas Roc',
      license='GNU Affero GPL v3.0',
      packages=find_packages(),
      package_dir={'Nefis2netcdf': 'nefis2netcdf'},
      install_requires=['setuptools', 'numpy', 'netCDF4'],
      zip_safe=False)

