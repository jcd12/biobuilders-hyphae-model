# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 14:13:46 2015

@author: jana
"""

from distutils.core import setup
from Cython.Build import cythonize
import numpy
from distutils.extension import Extension

setup(
    name = 'Cythonized functions',
    ext_modules=cythonize([Extension("C_net_functions",["C_net_functions.pyx"],
    include_dirs=[numpy.get_include()])])
)    
