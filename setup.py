#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
#
#   PyChemkin
#
################################################################################

import numpy
import sys

if __name__ == '__main__':
    
    # Use setuptools by default (requires Python>=2.6) if available
    # If not available, fall back to distutils
    # Using setuptools enables support for compiling wheels
    try:
        from setuptools import setup, Extension
    except ImportError:
        from distutils.core import setup
        from distutils.extension import Extension

    modules = ['pychemkin.chemkin', 
               'pychemkin.ignition'
               'pychemkin.sensitivity'
               ]

    # Read the version number
    exec(open('pychemkin/version.py').read())
    # Run the setup command
    setup(name='PyChemkin',
        version=__version__,
        description='A Python interface for CHEMKIN analysis',
        author='Connie W. Gao and the Reaction Mechanism Generator Team',
        author_email='rmg_dev@mit.edu',
        url='http://github.com/GreenGroup/PyChemkin',
        py_modules= modules,
        packages = ['pychemkin']
    )
