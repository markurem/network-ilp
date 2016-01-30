#!/usr/bin/env python2

from distutils.core import setup
from os.path import exists, abspath, dirname, join
import os
import sys

POSSIBLE_LIBS = ['libnetopt.so', 'netopt.dll', 'libnetopt.dll', 'libnetopt.dylib']

def find_path():
    lib_paths = [os.path.abspath('@CMAKE_LIBRARY_OUTPUT_DIRECTORY@'), 
    			 abspath(join(dirname(dirname(sys.argv[0])), 
    			 '../../../lib'))]
    
    found_path = None
    for path in lib_paths:
        for lib in POSSIBLE_LIBS:
            if exists(join(path,lib)):
                found_path = path
    if found_path is None:
    	return
    # create the __init__.py for netopt.lib
    open(os.path.join(found_path, '__init__.py'), 'a').close()
    	
    return found_path

setup(name='network-ilp',
      version='0.1',
      description='Network Optimization with Integer Programming',
      author='Markus Rempfler',
      author_email='markus.rempfler@tum.de',
      license='BSD',
      packages=['netopt', 'netopt.lib'],
      package_dir={'' : '@CMAKE_CURRENT_SOURCE_DIR@',
                  'netopt.lib': find_path()},
      package_data={'netopt.lib': POSSIBLE_LIBS}
)