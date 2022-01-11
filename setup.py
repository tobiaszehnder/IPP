from distutils.core import setup
from distutils.extension import Extension
import numpy as np

ipp_extension = Extension('ipp',
                          ['ippmodule.cpp', 'ipp.cpp'],
                          include_dirs=[np.get_include()],
                          extra_compile_args=['-std=c++17'])
setup(name='ipp',
      version='1.0',
      description='This is the IPP package',
      ext_modules=[ipp_extension])
