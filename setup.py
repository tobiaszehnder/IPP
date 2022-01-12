from distutils.core import setup
from distutils.extension import Extension
import numpy as np
import os

extra_compile_args = ['-std=c++17']
if os.environ.get("IPP_DEBUG"):
    extra_compile_args.append('-O0')

ipp_extension = Extension('ipp',
                          ['ippmodule.cpp', 'ipp.cpp'],
                          include_dirs=[np.get_include()],
                          extra_compile_args=extra_compile_args)
setup(name='ipp',
      version='1.0',
      description='This is the IPP package',
      ext_modules=[ipp_extension])
