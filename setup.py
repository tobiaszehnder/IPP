import os
from distutils.core import setup
from distutils.extension import Extension

import numpy as np


extra_compile_args = ["-std=c++17"]
if os.environ.get("IPP_DEBUG"):
    extra_compile_args.append("-O0")

ipp_extension = Extension(
    "ipp_cpp",
    sources=["cpp/src/ippmodule.cpp", "cpp/src/ipp.cpp"],
    include_dirs=[np.get_include(), "cpp/include"],
    extra_compile_args=extra_compile_args,
)
setup(
    name="ipp",
    version="1.0",
    description=(
        "Interspecies Point Projection - A tool for comparative genomics "
        "beyond direct sequence alignments"
    ),
    ext_modules=[ipp_extension],
)
