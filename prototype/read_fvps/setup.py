from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize("cwrite_fvfps.pyx","cread_fvfps.pyx"),
    include_dirs=[numpy.get_include()]
)