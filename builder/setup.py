import os

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

from lsst.ts.wep.Utility import getModulePath

# Get the path of module
modulePath = getModulePath()

extension = Extension(
           "cyMath",
           sources = [os.path.join(modulePath, "python", "lsst", "ts", "wep", "cwfs", "include", "cyMath.pyx")], 
           include_dirs = [numpy.get_include()], # Use numpy
)

setup(
        cmdclass = {"build_ext": build_ext},
        ext_modules = cythonize(extension),
)