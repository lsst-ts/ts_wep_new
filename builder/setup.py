import os

from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

from lsst.ts.wep.Utility import getModulePath

# Get the path of module
modulePath = getModulePath()

sourceFilePath = os.path.join(
    modulePath, "python", "lsst", "ts", "wep", "cwfs", "include", "cyMath.pyx"
)

# Use numpy
extension = Extension(
    "cyMath", sources=[sourceFilePath], include_dirs=[numpy.get_include()],
)

setup(
    cmdclass={"build_ext": build_ext}, ext_modules=cythonize(extension),
)
