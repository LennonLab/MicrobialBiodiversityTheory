from distutils.core import setup
from Cython.Build import cythonize

setup(ext_modules = cythonize("HMP_EMP_ObsPred.pyx"))
