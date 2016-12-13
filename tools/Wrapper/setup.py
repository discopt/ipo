###################################
# setup.py, 1.0
# 14.10.2016, Sandra Hicks
#	generiert den gcc Command, um die Cython files zu kompilieren
#
###################################
from setuptools import setup, Extension
from Cython.Build import cythonize
import os
from distutils.util import get_platform

ext = Extension(
    "IPO",                 # name of extension

    ["cppSoplexRational.pxd", "SoplexRational.pyx",  "cppIPOVector.pxd", "cppIPOSpace.pxd", "cppIPOLinearConstraint.pxd", "cppScipOracle.pxd", "IPOLinearConstraint.pyx","IPOVector.pyx",  "IPOScipOracle.pyx", "IPOSpace.pyx", "IPOErrors.pyx", "IPO.pyx"],           # filename of Cython source

    include_dirs = ["/usr/local/include/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/src", "/opt/scipoptsuite-3.2.1/scip-3.2.1/src"],

    libraries = ["ipo", "scipopt-3.2.1.linux.x86_64.gnu.opt", "soplex-2.2.1.linux.x86_64.gnu.opt", "scip-3.2.1.linux.x86_64.gnu.opt", "zimpl-3.3.3.linux.x86_64.gnu.opt", "nlpi.cppad.linux.x86_64.gnu.opt", "objscip-3.2.1.linux.x86_64.gnu.opt", "lpispx-3.2.1.linux.x86_64.gnu.opt", "gmp"],

    library_dirs = ["/usr/local/lib/", "/opt/scipoptsuite-3.2.1/lib/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/", "/opt/scipoptsuite-3.2.1/scip-3.2.1/lib", "/opt/scipoptsuite-3.2.1/zimpl-3.3.3/lib"],

    runtime_library_dirs = ["/usr/local/lib/", "/opt/scipoptsuite-3.2.1/lib/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/", "/opt/scipoptsuite-3.2.1/scip-3.2.1/lib", "/opt/scipoptsuite-3.2.1/zimpl-3.3.3/lib"],

    #needs C++11 standard to compile
    extra_compile_args=['-std=c++11', '-lz'],
    language="c++"
    )

setup(name = 'IPOVector', version = '0.1',
    ext_modules = cythonize([ext])
)
