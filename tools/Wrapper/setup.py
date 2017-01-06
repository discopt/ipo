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

scipExactInstalled = 0

scipLib = "scip-3.2.1.linux.x86_64.gnu.opt"
scipPath = "/opt/scipoptsuite-3.2.1/scip-3.2.1/lib/"
objscipLib = "objscip-3.2.1.linux.x86_64.gnu.opt"
objscipPath = scipPath
nlpiLib = "nlpi.cppad.linux.x86_64.gnu.opt"
nlpiPath = scipPath
lpiLib = "lpispx-3.2.1.linux.x86_64.gnu.opt"
lpiPath = scipPath
soplexLib = "soplex-2.2.1.linux.x86_64.gnu.opt"
soplexPath = "/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/"
zimplLib = "zimpl-3.3.3.linux.x86_64.gnu.opt"
zimplPath = "/opt/scipoptsuite-3.2.1/zimpl-3.3.3/lib"

if(scipExactInstalled == 0):
    ext = Extension(
        "IPO",                 # name of extension

        #["cppSoplexRational.pxd", "SoplexRational.pyx",  "cppIPOVector.pxd", "cppIPOSpace.pxd", "cppIPOLinearConstraint.pxd", "cppScipOracle.pxd", "IPOLinearConstraint.pyx","IPOVector.pyx",  "IPOScipOracle.pyx", "IPOSpace.pyx", "IPOErrors.pyx", "IPO.pyx"],           # filename of Cython source
        ["IPO.pyx"],           # filename of Cython source

        include_dirs = ["/usr/local/include/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/src/", "/opt/scipoptsuite-3.2.1/scip-3.2.1/src/"],

        libraries = ["ipo", scipLib, soplexLib, zimplLib, nlpiLib, objscipLib, lpiLib, "gmp", "gmpxx", "z"],

        library_dirs = ["/usr/local/lib/", scipPath, soplexPath, zimplPath],

        runtime_library_dirs = ["/usr/local/lib/", "/opt/scipoptsuite-3.2.1/lib/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/", "/opt/scipoptsuite-3.2.1/scip-3.2.1/lib/", "/opt/scipoptsuite-3.2.1/zimpl-3.3.3/lib/"],

        #needs C++11 standard to compile
        extra_compile_args=['-std=c++11'],
        language="c++"
    )

    setup(name = 'IPO', version = '0.1',
        ext_modules = cythonize([ext])
    )
