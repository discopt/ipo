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

#platform = get_platform()
#if platform in ["linux-x86_64"]:
#    os.environ["CFLAGS"] = "-fpic"

ext = Extension(
    "IPO",                 # name of extension
    ["/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPO.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppSoplexRational.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/SoplexRational.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOVector.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppIPOVector.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOLinearConstraint.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppIPOLinearConstraint.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/CppObjects/ScipOracleController.cpp", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppScipOracle.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOScipOracle.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOSpace.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppIPOSpace.pxd"],           # filename of our Pyrex/Cython source
    include_dirs = ["/usr/local/include/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/src", "/opt/scipoptsuite-3.2.1/scip-3.2.1/src"],
    libraries = ["ipo", "scipopt-3.2.1.linux.x86_64.gnu.opt", "soplex-2.2.1.linux.x86_64.gnu.opt", "scip-3.2.1.linux.x86_64.gnu.opt", "zimpl-3.3.3.linux.x86_64.gnu.opt", "nlpi.cppad.linux.x86_64.gnu.opt", "objscip-3.2.1.linux.x86_64.gnu.opt", "lpispx-3.2.1.linux.x86_64.gnu.opt", "lpispx-3.2.1.linux.x86_64.gnu.opt"],
    library_dirs = ["/usr/local/lib/", "/opt/scipoptsuite-3.2.1/lib/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/", "/opt/scipoptsuite-3.2.1/scip-3.2.1/lib", "/opt/scipoptsuite-3.2.1/zimpl-3.3.3/lib"],
    runtime_library_dirs = ["/usr/local/lib/", "/opt/scipoptsuite-3.2.1/lib/", "/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/", "/opt/scipoptsuite-3.2.1/scip-3.2.1/lib", "/opt/scipoptsuite-3.2.1/zimpl-3.3.3/lib"],
    extra_compile_args=['-std=c++11', '-ggdb'],
    #extra_objects=["/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/CppObjects/ScipOracleController.o"],
    language="c++"              # this causes Pyrex/Cython to create C++ source
    )

setup(name = 'IPOVector', version = '0.1',
    ext_modules = cythonize([ext])
)
