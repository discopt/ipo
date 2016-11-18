###################################
# setup.py, 1.0
# 14.10.2016, Sandra Hicks
#	generiert den gcc Command, um die Cython files zu kompilieren
#
###################################
from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    "IPO",                 # name of extension
    ["/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPO.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOVector.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppIPOVector.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/ipo/vectors.cpp", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOLinearConstraint.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppIPOLinearConstraint.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/ipo/linear_constraint.cpp", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppScipOracle.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOScipOracle.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/ScipOracleController.cpp"],           # filename of our Pyrex/Cython source
    libraries=["/opt/scipoptsuite-3.2.1/soplex-2.2.1/src","/home/sandra/Documents/HiWi/IPO/ipo/ipo/"],
    include_dirs=["/opt/scipoptsuite-3.2.1/soplex-2.2.1/src","/home/sandra/Documents/HiWi/IPO/ipo/"],
    language="c++"              # this causes Pyrex/Cython to create C++ source
    )

setup(name = 'IPOVector', version = '0.1',
    ext_modules = cythonize([ext])
    #ext_modules = [ext]
)
