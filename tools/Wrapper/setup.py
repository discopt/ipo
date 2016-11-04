###################################
# setup.py, 1.0
# 14.10.2016, Sandra Hicks
#	generiert den gcc Command, um die Cython files zu kompilieren
#
###################################
from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    "IPOVector",                 # name of extension
    ["/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/IPOVector.pyx", "/home/sandra/Documents/HiWi/IPO/ipo/tools/Wrapper/cppIPOVector.pxd", "/home/sandra/Documents/HiWi/IPO/ipo/ipo/vectors.cpp"],           # filename of our Pyrex/Cython source
    language="c++"              # this causes Pyrex/Cython to create C++ source
    )

setup(name = 'IPOVector', version = '0.1',
    ext_modules = cythonize([ext])
    #ext_modules = [ext]
)
