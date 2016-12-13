from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    "Example",                 # name of extension
    ["Foo.cpp", "cppExample.pxd", "Example.pyx"],
    extra_compile_args=['-std=c++11'],
    language="c++"              # this causes Pyrex/Cython to create C++ source
    )

setup(name = 'Example', version = '0.1',
    ext_modules = cythonize([ext])
    #ext_modules = [ext]
)
