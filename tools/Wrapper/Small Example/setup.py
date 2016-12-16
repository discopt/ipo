from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    "Example",                 # name of extension
    ["Foo.cpp", "Foo2.cpp", "Example.pyx"],
    include_dirs = ["/opt/scipoptsuite-3.2.1/soplex-2.2.1/src"],
    libraries=["soplex-2.2.1.linux.x86_64.gnu.opt", "gmp", "gmpxx","z"],
    library_dirs=["/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/"],
    runtime_library_dirs = ["/opt/scipoptsuite-3.2.1/soplex-2.2.1/lib/"],
    extra_compile_args=['-std=c++11'],
    language="c++"              # this causes Pyrex/Cython to create C++ source
    )

setup(name = 'Example', version = '0.1',
    ext_modules = cythonize([ext])
    #ext_modules = [ext]
)
