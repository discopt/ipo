###################################
# IPOSpace.pyx, 1.0
# 25.11.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse Vector aus der IPO Lib anzusprechen
#
###################################
cimport cppSoplexRational

from libc cimport stdio
from cython.operator cimport dereference as deref, address as ref
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE
    
cdef object CreateSoplexRational(cppSoplexRational.Rational *rational):
    py_rational = SoplexRational()
    py_rational.cpp_rational = rational
    return py_rational

cdef object CreateConstSoplexRational(const cppSoplexRational.Rational *rational):
    py_rational = SoplexRational()
    py_rational.const_rational = rational
    return py_rational



cdef class SoplexRational:
    cdef cppSoplexRational.Rational *cpp_rational
    cdef const cppSoplexRational.Rational *const_rational

    def __cinit__(self, isConst):
        if(not isConst):
            self.cpp_rational = new cppSoplexRational.Rational()
            if self.cpp_rational is NULL:
                raise MemoryError()

    def __dealloc__(self):
        del self.cpp_rational
        del self.const_rational
