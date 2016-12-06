###################################
# IPOSpace.pyx, 1.0
# 25.11.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse Vector aus der IPO Lib anzusprechen
#
###################################
cimport cppIPOSpace
cimport cppIPOVector
cimport cppIPOLinearConstraint

import IPOVector
import IPOLinearConstraint
import stdio
import IPOErrors

from libc cimport stdio
from cython.operator cimport dereference as deref, address as ref
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE
    
cdef object CreateIPOSpace(cppIPOSpace.Space *space):
    py_space = IPOSpace()
    py_space.cpp_space = space
    return py_space

cdef class IPOSpace:
    cdef cppIPOSpace.Space *cpp_space

    def __cinit__(self):
        self.cpp_space = new cppIPOSpace.Space()
        if self.cpp_space is NULL:
            raise MemoryError()

    def __dealloc__(self):
        del self.cpp_space

    def dimension(self):
        return self.space.dimension()

    def printVector(self, stream, cppIPOVector.IPOVector vector):
        if (not vector.isConstant()):
            raise IPOErrors.NonConstError('IPOVector')
    #cout, cerr, clog are ostreams
        if stream == "cout":
            #call function with cout
            print("debug")
        elif stream == "cerr":
            #call function with cerr
            print("debug")
        elif stream == "clog":
            #call function with clog
            print("debug")

    def printLinearForm(self, cppIPOVector.IPOVector vector):
        if (not vector.isConstant()):
            raise IPOErrors.NonConstError('IPOVector')
        print("lin")

    def printLinearConstraint(self, cppIPOLinearConstraint.IPOLinearConstraint lincons):
        if (not lincons.isConstant()):
            raise IPOErrors.NonConstError('IPOLinearConstraint')
        print("cons")

    def __getitem__(self, key):
        if type(key) is int:
            #cdef size_t var = key
            return self.space[key]
        else:
            raise TypeError()
        

    def __richcmp__(IPOSpace self, IPOSpace y not None, int op):
        if op == Py_EQ:
            return self.space==y.space
        elif op == Py_NE:
            return self.space!=y.space
        else:
            assert False

