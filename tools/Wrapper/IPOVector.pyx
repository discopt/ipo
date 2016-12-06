###################################
# IPOVector.pyx, 1.0
# 14.10.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse Vector aus der IPO Lib anzusprechen
#
###################################
cimport cppIPOVector
import IPOErrors

from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE
from cython.operator cimport dereference as deref, address as ref


cdef class IPOReferenceCountedVector:
    cdef cppIPOVector.ReferenceCountedVector *vec
    cdef const cppIPOVector.ReferenceCountedVector *const_vec
    
cdef object CreateIPOVector(cppIPOVector.Vector *vector):
    py_vector = IPOVector(False)
    py_vector.vec = vector
    return py_vector

cdef object CreateConstIPOVector(const cppIPOVector.Vector *vector):
    py_vector = IPOVector(True)
    py_vector.const_vec = vector
    return py_vector

cdef class IPOVector:
    cdef cppIPOVector.Vector *vec
    cdef const cppIPOVector.Vector *const_vec

    def __cinit__(self, isConst):
        if(not isConst):
            self.vec = new cppIPOVector.Vector()
            if (self.vec is NULL):
                raise MemoryError()
    def __dealloc__(self):
        del self.vec
        del self.const_vec

    def swap(self, IPOVector other):
        if(self.vec is not NULL):
            self.vec.swap(deref(other.vec))

    def isConstant(self):
        return (self.vec is NULL)

    #Parent Class methods
    def size(self):
        if(self.vec is not NULL):
            return (<cppIPOVector.ReferenceCountedVector*>(self.vec)).size()
        else:
            return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec)).size()

    def isSorted(self):
        if(self.vec is not NULL):
            return (<cppIPOVector.ReferenceCountedVector*>self.vec).isSorted()
        else:
            return (<cppIPOVector.ReferenceCountedVector*>self.const_vec).isSorted()

    def index(self, position):
        if(self.vec is not NULL):
            return (<cppIPOVector.ReferenceCountedVector*>self.vec).index(position)
        else:
            return (<cppIPOVector.ReferenceCountedVector*>self.const_vec).index(position)

    def approximation(self, position):
        if(self.vec is not NULL):
            return (<cppIPOVector.ReferenceCountedVector*>self.vec).approximation(position)
        else:
            return (<cppIPOVector.ReferenceCountedVector*>self.const_vec).approximation(position)

    def __richcmp__(IPOVector self, IPOReferenceCountedVector y not None, int op):
        if not y.isConstant():
            if not (self.vec is NULL):
                if op == Py_LT:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.vec))<y.vec
                elif op == Py_EQ:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.vec))==y.vec
                elif op == Py_NE:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.vec))!=y.vec
                else:
                    assert False
            else:
                if op == Py_LT:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec))<y.vec
                elif op == Py_EQ:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec))==y.vec
                elif op == Py_NE:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec))!=y.vec
                else:
                    assert False

        else:
            if not (self.vec is NULL):
                if op == Py_LT:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.vec))<y.const_vec
                elif op == Py_EQ:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.vec))==y.const_vec
                elif op == Py_NE:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.vec))!=y.const_vec
                else:
                    assert False
            else:
                if op == Py_LT:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec))<y.const_vec
                elif op == Py_EQ:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec))==y.const_vec
                elif op == Py_NE:
                    return (<cppIPOVector.ReferenceCountedVector*>(self.const_vec))!=y.const_vec
                else:
                    assert False
