###################################
# IPOVector.pyx, 1.0
# 14.10.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse Vector aus der IPO Lib anzusprechen
#
###################################
cimport cppIPOVector
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE

cdef class IPOReferenceCountedVector:
    cdef cppIPOVector.ReferenceCountedVector *vec
    
cdef object CreateIPOVector(cppIPOVector.Vector *vector):
    py_vector = IPOVector()
    py_vector.vec = vector
    return py_vector

cdef class IPOVector:
    cdef cppIPOVector.Vector *vec

    def __cinit__(self):
        self.vec = new cppIPOVector.Vector()
        if self.vec is NULL:
            raise MemoryError()
    def __dealloc__(self):
        del self.vec
    def __call__(self):
        return self._thisptr[0]()
    def swap(self):
        print("swap")

    #Parent Class methods
    def size(self):
        return (<cppIPOVector.ReferenceCountedVector*>(self.vec)).size()

    def isSorted(self):
        return (<cppIPOVector.ReferenceCountedVector*>self.vec).isSorted()

    def index(self, position):
        return (<cppIPOVector.ReferenceCountedVector*>self.vec).index(position)

    def approximation(self, position):
        return (<cppIPOVector.ReferenceCountedVector*>self.vec).approximation(position)

    def __richcmp__(IPOVector self, IPOReferenceCountedVector y not None, int op):
        if op == Py_LT:
            return (<cppIPOVector.ReferenceCountedVector*>(self.vec))<y.vec
        elif op == Py_EQ:
            return (<cppIPOVector.ReferenceCountedVector*>(self.vec))==y.vec
        elif op == Py_NE:
            return (<cppIPOVector.ReferenceCountedVector*>(self.vec))!=y.vec
        else:
            assert False

