###################################
# IPOLinearConstraint.pyx, 1.0
# 14.10.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse LinearConstraint aus der IPO Lib anzusprechen
#
###################################
cimport cppIPOLinearConstraint
cimport cppIPOVector
import IPOVector
from cpython.object cimport Py_EQ
from cython.operator cimport dereference as deref

cdef object CreateLinearConstraint(cppIPOLinearConstraint.LinearConstraint *linconst):
    py_linconst = IPOLinearConstraint()
    py_linconst.lin = linconst
    return py_linconst

cdef class IPOLinearConstraint:
    cdef cppIPOLinearConstraint.LinearConstraint *lin

    def __cinit__(self):
        self.lin = new cppIPOLinearConstraint.LinearConstraint()
        if self.lin is NULL:
            raise MemoryError()
    def __dealloc__(self):
        del self.lin

    def __richcmp__(IPOLinearConstraint self, IPOLinearConstraint y not None, int op):
        if op == Py_EQ:
            return self.lin==y.lin

    def isEquation(self):
        return self.lin.isEquation()

    def type(self):
        return self.lin.type()

    def normal(self):
        cdef cppIPOVector.Vector *c_vector
        c_vector = <cppIPOVector.Vector *>self.lin.normal()
        py_vector = cppIPOVector.CreateIPOVector(c_vector)
        
        return py_vector
        
    def definesCompleteFace(self):
        return self.lin.definesCompleteFace()

    def definesEmptyFace(self):
        return self.lin.definesEmptyFace()

    def definesTrivialFace(self):
        return self.lin.definesTrivialFace()

    def evaluatePoint(self, cppIPOVector.IPOVector point):
        return self.lin.evaluatePoint((<cppIPOVector.Vector &>point.vec))

    def evaluateRay(self, cppIPOVector.IPOVector ray):
        return self.lin.evaluateRay((<cppIPOVector.Vector &>ray.vec))

