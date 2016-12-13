###################################
# IPOLinearConstraint.pyx, 1.0
# 14.10.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse LinearConstraint aus der IPO Lib anzusprechen
#
###################################
cimport cppIPOLinearConstraint
cimport cppIPOVector
cimport cppSoplexRational
#import IPOVector
import IPOErrors
from cpython.object cimport Py_EQ, Py_LT
from cython.operator cimport dereference as deref, address as ref

cdef object CreateLinearConstraint(cppIPOLinearConstraint.LinearConstraint *linconst):
    py_linconst = IPOLinearConstraint()
    py_linconst.lin = linconst
    return py_linconst

cdef class IPOLinearConstraint:
    cdef cppIPOLinearConstraint.LinearConstraint *lin
    cdef const cppIPOLinearConstraint.LinearConstraint *const_lin

    def __cinit__(self, isConst):
        if (not isConst):
            self.lin = new cppIPOLinearConstraint.LinearConstraint()
            if self.lin is NULL:
                raise MemoryError()
    def __dealloc__(self):
        del self.lin
        del self.const_lin

    def isConstant(self):
        return (self.const_lin is not NULL)

    def __richcmp__(IPOLinearConstraint self, IPOLinearConstraint y not None, int op):
        if op == Py_EQ:
            return self.lin==y.lin
        elif op == Py_LT:
            return self.lin<y.lin

    def isEquation(self):
        if(self.lin is not NULL):
            return self.lin.isEquation()
        else:
            return self.const_lin.isEquation()

    def type(self):
        if(self.lin is not NULL):
            return self.lin.type()
        else:
            return self.const_lin.type()

    def normal(self):
        cdef const cppIPOVector.Vector *c_vector
        if(self.lin is not NULL):
            c_vector = ref(self.lin.normal())
            py_vector = cppIPOVector.CreateConstIPOVector(c_vector)
            return py_vector
        else:
            c_vector = ref(self.const_lin.normal())
            py_vector = cppIPOVector.CreateConstIPOVector(c_vector)
            return py_vector

    def rhs(self):
        cdef const cppSoplexRational.Rational *rational
        if(self.lin is not NULL):
            rational = ref(self.lin.rhs())
            py_rational = cppSoplexRational.CreateConstSoplexRational(rational)
            return py_rational
        else:
            rational = ref(self.const_lin.rhs())
            py_rational = cppSoplexRational.CreateConstSoplexRational(rational)
            return py_rational

    def getMaximumNorm(self):
        cdef cppSoplexRational.Rational *rational
        #if(self.lin is not NULL):
            #rational = &self.lin.getMaximumNorm()
            #py_rational = cppSoplexRational.CreateSoplexRational(rational)
            #return py_rational
        #else:
            #rational = &self.const_lin.getMaximumNorm()
            #py_rational = cppSoplexRational.CreateSoplexRational(rational)
            #return py_rational

    def definesCompleteFace(self):
        if(self.lin is not NULL):
            return self.lin.definesCompleteFace()
        else:
            return self.const_lin.definesCompleteFace()

    def definesEmptyFace(self):
        if(self.lin is not NULL):
            return self.lin.definesEmptyFace()
        else:
            return self.const_lin.definesEmptyFace()

    def definesTrivialFace(self):
        if(self.lin is not NULL):
            return self.lin.definesTrivialFace()
        else:
            return self.const_lin.definesTrivialFace()

    def evaluatePoint(self, cppIPOVector.IPOVector point):
        if (point.isConstant()):
            if(self.lin is not NULL):
                return self.lin.evaluatePoint(deref(point.const_vec))
            else:
                return self.const_lin.evaluatePoint(deref(point.const_vec))
        else:
            raise IPOErrors.NonConstError('IPOVector')

    def evaluateRay(self, cppIPOVector.IPOVector ray):
        if (ray.isConstant()):
            
            if(self.lin is not NULL):
                return self.lin.evaluateRay(deref(ray.const_vec))
            else:
                return self.const_lin.evaluateRay(deref(ray.const_vec))
        else:
            raise IPOErrors.NonConstError('IPOVector')
