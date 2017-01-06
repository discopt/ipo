## @package IPO
#  Documentation for this module.
#
#  More details.

####################################
#Imports
cimport cppIPOexact as cppIPO

cdef extern from "stdlib.h":
  void free(void* ptr)

import ctypes
from libc cimport stdio
from cython.operator cimport dereference as deref, address as ref
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE
from libcpp.memory cimport shared_ptr
####################################
#Errors

## NonConstError
# Error Class to throw errors if a pointer is used instead of a constant pointer to a C++ class.
#  
class NonConstError:
    def init(self, value):
        self.value = 'This is no const value: '+value

    def __str__(self):
        return repr(self.value)

####################################
#Soplex Rational

## Wrapper class for the Rational class from Soplex
#
#
cdef class SoplexRational:
    cdef cppIPO.Rational *cpp_rational
    cdef const cppIPO.Rational *const_rational
    ## Constructor
    #  @param isConst indicates if the C++ pointer for the object is a const pointer
    #
    def __cinit__(self, isConst):
        if(not isConst):
            self.cpp_rational = new cppIPO.Rational()
            if self.cpp_rational is NULL:
                raise MemoryError()

    def __dealloc__(self):
        del self.cpp_rational
        del self.const_rational

## Creates a Python SoplexRational from a pointer to a C++ Rational.
#
cdef object CreateSoplexRational(cppIPO.Rational *rational):
    py_rational = SoplexRational(False)
    py_rational.cpp_rational = rational
    return py_rational

## Creates a Python SoplexRational from a const pointer to a C++ Rational.
#
cdef object CreateConstSoplexRational(const cppIPO.Rational *rational):
    py_rational = SoplexRational(True)
    py_rational.const_rational = rational
    return py_rational

####################################
#IPO Vector

## Base class for IPO Vector class
#
cdef class IPOReferenceCountedVector:
    cdef cppIPO.ReferenceCountedVector *vec
    cdef const cppIPO.ReferenceCountedVector *const_vec

## Vector class
#
cdef class IPOVector:
    cdef cppIPO.Vector *vec
    cdef const cppIPO.Vector *const_vec

    ## Constructor
    #  @param isConst indicates if the C++ pointer for the object is a const pointer
    #
    def __cinit__(self, isConst):
        if(not isConst):
            self.vec = new cppIPO.Vector()
            if (self.vec is NULL):
                raise MemoryError()
    def __dealloc__(self):
        if(self.vec is not NULL):
            del self.vec
        if(self.const_vec is not NULL):
            del self.const_vec

    def swap(self, IPOVector other):
        if(self.vec is not NULL):
            self.vec.swap(deref(other.vec))

    def isConstant(self):
        return (self.vec is NULL)

    #Parent Class methods
    def size(self):
        if(self.vec is not NULL):
            return (<cppIPO.ReferenceCountedVector*>(self.vec)).size()
        else:
            return (<cppIPO.ReferenceCountedVector*>(self.const_vec)).size()

    def isSorted(self):
        if(self.vec is not NULL):
            return (<cppIPO.ReferenceCountedVector*>self.vec).isSorted()
        else:
            return (<cppIPO.ReferenceCountedVector*>self.const_vec).isSorted()

    def index(self, position):
        if(self.vec is not NULL):
            return (<cppIPO.ReferenceCountedVector*>self.vec).index(position)
        else:
            return (<cppIPO.ReferenceCountedVector*>self.const_vec).index(position)

    def approximation(self, position):
        if(self.vec is not NULL):
            return (<cppIPO.ReferenceCountedVector*>self.vec).approximation(position)
        else:
            return (<cppIPO.ReferenceCountedVector*>self.const_vec).approximation(position)

    ## Special Method for operator overloading
    # this method overloads the operators ==, != and <
    # @return boolean
    def __richcmp__(IPOVector self, IPOVector y not None, int op):
        if not y.isConstant():
            if not (self.vec is NULL):
                if op == Py_LT:
                    return (<cppIPO.ReferenceCountedVector*>(self.vec))<(<cppIPO.ReferenceCountedVector*>(y.vec))
                elif op == Py_EQ:
                    return (<cppIPO.ReferenceCountedVector*>(self.vec))==(<cppIPO.ReferenceCountedVector*>(y.vec))
                elif op == Py_NE:
                    return (<cppIPO.ReferenceCountedVector*>(self.vec))!=(<cppIPO.ReferenceCountedVector*>(y.vec))
                else:
                    assert False
            else:
                if op == Py_LT:
                    return (<cppIPO.ReferenceCountedVector*>(self.const_vec))<(<cppIPO.ReferenceCountedVector*>(y.vec))
                elif op == Py_EQ:
                    return (<cppIPO.ReferenceCountedVector*>(self.const_vec))==(<cppIPO.ReferenceCountedVector*>(y.vec))
                elif op == Py_NE:
                    return (<cppIPO.ReferenceCountedVector*>(self.const_vec))!=(<cppIPO.ReferenceCountedVector*>(y.vec))
                else:
                    assert False

        else:
            if not (self.vec is NULL):
                if op == Py_LT:
                    return (<cppIPO.ReferenceCountedVector*>(self.vec))<(<cppIPO.ReferenceCountedVector*>(y.const_vec))
                elif op == Py_EQ:
                    return (<cppIPO.ReferenceCountedVector*>(self.vec))==(<cppIPO.ReferenceCountedVector*>(y.const_vec))
                elif op == Py_NE:
                    return (<cppIPO.ReferenceCountedVector*>(self.vec))!=(<cppIPO.ReferenceCountedVector*>(y.const_vec))
                else:
                    assert False
            else:
                if op == Py_LT:
                    return (<cppIPO.ReferenceCountedVector*>(self.const_vec))<(<cppIPO.ReferenceCountedVector*>(y.const_vec))
                elif op == Py_EQ:
                    return (<cppIPO.ReferenceCountedVector*>(self.const_vec))==(<cppIPO.ReferenceCountedVector*>(y.const_vec))
                elif op == Py_NE:
                    return (<cppIPO.ReferenceCountedVector*>(self.const_vec))!=(<cppIPO.ReferenceCountedVector*>(y.const_vec))
                else:
                    assert False

## Creates a Python IPOVector from a pointer to a C++ ipo::Vector.
#
cdef object CreateIPOVector(cppIPO.Vector *vector):
    py_vector = IPOVector(False)
    py_vector.vec = vector
    return py_vector
## Creates a Python IPOVector from a const pointer to a C++ ipo::Vector.
#
cdef object CreateConstIPOVector(const cppIPO.Vector *vector):
    py_vector = IPOVector(True)
    py_vector.const_vec = vector
    return py_vector

####################################
#IPO Linear Constraint

cdef class IPOLinearConstraint:
    cdef cppIPO.LinearConstraint *lin
    cdef const cppIPO.LinearConstraint *const_lin

    ## Constructor
    #  @param isConst indicates if the C++ pointer for the object is a const pointer
    #
    def __cinit__(self, isConst):
        if (not isConst):
            self.lin = new cppIPO.LinearConstraint()
            if self.lin is NULL:
                raise MemoryError()
    def __dealloc__(self):
        del self.lin
        del self.const_lin

    def isConstant(self):
        return (self.const_lin is not NULL)

    ## Special Method for operator overloading
    # this method overloads the operators == and <
    # @return boolean
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
        cdef const cppIPO.Vector *c_vector
        if(self.lin is not NULL):
            c_vector = ref(self.lin.normal())
            py_vector = CreateConstIPOVector(c_vector)
            return py_vector
        else:
            c_vector = ref(self.const_lin.normal())
            py_vector = CreateConstIPOVector(c_vector)
            return py_vector

    def rhs(self):
        cdef const cppIPO.Rational *rational
        if(self.lin is not NULL):
            rational = ref(self.lin.rhs())
            py_rational = CreateConstSoplexRational(rational)
            return py_rational
        else:
            rational = ref(self.const_lin.rhs())
            py_rational = CreateConstSoplexRational(rational)
            return py_rational

    def getMaximumNorm(self):
        cdef cppIPO.Rational *rational
        #if(self.lin is not NULL):
            #rational = &self.lin.getMaximumNorm()
            #py_rational = CreateSoplexRational(rational)
            #return py_rational
        #else:
            #rational = &self.const_lin.getMaximumNorm()
            #py_rational = CreateSoplexRational(rational)
            #return py_rational
        print "Not yet implemented - Non-l-value Error.\n"

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

    def evaluatePoint(self, IPOVector point):
        if (point.isConstant()):
            if(self.lin is not NULL):
                return self.lin.evaluatePoint(deref(point.const_vec))
            else:
                return self.const_lin.evaluatePoint(deref(point.const_vec))
        else:
            raise NonConstError('IPOVector')

    def evaluateRay(self, IPOVector ray):
        if (ray.isConstant()):
            
            if(self.lin is not NULL):
                return self.lin.evaluateRay(deref(ray.const_vec))
            else:
                return self.const_lin.evaluateRay(deref(ray.const_vec))
        else:
            raise NonConstError('IPOVector')

## Creates a Python IPOLinearConstraint from a pointer to a C++ ipo::LinearConstraint.
#
cdef object CreateLinearConstraint(cppIPO.LinearConstraint *linconst):
    py_linconst = IPOLinearConstraint(False)
    py_linconst.lin = linconst
    return py_linconst
## Creates a Python IPOLinearConstraint from a const pointer to a C++ ipo::LinearConstraint.
#
cdef object CreateConstLinearConstraint(const cppIPO.LinearConstraint *linconst):
    py_linconst = IPOLinearConstraint(True)
    py_linconst.const_lin = linconst
    return py_linconst

####################################
#InnerDesciption/OuterDescription

## IPOInnerDescription is a tupel of lists containing IPOVectors
#
class IPOInnerDescription:
    def __init__(self):
        self.dpoints = []
        self.drays = []

    ## getPoints
    # @return points list of IPOVectors
    def getPoints(self):
        return self.dpoints

    ## getRays
    # @retrun rays list of IPOVectors
    def getRays(self):
        return self.drays

## IPOAffineOuterDescription is a list of IPOLinearConstraints
#
class IPOAffineOuterDescription:
    def __init__(self):
        self.dconstraints = []
    ## getConstraints
    # @return list of IPOLinearConstraints
    def getConstraints(self):
        return self.dconstraints

####################################
#IPO Space

cdef class IPOSpace:
    cdef cppIPO.Space *cpp_space
    cdef const cppIPO.Space *const_space

    ## Constructor
    #  @param isConst indicates if the C++ pointer for the object is a const pointer
    #
    def __cinit__(self, isConst, isNew):
        if(isConst is False and isNew):
            self.cpp_space = new cppIPO.Space()
            if self.cpp_space is NULL:
                raise MemoryError()

    def __dealloc__(self):
        if(self.cpp_space is not NULL):
            del self.cpp_space
        else:
            del self.const_space

    def isConstant(self):
        return (self.const_space is not NULL)

    ## Dimension of the Space
    # @return int dimension
    def dimension(self):
        if(self.cpp_space is not NULL):
            return self.cpp_space.dimension()
        else:
            return self.const_space.dimension()

    def printVector(self, IPOVector vector):
        self.cpp_space.printVector(cppIPO.cout, deref(vector.vec))

    def printLinearForm(self, IPOVector vector):
        self.cpp_space.printLinearForm(cppIPO.cout, deref(vector.vec))

    def printLinearConstraint(self, IPOLinearConstraint lincons):
        self.cpp_space.printLinearConstraint(cppIPO.cout, deref(lincons.lin))

    ## Special Method for operator overloading
    #  this method overloads the operator [] for indexing
    #  @return const string
    def __getitem__(self, int key):
        if type(key) is int:
            if(self.cpp_space is not NULL):
                return deref(self.cpp_space)[key]
            else:
                return deref(self.const_space)[key]
        else:
            raise TypeError()
        
    ## Special Method for operator overloading
    # this method overloads the operators == and !=
    # @return boolean
    def __richcmp__(IPOSpace self, IPOSpace y not None, int op):
        if(self.cpp_space is not NULL):
            if(y.isConstant()):
                if op == Py_EQ:
                    return self.cpp_space==y.const_space
                elif op == Py_NE:
                    return self.cpp_space!=y.const_space
                else:
                    assert False
            else:
                if op == Py_EQ:
                    return self.cpp_space==y.cpp_space
                elif op == Py_NE:
                    return self.cpp_space!=y.cpp_space
                else:
                    assert False
        else:
            if(y.isConstant()):
                if op == Py_EQ:
                    return self.const_space==y.const_space
                elif op == Py_NE:
                    return self.const_space!=y.const_space
                else:
                    assert False
            else:
                if op == Py_EQ:
                    return self.const_space==y.cpp_space
                elif op == Py_NE:
                    return self.const_space!=y.cpp_space
                else:
                    assert False

## Creates a Python IPOSpace from a pointer to a C++ ipo::Space.
#
cdef object CreateIPOSpace(cppIPO.Space *space):
    py_space = IPOSpace(False)
    py_space.cpp_space = space
    return py_space

## Creates a Python IPOSpace from aconst  pointer to a C++ ipo::Space.
#
cdef object CreateConstIPOSpace(const cppIPO.Space *space):
    py_space = IPOSpace(True)
    py_space.const_space = space
    return py_space

####################################
#IPO OracleBase
cdef class IPOOracleBase:
    cdef cppIPO.OracleControllerBase *base

####################################
#IPO ScipOracle
cdef class IPOScipOracle(IPOOracleBase):
    cdef cppIPO.ScipOracleController *oracle

    ## Constructor
    #  @param isNew indicates if the Oracle is a new one (True) or if it needs to set a next pointer (False)
    #
    def __cinit__(self, str name, isNew):
        self.oracle = NULL
        if(isNew == True):
            self.oracle = new cppIPO.ScipOracleController(name)

    def __dealloc__(self):
        if(self.oracle is not NULL):
            del self.oracle

    ## Name Property of the ScipOracle
    #
    def name(self):
        return self.oracle.name()

    ## Heuristic Level of the ScipOracle
    #
    def heuristicLevel(self):
        return self.oracle.heuristicLevel_ScipOracle()

    ## Heuristic Level of the CacheOracle
    #
    def heuristicLevel_CacheOracle(self):
        return self.oracle.heuristicLevel_CacheOracle()

    ## Affine Hull
    # @return (IPOInnerDescription, IPOAffineOuterDescription) Tupel of those classes
    def affineHull(self):
        #########-1-#########
        #convert inner description to 2-tupel IPOVector lists
        cdef cppIPO.InnerDescription c_inner = self.oracle.affineHullInner(1)

        #automagically convert to python list
        c_points = c_inner.points
        c_rays = c_inner.rays
        points = []
        rays = []
        cdef cppIPO.Vector *c_vector
        #convert points to python wrapperclass IPOVector
        if(c_points.size()>0):
            for i in range(0,c_points.size()):
                c_vector = ref(c_points[i])
                py_vector = CreateIPOVector(c_vector)
                points.append(py_vector)

        if(c_rays.size() > 0):
        #convert rays to python wrapperclass IPOVector
            for i in range(0,c_rays.size()):
                c_vector = ref(c_rays[i])
                py_vector = CreateIPOVector(c_vector)
                rays.append(py_vector)

        innerDescription = IPOInnerDescription()

        innerDescription.dpoints = points
        innerDescription.drays = rays
        #########-2-#########
        #convert outer description to IPOLinearConstraint list
        cdef cppIPO.AffineOuterDescription c_outer = self.oracle.affineHullOuter(1)
        outerDescription = IPOAffineOuterDescription()
        outer = []

        for i in range(0, c_outer.size()):
            c_linconst = ref(c_outer[i])
            py_linconst = CreateLinearConstraint(c_linconst)
            outer.append(py_linconst)

        outerDescription.dconstraints = outer

        return (innerDescription, outerDescription)

## Creates an IPOScipOracle from a filename and the pointer to another ScipOracleController
#  the ScipOracle of the ScipOracleController is set as the next pointer in the Oracle chain.
cdef object CreateScipOracle(str name, cppIPO.OracleControllerBase oracle):
    py_oracle = IPOScipOracle(name, 0)
    py_oracle.oracle = new cppIPO.ScipOracleController(name, oracle)
    return py_oracle

####################################
#IPO ExactScipOracle
cdef class IPOExactScipOracle(IPOOracleBase):
    cdef cppIPO.ExactScipOracleController *exactoracle

    ## Constructor
    #  @param isNew indicates if the Oracle is a new one (True) or if it needs to set a next pointer (False)
    #
    def __cinit__(self, str name, isNew):
        self.exactoracle = NULL
        cdef shared_ptr[cppIPO.MixedIntegerSet] mxip = cppIPO.getMixedIntegerSet(name)
        if(isNew == True):
            self.exactoracle = new cppIPO.ExactScipOracleController(name, mxip)

    def __dealloc__(self):
        if(self.exactoracle is not NULL):
            del self.exactoracle

    ## Name Property of the ScipOracle
    #
    def name(self):
        return self.exactoracle.name()

    ## Heuristic Level of the ScipOracle
    #
    def heuristicLevel(self):
        return self.exactoracle.heuristicLevel_ExactScipOracle()

    ## Heuristic Level of the CacheOracle
    #
    def heuristicLevel_CacheOracle(self):
        return self.exactoracle.heuristicLevel_CacheOracle()

cdef object CreateExactScipOracle(str name, cppIPO.OracleControllerBase oracle):
    cdef shared_ptr[cppIPO.MixedIntegerSet] mxip = cppIPO.getMixedIntegerSet(name)
    py_oracle = IPOExactScipOracle(name, 0)
    py_oracle.exactoracle = new cppIPO.ExactScipOracleController(name, mxip, oracle)
    return py_oracle
####################################
#IPO Polyhedron

cdef class IPOPolyhedron:
    cdef cppIPO.Polyhedron* poly

    ##Erstellung benoetigt ein Oracle
    ##affineHull benoetigt Handlers

    ## dimension
    # @return int dimension
    def dimension(self):
        return self.poly.dimension()

    def setAffineHullLastModerateHeuristicLevel(self, int level):
        self.poly.setAffineHullLastModerateHeuristicLevel(level)

    def setAffineHullLastCheapHeuristicLevel(self, int level):
        self.poly.setAffineHullLastCheapHeuristicLevel(level)

from cppIPO cimport Polyhedron

cdef class IPOFace:
    cdef Polyhedron.Face* face

####################################
#Mixed Integer Set

cdef class IPOMixedIntegerSet:
    cdef shared_ptr [cppIPO.MixedIntegerSet] mip
    def __cinit__(self, filename):
        mip = cppIPO.getMixedIntegerSet(filename)

    def space(self):
        cdef cppIPO.Space space = self.mip.get().space()
        py_space = IPOSpace(False, False)
        py_space.cpp_space = ref(space)
        return py_space

    def numVariables(self):
        return self.mip.get().numVariables()

    def numRows(self):
        return self.mip.get().numRows()

    def rowName(self, int row):
        return self.mip.get().rowName(row)

    def isIntegral(self, int param_variable):
        return self.mip.get().isIntegral(param_variable)

    def variable(self, int var):
        cdef cppIPO.Variable cpp_variable = self.mip.get().variable(var)
        py_upperBound = CreateSoplexRational(ref(cpp_variable.upperBound))
        py_lowerBound = CreateSoplexRational(ref(cpp_variable.lowerBound))
        py_variable = IPOVariable(cpp_variable.integral, py_upperBound, py_lowerBound)
        return py_variable

class IPOVariable:
    def __init__(self, _integral, SoplexRational _upperBound, SoplexRational _lowerBound):
        if(type(_integral) is not bool):
            raise TypeError
        integral = _integral
        upperBound = _upperBound
        lowerBound = _lowerBound

    def getIntegral(self):
        return self.integral
    def getUpperBound(self):
        return self.upperBound
    def getLowerBound(self):
        return self.lowerBound

####################################
#Affine Hull

def affineHull():
    return 0


####################################
#IPO Test

cdef class Test:
    cdef cppIPO.Foo* foo

    def __cinit__(self):
        self.foo = new cppIPO.Foo()
        if self.foo is NULL:
            raise MemoryError()

    def __dealloc__(self):
        del self.foo

    def print_Example(self):
        return self.foo.printFoo()
