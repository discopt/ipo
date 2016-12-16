####################################
#Imports
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector

####################################
#Soplex Rational
cdef extern from "rational.h" namespace "soplex":
    cdef cppclass Rational:
        Rational() except +

####################################
#IPO Vector


cdef extern from "ipo/vectors.h" namespace "ipo":
    cdef cppclass MutableVector:
        MutableVector() except +
    cdef cppclass Vector:

        Vector() except +
        void swap(Vector& other)

        #Vector& operator=(const Vector& other)
        #Vector& operator=(const MutableVector& other)
        #not possible in python to overload


    cdef cppclass ReferenceCountedVector:
        size_t size()
        size_t index(size_t position)
        bool_cpp isSorted()
        double approximation(size_t position)
        const Rational& value(size_t position)
        bool_cpp operator==(const ReferenceCountedVector& other)
        bool_cpp operator!=(const ReferenceCountedVector& other)
        bool_cpp operator<(const ReferenceCountedVector& other)

    cdef struct InnerDescription:
        vector[Vector] points;
        vector[Vector] rays;

####################################
#IPO Linear Constraint

cdef extern from "ipo/linear_constraint.h" namespace "ipo":
    cdef cppclass LinearConstraint:
        LinearConstraint() except +
        LinearConstraint(char type, const Vector& normal, const Rational& rhs) except +

        bool_cpp operator==(const LinearConstraint& other)
        bool_cpp operator<(const LinearConstraint& other)
        bool_cpp isEquation()
        char type()
        const Vector& normal()
        const Rational& rhs()
        Rational getMaximumNorm()

        bool_cpp definesCompleteFace()
        bool_cpp definesEmptyFace()
        bool_cpp definesTrivialFace()

        int evaluatePoint(const Vector& point)
        int evaluateRay(const Vector& ray)

    ctypedef vector[LinearConstraint] AffineOuterDescription

####################################
#IPO Space

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream()

cdef extern from "ipo/space.h" namespace "ipo":
    cdef cppclass Space:
        Space() except +
        Space(const Space& other) except +

        bool_cpp operator==(const Space& other)
        bool_cpp operator!=(const Space& other)
        const string& operator[](size_t variable)

        size_t dimension()

        void printVector(ostream& stream, const Vector& vector)
        void printLinearForm(ostream& stream, const Vector& linearForm)
        void printLinearConstraint(ostream& stream, const LinearConstraint& constraint)

####################################
#IPO ScipOracle Python Wrapper

cdef extern from "ipo/python_wrapper.h" namespace "ipo":
    cdef cppclass ScipOracleController:
        ScipOracleController(string name) except +
        ScipOracleController(string name, ScipOracleController prev) except +

        string name();
        int heuristicLevel_ScipOracle();
        int heuristicLevel_CacheOracle();

        InnerDescription affineHullInner(int outputMode);
        AffineOuterDescription affineHullOuter(int outputMode);
