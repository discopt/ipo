###################################
# cppLinearConstraint.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu LinearConstraint IPO C++ Klasse
#
###################################
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector
from cppIPOVector cimport Vector


cdef extern from "ipo/linear_constraint-pub.h" namespace "ipo":
    cdef cppclass LinearConstraint:
        LinearConstraint() except +
        #LinearConstraint(char type, const Vector& normal, const Rational& rhs) except +

        bool_cpp operator==(const LinearConstraint& other)
        bool_cpp isEquation()
        char type()
        Vector& normal()
        #Rational& rhs()
        #Rational getMaximumNorm()

        bool_cpp definesCompleteFace()
        bool_cpp definesEmptyFace()
        bool_cpp definesTrivialFace()

        int evaluatePoint(const Vector& point)
        int evaluateRay(const Vector& ray)

    ctypedef vector[LinearConstraint] AffineOuterDescription

## Wrapper ##
cdef object CreateLinearConstraint(LinearConstraint *linconst)

cdef class IPOLinearConstraint:
    cdef LinearConstraint *lin
