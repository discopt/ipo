###################################
# cppIPOSpace.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu LinearConstraint IPO C++ Klasse
#
###################################
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector
from cppIPOVector cimport Vector
from cppIPOLinearConstraint cimport LinearConstraint


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

## Wrapper ##
cdef object CreateSpace(Space *linconst)

cdef class IPOSpace:
    cdef Space *cpp_space

cdef extern from "<iostream>" namespace "std":
    cdef cppclass ostream:
        ostream()
