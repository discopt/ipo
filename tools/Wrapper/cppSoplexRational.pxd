###################################
# cppIPOLinearConstraint.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu Rational Soplex C++ Klasse
#
###################################
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector

cdef extern from "rational.h" namespace "soplex":
    cdef cppclass Rational:
        Rational() except +

cdef object CreateSoplexRational(Rational *rational)
cdef object CreateConstSoplexRational(const Rational *rational)

cdef class SoplexRational:
    cdef Rational *cpp_rational
