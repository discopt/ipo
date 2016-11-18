###################################
# cppIPOVector.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu IPO Vector C++ Klasse
#
###################################
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector

cdef extern from "ipo/vectors-pub.h" namespace "ipo":
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
        bool_cpp operator==(const ReferenceCountedVector& other)
        bool_cpp operator!=(const ReferenceCountedVector& other)
        bool_cpp operator<(const ReferenceCountedVector& other)

    cdef struct InnerDescription:
        vector[Vector] points;
        vector[Vector] rays;

## Wrapper ##
cdef object CreateIPOVector(Vector *vector)

cdef class IPOVector:
    cdef Vector *vec
