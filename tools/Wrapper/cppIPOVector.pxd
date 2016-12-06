###################################
# cppIPOVector.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu IPO Vector C++ Klasse
#
###################################
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp
from libcpp.vector cimport vector
from cppSoplexRational cimport Rational

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

## Wrapper ##
cdef object CreateIPOVector(Vector *vector)
cdef object CreateConstIPOVector(const Vector *vector)

cdef class IPOVector:
    cdef Vector *vec
    cdef const Vector *const_vec

cdef class IPOReferenceCountedVector:
    cdef ReferenceCountedVector *vec
    cdef const ReferenceCountedVector *const_vec
