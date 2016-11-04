###################################
# Vector.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu Matrix C++ Klasse
#
###################################
from libcpp.string cimport string
from libcpp cimport bool as bool_cpp

cdef extern from "ipo/vectors-pub.h":
    cdef cppclass MutableVector:
        MutableVector() except +
    cdef cppclass Vector:

        Vector() except +
        void swap(Vector& other)
        Vector& operator=(const Vector& other)
        Vector& operator=(const MutableVector& other)


    cdef cppclass ReferenceCountedVector:
        size_t size()
        size_t index(size_t position)
        bool_cpp isSorted()
        double approximation(size_t position)
        bool_cpp operator==(const ReferenceCountedVector& other)
        bool_cpp operator!=(const ReferenceCountedVector& other)
        bool_cpp operator<(const ReferenceCountedVector& other)
