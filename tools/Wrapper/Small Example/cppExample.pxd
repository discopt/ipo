from libcpp.string cimport string

cdef extern from "Foo.h":
    cdef cppclass Foo:
        Foo() except +
        string printFoo()

cdef extern from "Foo2.h":
    cdef cppclass Foo2:
        Foo2() except+
        string printFoo2()

cdef extern from "rational.h" namespace "soplex":
    cdef cppclass Rational:
        Rational() except +
