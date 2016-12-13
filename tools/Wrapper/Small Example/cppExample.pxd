from libcpp.string cimport string

cdef extern from "Foo.h":
    cdef cppclass Foo:
        Foo() except +

        string printFoo()
