cimport cppExample
from libc cimport stdio
from cython.operator cimport dereference as deref, address as ref
from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_GE, Py_GT, Py_NE

cdef class Example:
    cdef cppExample.Foo* foo

    def __cinit__(self):
        self.foo = new cppExample.Foo()
        if self.foo is NULL:
            raise MemoryError()

    def __dealloc__(self):
        del self.foo

    def print_Example(self):
        return self.foo.printFoo()

cdef class Example2:
    cdef cppExample.Foo2* foo2
    def __cinit__(self):
        self.foo2 = new cppExample.Foo2()
        if self.foo2 is NULL:
            raise MemoryError()

    def __dealloc__(self):
        del self.foo2

    def print_Example(self):
        return self.foo2.printFoo2()

cdef class SoplexRational:
    cdef cppExample.Rational *cpp_rational
    cdef const cppExample.Rational *const_rational

    def __cinit__(self, isConst):
        if(not isConst):
            self.cpp_rational = new cppExample.Rational()
            if self.cpp_rational is NULL:
                raise MemoryError()

    def __dealloc__(self):
        del self.cpp_rational
        del self.const_rational
