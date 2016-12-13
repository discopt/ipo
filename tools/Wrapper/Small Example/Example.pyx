cimport cppExample

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
