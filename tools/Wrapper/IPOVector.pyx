###################################
# IPOVector.pyx, 1.0
# 14.10.2016, Sandra Hicks
#	python wrapper Klasse, um die C++ Klasse Vector aus der IPO Lib anzusprechen
#
###################################
cimport cppIPOVector

cdef class IPOVector:
    cdef cppIPOVector.Vector* vec

    def __cinit__(self):
        self.vec = new cppIPOVector.Vector()
        if self.vec is NULL:
            raise MemoryError()
    def __dealloc__(self):
        print "dealloc C++ Vector"
        del self.vec
    def __call__(self):
        return deref(self._thisptr)()

    #def size(self):
        #(<ReferenceCountedVector *>(self.vec)).size()

