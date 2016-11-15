###################################
# ScipOracle.pyx, 1.0
# 14.10.2016, Sandra Hicks
#	python wrapper Klasse, um den ScipOracleController-Wrapper in C++ anzusprechen
# besteht aus einem ScipOracle, CacheOracle und StatisticsOracle
#
###################################
cimport cppScipOracle
cimport cppIPOVector
cimport cppIPOLinearConstraint
import IPOVector
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

cdef object CreateScipOracle(str name, cppScipOracle.ScipOracleController oracle):
    py_oracle = IPOScipOracle(name, 0)
    py_oracle.oracle = new cppScipOracle.ScipOracleController(name, oracle)
    return py_oracle

cdef class IPOScipOracle:
    cdef cppScipOracle.ScipOracleController *oracle

    def __cinit__(self, str name, int isNew):
        if(isNew == 1):
            self.oracle = new cppScipOracle.ScipOracleController(name)

    def __dealloc__(self):
        del self.oracle

    def name(self):
        return self.oracle.name()

    def heuristicLevel(self):
        return self.oracle.heuristicLevel_ScipOracle()

    def heuristicLevel_CacheOracle(self):
        return self.oracle.heuristicLevel_CacheOracle()

    def affineHull(self, outputMode):
        #########-1-#########
        #convert inner description to 2-tupel IPOVector lists
        cdef cppIPOVector.InnerDescription c_inner = self.oracle.affineHullInner(1)

        #automagically convert to python list
        c_points = c_inner.points
        c_rays = c_inner.rays
        points = ()
        rays = ()
        cdef cppIPOVector.Vector *c_vector
        #convert points to python wrapperclass IPOVector
        for i in range(0,c_points.size()):
            c_vector = <cppIPOVector.Vector *>c_points[i]
            py_vector = cppIPOVector.CreateIPOVector(c_vector)
            points.append(py_vector)

        #convert rays to python wrapperclass IPOVector
        for i in range(0,c_rays.size()):
            c_vector = <cppIPOVector.Vector *>c_rays[i]
            py_vector = cppIPOVector.CreateIPOVector(c_vector)
            rays.append(py_vector)

        inner = (points, rays)
        #########-2-#########
        #convert outer description to IPOLinearConstraint list
        cdef cppIPOLinearConstraint.AffineOuterDescription c_outer = self.oracle.affineHullOuter(1)
        outer = ()

        for i in range(0, c_outer.size()):
            c_linconst = <cppIPOLinearConstraint.LinearConstraint *>c_outer[i]
            py_linconst = cppIPOLinearConstraint.CreateLinearConstraint(c_linconst)
            outer.append(py_linconst)

        return (inner, outer)
