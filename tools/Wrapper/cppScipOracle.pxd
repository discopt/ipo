###################################
# ScipOracle.pxd, 1.0
# 14.10.2016, Sandra Hicks
#	Definition des Interface zu Matrix C++ Klasse
#
###################################
from libcpp.string cimport string
from cppIPOVector cimport InnerDescription
from cppIPOLinearConstraint cimport AffineOuterDescription

cdef extern from "scip_oracles.h":
    cdef cppclass ScipOracleController:
        ScipOracleController(string name) except +
        ScipOracleController(string name, ScipOracleController prev) except +

        string name();
        int heuristicLevel_ScipOracle();
        int heuristicLevel_CacheOracle();

        InnerDescription affineHullInner(int outputMode);
        AffineOuterDescription affineHullOuter(int outputMode);

## Wrapper ##
cdef object CreateScipOracle(string name, ScipOracleController *oracle)

cdef class ScipOracle:
    cdef ScipOracleController *oracle
