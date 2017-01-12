from libcpp.string cimport string

def bar():
  foo()

cdef foo():
  print('hello world!')

cdef extern from "ipo/scip_oracle.h" namespace "ipo":
  cdef cppclass SCIPOracle:
    SCIPOracle(string name) except +
    string name()

def test(fileName):
  ctest(fileName)

cdef ctest(fileName):
  print 'before'
  cdef SCIPOracle* oracle
  oracle = new SCIPOracle(fileName)
#  print oracle.name()
  del oracle
  print 'after'

#    cdef cppclass ScipOracleController(OracleControllerBase):
#        ScipOracleController(string name) except +
#        ScipOracleController(string name, OracleControllerBase* next) except +
#        ScipOracleController(string name, OracleControllerBase* next, const shared_ptr[MixedIntegerSet] mixedIntegerSet) except +
#        #ScipOracleController(string name, const shared_ptr[MixedIntegerSet] mixedIntegerSet) except +
#
#
#        string name()
#        int heuristicLevel_ScipOracle()
#        int heuristicLevel_CacheOracle()
#        shared_ptr[OracleBase] getCacheOracle()
#
#        InnerDescription affineHullInner(int outputMode)
#        AffineOuterDescription affineHullOuter(int outputMode)
