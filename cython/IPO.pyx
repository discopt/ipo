from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp cimport bool

## Rational ##

cdef extern from "ipo/rational.h" namespace "ipo":
  cdef cppclass IPORational "ipo::Rational":
    pass

  string IPOrationalToString "rationalToString" (const IPORational& rational)

cdef class Rational:
  cdef IPORational _rational

  def __str__(self):
    return IPOrationalToString(self._rational)

cdef _createRational(const IPORational& number):
  result = Rational()
  result._rational = number
  return result

## Space ##

cdef extern from "ipo/space.h" namespace "ipo":
  cdef cppclass IPOSpace "ipo::Space":

    const string& operator[](long index) const
    size_t dimension() const
    string vectorToString(const IPOVector& vector) const
    string linearConstraintToString(const IPOLinearConstraint& constraint) const

cdef class Space:
  cdef IPOSpace _space

  def __getitem__(self, index):
    assert index < self._space.dimension()
    return self._space[<long>(index)]

  @property
  def dimension(self):
    return self._space.dimension()

cdef _createSpace(const IPOSpace& space):
  result = Space()
  result._space = space
  return result

## Vector ##

cdef extern from "ipo/vectors-pub.h" namespace "ipo":
  cdef cppclass IPOReferenceCountedVector "ipo::Vector":
    size_t size() const
    size_t index(size_t position) const
    const IPORational& value(size_t position) const

  cdef cppclass IPOVector "ipo::Vector" (IPOReferenceCountedVector):
    pass

cdef class Vector:
  cdef IPOVector _vector
  cdef IPOSpace _space

  @property
  def space(self):
    return _createSpace(self._space)

  @property
  def size(self):
    return self._vector.size()

  def __str__(self):
    return self._space.vectorToString(self._vector)

  def index(self, position):
    assert position < self._vector.size()
    return self._vector.index(<long>(position))

  def value(self, position):
    assert position < self._vector.size()
    return _createRational(self._vector.value(<long>(position)))

cdef _createVector(const IPOVector& vector, const IPOSpace& space):
  result = Vector()
  result._vector = vector
  result._space = space
  return result

## LinearConstraint ##

cdef extern from "ipo/linear_constraint-pub.h" namespace "ipo":
  cdef cppclass IPOLinearConstraint "ipo::LinearConstraint":
    const IPOVector& normal() const

cdef class LinearConstraint:
  cdef IPOLinearConstraint _constraint
  cdef IPOSpace _space

  @property
  def space(self):
    return _createSpace(self._space)

  @property
  def normal(self):
    return _createVector(self._constraint.normal(), self._space)

  def __str__(self):
    return self._space.linearConstraintToString(self._constraint)

cdef _createLinearConstraint(const IPOLinearConstraint& constraint, const IPOSpace& space):
  result = LinearConstraint()
  result._constraint = constraint
  result._space = space
  return result

## InnerDescription ##

cdef extern from "ipo/vectors-pub.h" namespace "ipo":
  cdef cppclass IPOInnerDescription "ipo::InnerDescription":
    vector[IPOVector] points
    vector[IPOVector] rays

## Oracles ## 

cdef extern from "ipo/oracles.h" namespace "ipo":
  cdef cppclass IPOOracleBase "ipo::OracleBase":
    IPOOracleBase(const string& name) except +
    const IPOSpace& space() const
ctypedef shared_ptr[IPOOracleBase] PtrIPOOracleBase

cdef extern from "ipo/oracle_wrapper.h" namespace "ipo":
  cdef cppclass IPODefaultOracleWrapper "ipo::DefaultOracleWrapper":
    IPODefaultOracleWrapper(PtrIPOOracleBase& mainOracle) except +
    PtrIPOOracleBase queryOracle()
ctypedef shared_ptr[IPODefaultOracleWrapper] PtrIPODefaultOracleWrapper

cdef extern from "ipo/scip_oracle.h" namespace "ipo":
  cdef cppclass IPOSCIPOracle "ipo::SCIPOracle" (IPOOracleBase):
    IPOSCIPOracle(string name) except +
    string name()
    int heuristicLevel()
ctypedef shared_ptr[IPOSCIPOracle] PtrIPOSCIPOracle


cdef class SCIPOracle:
  cdef PtrIPOSCIPOracle _scipOracle
  cdef PtrIPODefaultOracleWrapper _wrappedOracle

  def __cinit__(self, fileName):
    self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(fileName))
    cdef PtrIPOOracleBase mainOracle = <PtrIPOOracleBase>(self._scipOracle)
    self._wrappedOracle = PtrIPODefaultOracleWrapper(new IPODefaultOracleWrapper(mainOracle))

  def __dealloc__(self):
    self._wrappedOracle.reset()
    self._scipOracle.reset()

  @property
  def name(self):
    return deref(self._scipOracle).name()

  @property
  def heuristicLevel(self):
    return deref(self._scipOracle).heuristicLevel()

  @property
  def space(self):
    cdef Space result = Space()
    result._space = deref(self._scipOracle).space()
    return result

## affineHull ##

ctypedef vector[IPOLinearConstraint] IPOAffineOuterDescription

cdef extern from "ipo/affine_hull.h" namespace "ipo":
  cdef cppclass IPOAffineHullHandler "ipo::AffineHullHandler":
    pass


  void IPOaffineHull "affineHull" (const PtrIPOOracleBase&, IPOInnerDescription&, IPOAffineOuterDescription&, vector[IPOAffineHullHandler*]&,
    size_t, size_t, vector[IPOLinearConstraint]&, bool)

cdef _affineHull(SCIPOracle oracle, lastCheapHeuristic, lastModerateHeuristic, givenEquations, approximateDirections):
  cdef IPOInnerDescription inner
  cdef IPOAffineOuterDescription outer
  cdef vector[IPOAffineHullHandler*] handlers
  cdef vector[IPOLinearConstraint] equations
  cdef bool approx = approximateDirections
  cdef PtrIPOOracleBase queryOracle = deref(oracle._wrappedOracle).queryOracle()
  IPOaffineHull(queryOracle, inner, outer, handlers, <size_t>(lastCheapHeuristic), <size_t>(lastModerateHeuristic), equations, approx)
  points = [None] * <long>(inner.points.size())
  for i in xrange(inner.points.size()):
    points[i] = _createVector(inner.points[i], deref(oracle._scipOracle).space())
  rays = [None] * <long>(inner.rays.size())
  for i in xrange(inner.rays.size()):
    rays[i] = _createVector(inner.rays[i], deref(oracle._scipOracle).space())
  return (points,rays)

def affineHull(oracle, lastCheapHeuristic, lastModerateHeuristic, givenEquations = [], approximateDirections = True):
  return _affineHull(oracle, lastModerateHeuristic, lastModerateHeuristic, givenEquations, approximateDirections)
  


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
