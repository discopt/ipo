from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp cimport bool

## Rational ##

cdef extern from "ipo/rational.h" namespace "ipo":
  cdef cppclass IPORational "ipo::Rational":
    pass

  string IPOrationalToString "ipo::rationalToString" (const IPORational& rational)

cdef class Rational:
  cdef IPORational _rational

  def __repr__(self):
    return '{' + IPOrationalToString(self._rational) + '}'

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

  def __repr__(self):
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
    const IPORational& rhs() const
    bool isEquation() const
    char type() const

cdef class LinearConstraint:
  cdef IPOLinearConstraint _constraint
  cdef IPOSpace _space

  @property
  def space(self):
    return _createSpace(self._space)

  @property
  def normal(self):
    return _createVector(self._constraint.normal(), self._space)

  @property
  def rhs(self):
    return _createRational(self._constraint.rhs())

  @property
  def isEquation(self):
    return self._constraint.isEquation()

  @property
  def type(self):
    return self._constraint.type()

  def __repr__(self):
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

ctypedef size_t HeuristicLevel

cdef extern from "ipo/oracles.h" namespace "ipo":
  cdef cppclass IPOOracleBase "ipo::OracleBase":
    IPOOracleBase(const string& name) except +
    string name()
    const IPOSpace& space() const
    HeuristicLevel heuristicLevel() const
ctypedef shared_ptr[IPOOracleBase] PtrIPOOracleBase

cdef extern from "ipo/oracle_wrapper.h" namespace "ipo":
  cdef cppclass IPODefaultOracleWrapper "ipo::DefaultOracleWrapper":
    IPODefaultOracleWrapper(PtrIPOOracleBase& mainOracle) except +
    PtrIPOOracleBase queryOracle()
    PtrIPOOracleBase linkOracle()
ctypedef shared_ptr[IPODefaultOracleWrapper] PtrIPODefaultOracleWrapper

cdef class OracleBase:
  cdef PtrIPODefaultOracleWrapper _wrappedOracle

  @property
  def name(self):
    return deref(deref(self._wrappedOracle).linkOracle()).name()

  @property
  def heuristicLevel(self):
    return deref(deref(self._wrappedOracle).linkOracle()).heuristicLevel()

  @property
  def cacheHeuristicLevel(self):
    return deref(deref(self._wrappedOracle).queryOracle()).heuristicLevel()

  @property
  def space(self):
    cdef Space result = Space()
    result._space = deref(deref(self._wrappedOracle).linkOracle()).space()
    return result

  cdef PtrIPOOracleBase getLinkOraclePtr(self):
    return deref(self._wrappedOracle).linkOracle()

## MixedIntegerSet ##

cdef extern from "ipo/mip.h" namespace "ipo":
  cdef cppclass IPOMixedIntegerSetVariable "ipo::MixedIntegerSet::Variable":
    bool integral
    IPORational upperBound
    IPORational lowerBound

  cdef cppclass IPOMixedIntegerSet "ipo::MixedIntegerSet":
    IPOMixedIntegerSet(const string& fileName)
    const IPOSpace& space() const
    size_t numVariables() const
    size_t numRows() const
    const IPOMixedIntegerSetVariable& variable(size_t) const
    IPOLinearConstraint upperBoundConstraint(size_t) const
    IPOLinearConstraint lowerBoundConstraint(size_t) const
    const IPOLinearConstraint& rowConstraint(size_t row) const
    const vector[IPOLinearConstraint]& rowConstraints() const
    const string& rowName(size_t) const
ctypedef shared_ptr[IPOMixedIntegerSet] PtrIPOMixedIntegerSet

class Variable:
  pass

cdef class MixedIntegerSet:
  cdef PtrIPOMixedIntegerSet _mixedIntegerSet

  def __cinit__(self, const string& fileName):
    self._mixedIntegerSet = PtrIPOMixedIntegerSet(new IPOMixedIntegerSet(fileName))

  @property
  def space(self):
    return _createSpace(deref(self._mixedIntegerSet).space())

  @property
  def numVariables(self):
    return deref(self._mixedIntegerSet).numVariables()

  @property
  def numRows(self):
    return deref(self._mixedIntegerSet).numRows()

  def variable(self, long variableIndex):
    cdef IPOMixedIntegerSetVariable var = deref(self._mixedIntegerSet).variable(variableIndex)
    result = Variable()
    result.integral = var.integral
    result.upperBound = _createRational(var.upperBound)
    result.lowerBound = _createRational(var.lowerBound)
    return result

  @property
  def variables(self):
    result = [None] * self.numVariables
    for i in xrange(self.numVariables):
      result[i] = self.variable(i)
    return result

  def upperBoundConstraint(self, long variableIndex):
    return _createLinearConstraint(deref(self._mixedIntegerSet).upperBoundConstraint(variableIndex), deref(self._mixedIntegerSet).space())

  def lowerBoundConstraint(self, long variableIndex):
    return _createLinearConstraint(deref(self._mixedIntegerSet).lowerBoundConstraint(variableIndex), deref(self._mixedIntegerSet).space())

  def rowConstraint(self, long rowIndex):
    return _createLinearConstraint(deref(self._mixedIntegerSet).rowConstraint(rowIndex), deref(self._mixedIntegerSet).space())

  @property
  def rowConstraints(self):
    result = [None] * self.numRows
    for i in xrange(self.numRows):
      result[i] = self.rowConstraint(i)
    return result
  
  def rowName(self, long rowIndex):
    return deref(self._mixedIntegerSet).rowName(rowIndex) 

cdef _createMixedIntegerSet(const PtrIPOMixedIntegerSet& mixedIntegerSet):
  result = MixedIntegerSet()
  result._mixedIntegerSet = mixedIntegerSet
  return result

## SCIPOracle ##

cdef extern from "ipo/scip_oracle.h" namespace "ipo":
  cdef cppclass IPOSCIPOracle "ipo::SCIPOracle" (IPOOracleBase):
    IPOSCIPOracle(const string&) except +
    IPOSCIPOracle(const string&, const PtrIPOOracleBase&) except +
    IPOSCIPOracle(const string&, const PtrIPOMixedIntegerSet& mixedIntegerSet) except +
    IPOSCIPOracle(const string&, const PtrIPOMixedIntegerSet&, const PtrIPOOracleBase&) except +
    const PtrIPOMixedIntegerSet& mixedIntegerSet() const
    double setTimeLimit(double)
    double getTimeLimit()
ctypedef shared_ptr[IPOSCIPOracle] PtrIPOSCIPOracle

cdef class SCIPOracle (OracleBase):
  cdef PtrIPOSCIPOracle _scipOracle

  def __cinit__(self, const string& name, MixedIntegerSet mixedIntegerSet = None, OracleBase nextOracle = None):
    if mixedIntegerSet is None:
      if nextOracle is None:
        self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name))
      else:
        self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name, nextOracle.getLinkOraclePtr()))
    else:
      if nextOracle is None:
        self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name, mixedIntegerSet._mixedIntegerSet))
      else:
        self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name, mixedIntegerSet._mixedIntegerSet, nextOracle.getLinkOraclePtr()))
    cdef PtrIPOOracleBase mainOracle = <PtrIPOOracleBase>(self._scipOracle)
    self._wrappedOracle = PtrIPODefaultOracleWrapper(new IPODefaultOracleWrapper(mainOracle))

  @property
  def mixedIntegerSet(self):
    return _createMixedIntegerSet(deref(self._scipOracle).mixedIntegerSet())

  @property
  def timeLimit(self):
    return deref(self._scipOracle).getTimeLimit()

  @timeLimit.setter
  def timeLimit(self, limit):
    assert limit >= 0
    deref(self._scipOracle).setTimeLimit(limit)

## ExactSCIPOracle ##

cdef extern from "ipo/exactscip_oracle.h" namespace "ipo":
  cdef cppclass IPOExactSCIPOracle "ipo::ExactSCIPOracle" (IPOOracleBase):
    IPOExactSCIPOracle(const string& name, const PtrIPOMixedIntegerSet& mixedIntegerSet) except +
    const PtrIPOMixedIntegerSet& mixedIntegerSet() const
    double setTimeLimit(double)
    double getTimeLimit()
ctypedef shared_ptr[IPOExactSCIPOracle] PtrIPOExactSCIPOracle

cdef class ExactSCIPOracle (OracleBase):
  cdef PtrIPOExactSCIPOracle _exactscipOracle

  def _init(self, name, MixedIntegerSet mixedIntegerSet):
    cdef PtrIPOMixedIntegerSet mis = mixedIntegerSet._mixedIntegerSet
    self._exactscipOracle = PtrIPOExactSCIPOracle(new IPOExactSCIPOracle(name, mis))
    cdef PtrIPOOracleBase mainOracle = <PtrIPOOracleBase>(self._exactscipOracle)
    self._wrappedOracle = PtrIPODefaultOracleWrapper(new IPODefaultOracleWrapper(mainOracle))

  def __cinit__(self, name, mixedIntegerSet):
    self._init(name, mixedIntegerSet)

  @property
  def mixedIntegerSet(self):
    return _createMixedIntegerSet(deref(self._exactscipOracle).mixedIntegerSet())

  @property
  def timeLimit(self):
    return deref(self._exactscipOracle).getTimeLimit()

  @timeLimit.setter
  def timeLimit(self, limit):
    assert limit >= 0
    deref(self._exactscipOracle).setTimeLimit(limit)

## affineHull ##

ctypedef vector[IPOLinearConstraint] IPOAffineOuterDescription

cdef extern from "ipo/affine_hull.h" namespace "ipo":
  cdef cppclass IPOAffineHullHandler "ipo::AffineHullHandler":
    pass


  void IPOaffineHull "affineHull" (const PtrIPOOracleBase&, IPOInnerDescription&, IPOAffineOuterDescription&, vector[IPOAffineHullHandler*]&,
    size_t, size_t, vector[IPOLinearConstraint]&, bool)

def affineHull(OracleBase oracle, HeuristicLevel lastCheapHeuristic = -1, HeuristicLevel lastModerateHeuristic = -1,
  givenEquations = [], approximateDirections = True):
  # Setup parameters.
  if lastCheapHeuristic < 0:
    lastCheapHeuristic = oracle.heuristicLevel + 1
  if lastModerateHeuristic < 0:
    lastModerateHeuristic = oracle.heuristicLevel
  cdef PtrIPOOracleBase queryOracle = deref(oracle._wrappedOracle).queryOracle()
  cdef IPOInnerDescription inner
  cdef IPOAffineOuterDescription outer
  cdef vector[IPOAffineHullHandler*] handlers
  cdef vector[IPOLinearConstraint] givenEqns
  cdef bool approx = approximateDirections

  # Call routine.

  IPOaffineHull(queryOracle, inner, outer, handlers, <size_t>(lastCheapHeuristic), <size_t>(lastModerateHeuristic),
    givenEqns, approx)

  # Compute dimension.
  dim = <long>(inner.points.size()) + <long>(inner.rays.size()) - 1

  # Extract inner description.
  points = [None] * <long>(inner.points.size())
  for i in xrange(inner.points.size()):
    points[i] = _createVector(inner.points[i], deref(queryOracle).space())
  rays = [None] * <long>(inner.rays.size())
  for i in xrange(inner.rays.size()):
    rays[i] = _createVector(inner.rays[i], deref(queryOracle).space())

  # Extract outer description
  equations = [None] * <long>(outer.size())
  for i in xrange(outer.size()):
    equations[i] = _createLinearConstraint(outer[i], deref(queryOracle).space())

  return (dim, (points, rays), equations)

## Polyhedron::Face ##

cdef extern from "ipo/polyhedron.h" namespace "ipo":
  cdef cppclass IPOPolyhedronFace "ipo::Polyhedron::Face":
    IPOPolyhedronFace(const IPOLinearConstraint&)
    const IPOLinearConstraint& inequality() const
    bool hasDimension() const
    const IPOAffineOuterDescription& outerDescription() const
    const IPOInnerDescription& innerDescription() const
    const int dimension() const
ctypedef shared_ptr[IPOPolyhedronFace] PtrIPOPolyhedronFace
  
cdef class Face:
  cdef PtrIPOPolyhedronFace _face
  cdef IPOSpace _space

  def __cinit__(self):
    pass

  @property
  def inequality(self):
    return _createLinearConstraint(deref(self._face).inequality(), self._space)

  @property
  def hasDimension(self):
    return deref(self._face).hasDimension()

  @property
  def affineHullOuterDescription(self):
    cdef IPOAffineOuterDescription outer = deref(self._face).outerDescription()
    equations = [None] * <long>(outer.size())
    for i in xrange(outer.size()):
      equations[i] = _createLinearConstraint(outer[i], self._space)
    return equations

  @property
  def affineHullInnerDescription(self):
    cdef IPOInnerDescription inner = deref(self._face).innerDescription()
    points = [None] * <long>(inner.points.size())
    for i in xrange(inner.points.size()):
      points[i] = _createVector(inner.points[i], self._space)
    rays = [None] * <long>(inner.rays.size())
    for i in xrange(inner.rays.size()):
      rays[i] = _createVector(inner.rays[i], self._space)
    return (points, rays)

  @property
  def dimension(self):
    return deref(self._face).dimension()

cdef _createFace(const PtrIPOPolyhedronFace& face, const IPOSpace& space):
  result = Face()
  result._face = face
  result._space = space
  return result

cdef extern from "ipo/polyhedron.h" namespace "ipo":
  cdef cppclass IPOPolyhedron "ipo::Polyhedron":
    IPOPolyhedron(PtrIPODefaultOracleWrapper&) except +
    IPOSpace space() const
    size_t numPoints() const
    size_t numRays() const
    size_t numInequalities() const
    int dimension()
    void affineHull(PtrIPOPolyhedronFace& face)
    const IPOAffineOuterDescription& affineHullOuterDescription()
    const IPOInnerDescription& affineHullInnerDescription()
    void setAffineHullLastCheapHeuristicLevel(HeuristicLevel)
    HeuristicLevel getAffineHullLastCheapHeuristicLevel() const
    void setAffineHullLastModerateHeuristicLevel(HeuristicLevel)
    HeuristicLevel getAffineHullLastModerateHeuristicLevel() const
    PtrIPOPolyhedronFace inequalityToFace(const IPOLinearConstraint&)
    void addConstraint(const IPOLinearConstraint&, bool);
    void getFaces(vector[PtrIPOPolyhedronFace]& constraints, bool, bool);
ctypedef shared_ptr[IPOPolyhedron] PtrIPOPolyhedron

## Polyhedron ##

cdef class Polyhedron:
  cdef PtrIPOPolyhedron _polyhedron

  def __cinit__(self, OracleBase oracle):
    self._polyhedron = PtrIPOPolyhedron(new IPOPolyhedron(oracle._wrappedOracle))

  @property
  def space(self):
    cdef Space result = Space()
    result._space = deref(self._polyhedron).space()
    return result

  @property
  def numPoints(self):
    return deref(self._polyhedron).numPoints()

  @property
  def numRays(self):
    return deref(self._polyhedron).numRays()

  @property
  def numInequalities(self):
    return deref(self._polyhedron).numInequalities()

  def inequalityToFace(self, LinearConstraint constraint):
    return _createFace(deref(self._polyhedron).inequalityToFace(constraint._constraint), deref(self._polyhedron).space())

  @property
  def dimension(self):
    return deref(self._polyhedron).dimension()

  def affineHull(self, Face face):
    deref(self._polyhedron).affineHull(face._face)

  @property
  def affineHullOuterDescription(self):
    cdef IPOAffineOuterDescription outer = deref(self._polyhedron).affineHullOuterDescription()
    equations = [None] * <long>(outer.size())
    for i in xrange(outer.size()):
      equations[i] = _createLinearConstraint(outer[i], deref(self._polyhedron).space())
    return equations

  @property
  def affineHullInnerDescription(self):
    cdef IPOInnerDescription inner = deref(self._polyhedron).affineHullInnerDescription()
    points = [None] * <long>(inner.points.size())
    for i in xrange(inner.points.size()):
      points[i] = _createVector(inner.points[i], deref(self._polyhedron).space())
    rays = [None] * <long>(inner.rays.size())
    for i in xrange(inner.rays.size()):
      rays[i] = _createVector(inner.rays[i], deref(self._polyhedron).space())
    return (points, rays)

  @property
  def affineHullLastCheapHeuristicLevel(self):
    return deref(self._polyhedron).getAffineHullLastCheapHeuristicLevel()

  @affineHullLastCheapHeuristicLevel.setter
  def affineHullLastCheapHeuristicLevel(self, HeuristicLevel heuristicLevel):
    deref(self._polyhedron).setAffineHullLastCheapHeuristicLevel(heuristicLevel)

  @property
  def affineHullLastModerateHeuristicLevel(self):
    return deref(self._polyhedron).getAffineHullLastModerateHeuristicLevel()

  @affineHullLastModerateHeuristicLevel.setter
  def affineHullLastModerateHeuristicLevel(self, HeuristicLevel heuristicLevel):
    deref(self._polyhedron).setAffineHullLastModerateHeuristicLevel(heuristicLevel)
    
  def addConstraint(self, LinearConstraint constraint, bool normalize = True):
    deref(self._polyhedron).addConstraint(constraint._constraint, normalize)

  def getFaces(self, bool onlyInequalities = True, bool onlyWithDimension = False):
    cdef vector[PtrIPOPolyhedronFace] faces
    deref(self._polyhedron).getFaces(faces, onlyInequalities, onlyWithDimension)
    result = [None] * <long>(faces.size())
    for i in xrange(faces.size()):
      result[i] = _createFace(faces[i], deref(self._polyhedron).space())
    return result

## SoPlex ##


