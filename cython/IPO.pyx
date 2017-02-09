from cython.operator cimport dereference as deref
from libcpp.string cimport string
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector
from libcpp cimport bool

include "config.pxi" # Compile-time definitions

## Rational ##

cdef extern from "ipo/rational.h" namespace "ipo":
  cdef cppclass IPORational "ipo::Rational":
    IPORational()
    IPORational(double x)
    IPORational(int x)

  string IPOrationalToString "ipo::rationalToString" (const IPORational& rational)

cdef class Rational:
  cdef IPORational _rational

  def __init__(self, init = None):
    if isinstance(init, float):
      self._initFromFloat(init)
    elif isinstance(init, int):
      self._initFromInt(init)

  def _initFromFloat(self, double x):
    self._rational = IPORational(x)

  def _initFromInt(self, int x):
    self._rational = IPORational(x)

  def __repr__(self):
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

  def __repr__(self):
    return '(' + ','.join([ self.__getitem__(v) for v in xrange(self.dimension) ]) + ')'

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

cdef extern from "ipo/linear_constraint.h" namespace "ipo":
  void IPOscaleIntegral "ipo::scaleIntegral" (IPOLinearConstraint& vector)

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

## LinearSet ##

cdef extern from "ipo/lp.h" namespace "ipo":
  cdef cppclass IPOLinearSet "ipo::LinearSet":
    IPOLinearSet(const string& fileName)
    IPOLinearSet(const IPOLinearSet& linearSet)
    const IPOSpace& space() const
    size_t numVariables() const
    size_t numRows() const
    const IPORational& lowerBound(size_t) const
    IPOLinearConstraint lowerBoundConstraint(size_t) const
    const IPORational& upperBound(size_t) const
    IPOLinearConstraint upperBoundConstraint(size_t) const
    const string& variableName(size_t) const
    const string& rowName(size_t) const
    const IPOLinearConstraint& rowConstraint(size_t) const
    void getConstraints(vector[IPOLinearConstraint]&, bool, bool, bool) const
    void changeLowerBound(size_t, const IPORational&)
    void changeUpperBound(size_t, const IPORational&)
    void changeBounds(size_t, const IPORational&, const IPORational&)
    void addConstraint(const IPOLinearConstraint&, const string&)
    void removeConstraint(size_t)
ctypedef shared_ptr[IPOLinearSet] PtrIPOLinearSet

cdef class LinearSet:
  cdef PtrIPOLinearSet _linearSet

  def __init__(self, init):
    if isinstance(init, str):
      self._initFromString(init)
    elif isinstance(init, LinearSet):
      self._initFromLinearSet(init)
    else:
      raise Exception('Cannot initialize IPO.LinearSet from %s' % (type(init)))

  def _initFromLinearSet(self, LinearSet linearSet):
    self._linearSet = PtrIPOLinearSet(new IPOLinearSet(deref(linearSet._getIPOLinearSet())))

  def _initFromString(self, const string& fileName):
    self._linearSet = PtrIPOLinearSet(new IPOLinearSet(fileName))

  cdef IPOLinearSet* _getIPOLinearSet(self):
    return self._linearSet.get()

  cdef IPOSpace _getIPOSpace(self):
    return deref(self._getIPOLinearSet()).space()

  @property
  def space(self):
    return _createSpace(self._getIPOLinearSet().space())

  @property
  def numVariables(self):
    return deref(self._getIPOLinearSet()).numVariables()

  @property
  def numRows(self):
    return deref(self._getIPOLinearSet()).numRows()

  def lowerBound(self, long variableIndex):
    return _createRational(deref(self._getIPOLinearSet()).lowerBound(variableIndex))

  def lowerBoundConstraint(self, long variableIndex):
    return _createLinearConstraint(deref(self._getIPOLinearSet()).lowerBoundConstraint(variableIndex), self._getIPOSpace())

  def upperBound(self, long variableIndex):
    return _createRational(deref(self._getIPOLinearSet()).upperBound(variableIndex))

  def upperBoundConstraint(self, long variableIndex):
    return _createLinearConstraint(deref(self._getIPOLinearSet()).upperBoundConstraint(variableIndex), self._getIPOSpace())

  def variableName(self, long variableIndex):
    return deref(self._getIPOLinearSet()).variableName(variableIndex)

  def rowName(self, long rowIndex):
    return deref(self._getIPOLinearSet()).rowName(rowIndex)

  def rowConstraint(self, long rowIndex):
    return _createLinearConstraint(deref(self._getIPOLinearSet()).rowConstraint(rowIndex), self._getIPOSpace())

  def getConstraints(self, bool excludeEquations = False, bool includeBounds = True, bool includeRows = True):
    cdef vector[IPOLinearConstraint] cons
    deref(self._getIPOLinearSet()).getConstraints(cons, excludeEquations, includeBounds, includeRows)
    result = [None] * <long>(cons.size())
    for i in xrange(len(result)):
      result[i] = _createLinearConstraint(cons[i], self._getIPOSpace())
    return result

  def changeLowerBound(self, long variableIndex, Rational newLower):
    deref(self._getIPOLinearSet()).changeLowerBound(variableIndex, newLower._rational)

  def changeUpperBound(self, long variableIndex, Rational newUpper):
    deref(self._getIPOLinearSet()).changeUpperBound(variableIndex, newUpper._rational)

  def changeBounds(self, long variableIndex, Rational newLower, Rational newUpper):
    deref(self._getIPOLinearSet()).changeBounds(variableIndex, newLower._rational, newUpper._rational)

  def addConstraint(self, LinearConstraint constraint, const string& name = ''):
    deref(self._getIPOLinearSet()).addConstraint(constraint._constraint, name)

  def removeConstraint(self, long rowIndex):
    deref(self._getIPOLinearSet()).removeConstraint(rowIndex)

cdef _createLinearSet(const PtrIPOLinearSet& linearSet):
  result = LinearSet()
  result._linearSet = linearSet
  return result

## MixedIntegerLinearSet ##

cdef extern from "ipo/lp.h" namespace "ipo":
  cdef cppclass IPOMixedIntegerLinearSet "ipo::MixedIntegerLinearSet" (IPOLinearSet):
    IPOMixedIntegerLinearSet(const string& fileName)
    IPOMixedIntegerLinearSet(const IPOMixedIntegerLinearSet& mixedIntegerLinearSet)
    IPOMixedIntegerLinearSet(const IPOLinearSet& linearSet)
    bool isIntegral(size_t) const
ctypedef shared_ptr[IPOMixedIntegerLinearSet] PtrIPOMixedIntegerLinearSet

cdef class MixedIntegerLinearSet (LinearSet):
  cdef PtrIPOMixedIntegerLinearSet _mixedIntegerLinearSet

  def __init__(self, init):
    if isinstance(init, str):
      self._initFromString(init)
    elif isinstance(init, LinearSet):
      self._initFromLinearSet(init)
    else:
      raise Exception('Cannot initialize IPO.MixedIntegerLinearSet from %s' % (type(init)))

  def _initFromLinearSet(self, LinearSet linearSet):
    self._mixedIntegerLinearSet = PtrIPOMixedIntegerLinearSet(new IPOMixedIntegerLinearSet(deref(linearSet._getIPOLinearSet())))

  def _initFromString(self, const string& fileName):
    self._mixedIntegerLinearSet = PtrIPOMixedIntegerLinearSet(new IPOMixedIntegerLinearSet(fileName))

  cdef IPOLinearSet* _getIPOLinearSet(self):
    return self._mixedIntegerLinearSet.get()

  def isIntegral(self, long variableIndex):
    return deref(self._mixedIntegerLinearSet).isIntegral(variableIndex)

cdef _createMixedIntegerLinearSet(const PtrIPOMixedIntegerLinearSet& mixedIntegerLinearSet):
  result = MixedIntegerLinearSet()
  result._mixedIntegerLinearSet = mixedIntegerLinearSet
  return result

## LinearProgram ##

cdef extern from "ipo/lp.h" namespace "ipo":
  cdef enum IPOLinearProgramResult "ipo::LinearProgram::Result":
    INFEASIBLE,
    OPTIMAL,
    UNBOUNDED

cdef extern from "ipo/lp.h" namespace "ipo":
  cdef cppclass IPOLinearProgram "ipo::LinearProgram" (IPOLinearSet):
    IPOLinearProgram(const string& fileName)
    IPOLinearProgram(const PtrIPOLinearSet& linearSet)
    IPOLinearProgram(const IPOLinearSet& other)
    void changeObjective(const vector[IPORational]&)
    IPOLinearProgramResult solve(IPOVector&, IPORational&)
ctypedef shared_ptr[IPOLinearProgram] PtrIPOLinearProgram

cdef class LinearProgram (LinearSet):
  cdef PtrIPOLinearProgram _linearProgram

  INFEASIBLE = -1
  OPTIMAL = 0
  UNBOUNDED = 1

  def __init__(self, init):
    if isinstance(init, str):
      self._initFromString(init)
    if isinstance(init, LinearProgram):
      self._initFromLinearProgram(init)
    elif isinstance(init, LinearSet):
      self._initFromLinearSet(init)
    else:
      raise Exception('Cannot initialize IPO.LinearProgram from %s' % (type(init)))

  def _initFromLinearProgram(self, LinearProgram linearProgram):
    self._linearProgram = PtrIPOLinearProgram(new IPOLinearProgram(deref(linearProgram._getIPOLinearProgram())))

  def _initFromLinearSet(self, LinearSet linearSet):
    self._linearProgram = PtrIPOLinearProgram(new IPOLinearProgram(deref(linearSet._getIPOLinearSet())))

  def _initFromString(self, const string& fileName):
    self._linearProgram = PtrIPOLinearProgram(new IPOLinearProgram(fileName))

  cdef IPOLinearProgram* _getIPOLinearProgram(self):
    return self._linearProgram.get()

  cdef IPOLinearSet* _getIPOLinearSet(self):
    return self._linearProgram.get()

  def changeObjective(self, objective):
    assert len(objective) == self.numVariables
    cdef Rational x
    cdef vector[IPORational] obj
    obj.resize(self.numVariables)
    for i in xrange(self.numVariables):
      c = objective[i]
      if isinstance(c, Rational):
        x = objective[i]
      else:
        x = Rational(c)
      obj[i] = x._rational
    deref(self._linearProgram).changeObjective(obj)

  def solve(self):
    cdef IPOVector vector
    cdef IPORational objectiveValue
    cdef IPOLinearProgramResult result
    result = deref(self._linearProgram).solve(vector, objectiveValue)
    return (result, _createVector(vector, self._getIPOSpace()), _createRational(objectiveValue))

cdef _createLinearProgram(const PtrIPOLinearProgram& linearProgram):
  result = LinearProgram()
  result._linearProgram = linearProgram
  return result

IF IPO_WITH_SCIP:

  ## SCIPOracle ##

  cdef extern from "ipo/scip_oracle.h" namespace "ipo":
    cdef cppclass IPOSCIPOracle "ipo::SCIPOracle" (IPOOracleBase):
      IPOSCIPOracle(const string&) except +
      IPOSCIPOracle(const string&, const PtrIPOOracleBase&) except +
      IPOSCIPOracle(const string&, const PtrIPOMixedIntegerLinearSet&) except +
      IPOSCIPOracle(const string&, const PtrIPOMixedIntegerLinearSet&, const PtrIPOOracleBase&) except +
      const PtrIPOMixedIntegerLinearSet& mixedIntegerLinearSet() const
      double setTimeLimit(double)
      double getTimeLimit()
  ctypedef shared_ptr[IPOSCIPOracle] PtrIPOSCIPOracle

  cdef class SCIPOracle (OracleBase):
    cdef PtrIPOSCIPOracle _scipOracle

    def __cinit__(self, const string& name, MixedIntegerLinearSet mixedIntegerLinearSet = None, OracleBase nextOracle = None):
      if mixedIntegerLinearSet is None:
        if nextOracle is None:
          self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name))
        else:
          self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name, nextOracle.getLinkOraclePtr()))
      else:
        if nextOracle is None:
          self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name, mixedIntegerLinearSet._mixedIntegerLinearSet))
        else:
          self._scipOracle = PtrIPOSCIPOracle(new IPOSCIPOracle(name, mixedIntegerLinearSet._mixedIntegerLinearSet, nextOracle.getLinkOraclePtr()))
      cdef PtrIPOOracleBase mainOracle = <PtrIPOOracleBase>(self._scipOracle)
      self._wrappedOracle = PtrIPODefaultOracleWrapper(new IPODefaultOracleWrapper(mainOracle))

    @property
    def mixedIntegerLinearSet(self):
      return _createMixedIntegerLinearSet(deref(self._scipOracle).mixedIntegerLinearSet())

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

    IPOExactSCIPOracle(const string& name, const PtrIPOMixedIntegerLinearSet&) except +
    IPOExactSCIPOracle(const string& binary, const string& name, const PtrIPOMixedIntegerLinearSet&) except +
    const PtrIPOMixedIntegerLinearSet& mixedIntegerLinearSet() const
    double setTimeLimit(double)
    double getTimeLimit()
ctypedef shared_ptr[IPOExactSCIPOracle] PtrIPOExactSCIPOracle

cdef class ExactSCIPOracle (OracleBase):
  cdef PtrIPOExactSCIPOracle _exactscipOracle

  IF IPO_WITH_EXACT_SCIP:

    def _init(self, name, MixedIntegerLinearSet mixedIntegerLinearSet):
      cdef PtrIPOMixedIntegerLinearSet mis = mixedIntegerLinearSet._mixedIntegerLinearSet
      self._exactscipOracle = PtrIPOExactSCIPOracle(new IPOExactSCIPOracle(name, mis))
      cdef PtrIPOOracleBase mainOracle = <PtrIPOOracleBase>(self._exactscipOracle)
      self._wrappedOracle = PtrIPODefaultOracleWrapper(new IPODefaultOracleWrapper(mainOracle))

    def __cinit__(self, name, mixedIntegerLinearSet):
      self._init(name, mixedIntegerLinearSet)

  ELSE:

    def _init(self, binary, name, MixedIntegerLinearSet mixedIntegerLinearSet):
      cdef PtrIPOMixedIntegerLinearSet mis = mixedIntegerLinearSet._mixedIntegerLinearSet
      self._exactscipOracle = PtrIPOExactSCIPOracle(new IPOExactSCIPOracle(binary, name, mis))
      cdef PtrIPOOracleBase mainOracle = <PtrIPOOracleBase>(self._exactscipOracle)
      self._wrappedOracle = PtrIPODefaultOracleWrapper(new IPODefaultOracleWrapper(mainOracle))

    def __cinit__(self, binary, name, mixedIntegerLinearSet):
      self._init(binary, name, mixedIntegerLinearSet)

  @property
  def mixedIntegerLinearSet(self):
    return _createMixedIntegerLinearSet(deref(self._exactscipOracle).mixedIntegerLinearSet())

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
    size_t numConstraints() const
    int dimension()
    void affineHull(PtrIPOPolyhedronFace& face)
    const IPOAffineOuterDescription& affineHullOuterDescription()
    const IPOInnerDescription& affineHullInnerDescription()
    void setAffineHullLastCheapHeuristicLevel(HeuristicLevel)
    HeuristicLevel getAffineHullLastCheapHeuristicLevel() const
    void setAffineHullLastModerateHeuristicLevel(HeuristicLevel)
    HeuristicLevel getAffineHullLastModerateHeuristicLevel() const
    PtrIPOPolyhedronFace constraintToFace(const IPOLinearConstraint&)
    void addConstraint(const IPOLinearConstraint&, bool);
    void getFaces(vector[PtrIPOPolyhedronFace]& constraints, bool, bool);
    bool separatePoint(const IPOVector&, IPOLinearConstraint&, IPOInnerDescription*)
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
  def numConstraints(self):
    return deref(self._polyhedron).numConstraints()

  def constraintToFace(self, LinearConstraint constraint):
    return _createFace(deref(self._polyhedron).constraintToFace(constraint._constraint), deref(self._polyhedron).space())

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

  def separatePoint(self, point, withCertificate = False, scaleIntegral = True):
    cdef Vector cppPoint = point
    cdef IPOVector ipoPoint = cppPoint._vector
    cdef IPOLinearConstraint ipoConstraint
    cdef IPOInnerDescription inner
    cdef bool separated
    if withCertificate:
      separated = deref(self._polyhedron).separatePoint(ipoPoint, ipoConstraint, &inner)
    else:
      separated = deref(self._polyhedron).separatePoint(ipoPoint, ipoConstraint, NULL)
    if separated:
      if scaleIntegral:
        IPOscaleIntegral(ipoConstraint)
      constraint = _createLinearConstraint(ipoConstraint, deref(self._polyhedron).space())
    else:
      constraint = None
    if withCertificate:
      return (separated, constraint, None)
      # TODO: Convert inner...
    else:
      return constraint

