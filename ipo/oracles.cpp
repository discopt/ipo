#include "oracles.h"

#include <cassert>

using namespace soplex;

namespace ipo {

  Face::Face(const Space& space) 
    : _synced(false), _denseNormal(space.dimension()), _sparseNormal(space.dimension())
  {

  }

  Face::Face(const Space& space, const LPRowRational& inequality) : _synced(false),
    _denseNormal(space.dimension()), _sparseNormal(space.dimension())
  {
    add(inequality);
  }

  Face::Face(const Space& space, const LPRowSetRational& inequalities) : _synced(false),
    _denseNormal(space.dimension()), _sparseNormal(space.dimension())
  {
    add(inequalities);
  }

  Face::~Face()
  {

  }

  void Face::add(const LPRowRational& inequality)
  {
    if (inequality.lhs() > -infinity && inequality.rhs() < infinity && inequality.lhs() < inequality.rhs())
    {
      throw std::runtime_error("Cannot add a proper ranged row as face-defining inequality.");
    }
    _inequalities.add(inequality);
    _synced = false;
  }

  void Face::add(const LPRowSetRational& inequalities)
  {
    for (std::size_t i = 0; i < inequalities.num(); ++i)
    {
      if (inequalities.lhs(i) > -infinity && inequalities.rhs(i) < infinity
          && inequalities.lhs(i) < inequalities.rhs(i))
      {
        throw std::runtime_error("Cannot add a proper ranged row as face-defining inequality.");
      }
    }

    if (_inequalities.max() < _inequalities.num() + inequalities.num())
      _inequalities.reMax(_inequalities.num() + inequalities.num());
    _inequalities.add(inequalities);
    _synced = false;
  }

  bool Face::containsPoint(const SparseVector& point)
  {
    ensureSync();
    return scalarProduct(_denseNormal, point) == _rhs;
  }

  bool Face::containsDirection(const SparseVector& direction)
  {
    ensureSync();
    return scalarProduct(_denseNormal, direction) == 0;
  }

  void Face::ensureSync()
  {
    if (_synced)
      return;

    _denseNormal.clear();
    _rhs = 0;

    for (std::size_t i = 0; i < _inequalities.num(); ++i)
    {
      const Rational& lhs = _inequalities.lhs(i);
      const Rational& rhs = _inequalities.rhs(i);
      if (lhs <= -infinity)
      {
        _rhs += rhs;
        const SVectorRational& row = _inequalities.rowVector(i);
        for (int p = row.size() - 1; p >= 0; --p)
          _denseNormal[row.index(p)] += row.value(p);
      }
      else
      {
        assert(rhs >= -infinity);
        _rhs -= lhs;
        const SVectorRational& row = _inequalities.rowVector(i);
        for (int p = row.size() - 1; p >= 0; --p)
          _denseNormal[row.index(p)] -= row.value(p);
      }
    }

    _sparseNormal = SparseVector(_denseNormal.dim());
    assign(_sparseNormal, _denseNormal);
    _maximumNorm = 0;
    for (std::size_t p = 0; p < _sparseNormal.size(); ++p)
    {
      const Rational& x = _sparseNormal.value(p);
      if (x > 0 && x > _maximumNorm)
        _maximumNorm = x;
      if (x < 0 && x < -_maximumNorm)
        _maximumNorm = -x;
    }
  }

  OracleResult::Point::Point(SparseVector& vec)
    : vector(vec)
  {
    this->objectiveValue = minusInfinity;
  }
  
  OracleResult::Point::Point(SparseVector& vec, const Rational& value)
    : vector(vec)
  {
    this->objectiveValue = value;
  }

  OracleResult::Direction::Direction(SparseVector& vec)
    : vector(vec)
  {
    
  }

  OracleResult::OracleResult() : _heuristicLevel(0)
  {

  }

  OracleResult::~OracleResult()
  {

  }

  void OracleResult::checkConsistent()
  {
    for (std::size_t i = 1; i < points.size(); ++i)
    {
      if (points[i-1].objectiveValue < points[i].objectiveValue)
        throw std::runtime_error("Inconsistent OracleResult: Points are not sorted.");
    }
  }
  
  void OracleResult::computeMissingObjectiveValues()
  {
    for (std::size_t i = 0; i < points.size(); ++i)
    {
      if (points[i].objectiveValue == minusInfinity)
      {
        points[i].objectiveValue = scalarProduct(*_objective, points[i].vector);
      }
    }
  }

  void OracleResult::removeDuplicates()
  {
    // Check points.
    
    UniqueSparseVectors usv(_objective->dim());
    std::size_t write = 0;
    for (std::size_t read = 0; read < points.size(); ++read)
    {
      if (usv.insert(points[read].vector))
      {
        if (write != read)
          points[write] = points[read];
        ++write;
      }
    }
    while (points.size() > write)
      points.pop_back();

    // Check rays.

    usv.clear();
    write = 0;
    for (std::size_t read = 0; read < directions.size(); ++read)
    {
      if (usv.insert(directions[read].vector))
      {
        if (write != read)
          directions[write] = directions[read];
        ++write;
      }
    }
    while (directions.size() > write)
      directions.pop_back();
  }

  OracleBase::OracleBase(const std::string& name, const Space& space) :
      _name(name), _space(space), _nextOracle(NULL), _heuristicLevel(0), _currentFace(NULL)
  {

  }

  OracleBase::OracleBase(const std::string& name,
    OracleBase* nextOracle) : _name(name), _space(nextOracle->space()),
    _nextOracle(nextOracle), _heuristicLevel(nextOracle->heuristicLevel() + 1), _currentFace(NULL)
  {

  }

  void OracleBase::initializedSpace()
  {
    _tempObjective.reDim(_space.dimension(), false);
  }

  OracleBase::~OracleBase()
  {

  }

  void OracleBase::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;
    _currentFace = newFace;
    if (_nextOracle != NULL)
      _nextOracle->setFace(newFace);
  }

  void OracleBase::maximize(OracleResult& result, const SparseVector& objective, const ObjectiveBound& objectiveBound,
    std::size_t minHeuristic, std::size_t maxHeuristic)
  {
    _tempObjective.clear();
    assign(_tempObjective, objective);
    return maximize(result, _tempObjective, objectiveBound, maxHeuristic, minHeuristic);
  }

  void OracleBase::maximize(OracleResult& result, const DenseVector& objective, const ObjectiveBound& objectiveBound, 
    std::size_t minHeuristic, std::size_t maxHeuristic)
  {
    // Initialize result.

    result._objective = &objective;
    result.points.clear();
    result.directions.clear();

    bool sort = false;
    bool checkDuplicates = false;

//     std::cerr << "Oracle(" << heuristicLevel() <<") <" << name() << "> called for objective: " << objective << std::endl;

    result._heuristicLevel = maximizeController(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, 
      checkDuplicates);

    // Compute missing objective values.
    // TODO: Code duplication?

    for (std::size_t i = 0; i < result.points.size(); ++i)
    {
      if (result.points[i].objectiveValue == minusInfinity)
      {
        result.points[i].objectiveValue = scalarProduct(objective, result.points[i].vector);
      }
    }

    // If requested, sort all points by objective value.

    if (sort)
      std::sort(result.points.begin(), result.points.end());

    // If requested, remove duplicates.

    if (checkDuplicates)
      result.removeDuplicates();

//     std::cerr << "Oracle(" << heuristicLevel() <<") <" << name() << ">";
//     if (result.isFeasible())
//       std::cerr << " returned " << result.points.size() << " points, maximum = " << result.points.front().objectiveValue;
//     else if (result.isUnbounded())
//       std::cerr << " returned " << result.directions.size() << " rays.";
//     else
//       std::cerr << " claims infeasible." << std::endl;
//     std::cerr << std::endl;
  }

  std::size_t OracleBase::maximizeController(OracleResult& result, const DenseVector& objective,
    const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic, bool& sort, bool& checkDups)
  {
    assert((heuristicLevel() == 0 && _nextOracle == NULL)
      || heuristicLevel() > 0 && _nextOracle != NULL);

    // If requested, forward to next oracle.
    
    if (heuristicLevel() > maxHeuristic)
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, sort, checkDups);
    }

    // Call implementation and check whether results are satisfactory.

    maximizeImplementation(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    if (heuristicLevel() == 0 || !result.directions.empty())
      return heuristicLevel();

    if (!result.points.empty())
    {
      result.computeMissingObjectiveValues();

      // Check if one of them satisfies the objective bound.
      
      bool satisfied = false;
      for (std::size_t i = 0; i < result.points.size(); ++i)
      {
        if (objectiveBound.satisfiedBy(result.points[i].objectiveValue))
        {
          satisfied = true;
          break;
        }
      }

      if (satisfied)
        return heuristicLevel();
    }

    // Answer not satisfactory.

    if (heuristicLevel() > minHeuristic || (result.points.empty() && result.directions.empty()))
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, sort, checkDups);
    }
    else
    {
      return heuristicLevel();
    }
  }

  Face* OracleBase::currentFace()
  {
    return _currentFace;
  }

  FaceOracleBase::FaceOracleBase(const std::string& name, OracleBase* nextOracle, std::size_t maxInfeasibleIterations,
    double initialM)
    : OracleBase(name, nextOracle), _maxInfeasibleIterations(maxInfeasibleIterations), _M(initialM)
  {

  }

  FaceOracleBase::FaceOracleBase(const std::string& name, const Space& space, std::size_t maxInfeasibleIterations,
    double initialM)
    : OracleBase(name, space), _maxInfeasibleIterations(maxInfeasibleIterations), _M(initialM)
  {

  }

  FaceOracleBase::~FaceOracleBase()
  {

  }

  void FaceOracleBase::initializedSpace()
  {
    OracleBase::initializedSpace();

    _modifiedObjective.reDim(_space.dimension(), false);
  }

  void FaceOracleBase::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);
    assert(newFace == OracleBase::currentFace());

    if (newFace != NULL)
      _factor = _M / newFace->maxNorm();
  }

  Face* FaceOracleBase::currentFace()
  {
    return NULL;
  }

  std::size_t FaceOracleBase::maximizeController(OracleResult& result, const DenseVector& objective,
    const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic, bool& sort, bool& checkDups)
  {
    assert((heuristicLevel() == 0 && _nextOracle == NULL)
      || heuristicLevel() > 0 && _nextOracle != NULL);

    // If requested, forward to next oracle.

    if (heuristicLevel() > maxHeuristic)
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, sort, checkDups);
    }
    
    // Get face (note that currentFace() always yields NULL.
    
    Face* face = OracleBase::currentFace();
    if (face == NULL)
    {
      maximizeImplementation(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
    }
    else // We have a face and must tilt the objective vector.
    {
      // Compute max-norm of c.

      Rational maxObjectiveCoefficient = 0;
      for (int v = 0; v < space().dimension(); ++v)
      {
        if (objective[v] >= 0)
        {
          if (objective[v] > maxObjectiveCoefficient)
            maxObjectiveCoefficient = objective[v];
        }
        else
        {
          if (-objective[v] > maxObjectiveCoefficient)
            maxObjectiveCoefficient = -objective[v];
        }
      }
      
      // Check if c is the zero vector.

      if (maxObjectiveCoefficient == 0) 
      {
        maximizeImplementation(result, objective, objectiveBound, minHeuristic, maxHeuristic, sort, checkDups);
      }
      else
      {
        // Copy objective since modification might only affect a subset of the coordinates.

        for (std::size_t v = 0; v < space().dimension(); ++v)
        {
          _modifiedObjective[v] = objective[v];
        }
        const DenseVector& denseNormal = face->denseNormal();
        const SparseVector& sparseNormal = face->sparseNormal();

        // Loop that increases M as long as results are not feasible.

        OracleResult unrestrictedResult;
        bool verifiedNonEmpty = false;
        ObjectiveBound modifiedObjectiveBound;
        Rational modifiedObjectiveShift;
        for (std::size_t iteration = 0; iteration <= _maxInfeasibleIterations; ++iteration)
        {
          if (iteration == _maxInfeasibleIterations)
          {
            // Setup objective via face's normal vector.

            for (std::size_t v = 0; v < space().dimension(); ++v)
              _modifiedObjective[v] = 0;
            for (int p = sparseNormal.size() - 1; p >= 0; --p)
              _modifiedObjective[sparseNormal.index(p)] = sparseNormal.value(p);
            modifiedObjectiveBound.value = face->rhs();
            modifiedObjectiveBound.strict = false;
          }
          else
          {
            // Setup objective via tilting.

            Rational scaling = _factor / maxObjectiveCoefficient;
            modifiedObjectiveShift = scaling * face->rhs();
            modifiedObjectiveBound.value = objectiveBound.value + modifiedObjectiveShift;
            modifiedObjectiveBound.strict = objectiveBound.strict;
            for (int p = sparseNormal.size() - 1; p >= 0; --p)
            {
              int v = sparseNormal.index(p);
              _modifiedObjective[v] = objective[v] + scaling * sparseNormal.value(p);
            }
          }
          
          // Carry out the optimization over P.
          
          maximizeImplementation(unrestrictedResult, _modifiedObjective, modifiedObjectiveBound, minHeuristic, maxHeuristic, 
            sort, checkDups);

          // Check results.
          
          if (!unrestrictedResult.points.empty())
          {
            bool satisfied = false;
            std::size_t bestIndex = std::numeric_limits<std::size_t>::max();
            Rational bestValue = -infinity;
            for (std::size_t i = 0; i < unrestrictedResult.points.size(); ++i)
            {
              SparseVector& vector = unrestrictedResult.points[i].vector;
              if (!face->containsPoint(vector))
                continue;
              
              Rational objectiveValue = scalarProduct(objective, vector);
              if (objectiveBound.satisfiedBy(objectiveValue))
              {
                result.points.push_back(OracleResult::Point(vector, objectiveValue));
                satisfied = true;
              }
              else if (objectiveValue > bestValue)
              {
                bestIndex = i;
                bestValue = objectiveValue;
              }
            }

            // If we found a good feasible point, we can return.

            if (satisfied)
              return heuristicLevel();

            // If feasible points were found, but they were too bad, we still add the best for feasibility, aborting the loop.

            if (bestIndex < std::numeric_limits<std::size_t>::max())
            {
              result.points.push_back(OracleResult::Point(unrestrictedResult.points[bestIndex].vector, bestValue));
              break;
            }
          }
          else if (!unrestrictedResult.directions.empty())
          {
            for (std::size_t i = 0; i < unrestrictedResult.directions.size(); ++i)
            {
              SparseVector& vector = unrestrictedResult.directions[i].vector;
              if (!face->containsDirection(vector))
                continue;

              result.directions.push_back(OracleResult::Direction(vector));
            }
            
            if (!result.directions.empty())
              return heuristicLevel();
          }
          else
          {
            // Infeasible: Forward to next oracle.
            
            break;
          }

          _factor *= 2; // Increase current scaling.
          _M *= 2; // Increase scaling for future calls.
        }
      }
    }

    // Answer not satisfactory.

    if (heuristicLevel() > minHeuristic || (result.points.empty() && result.directions.empty()))
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, sort, checkDups);
    }
    else
    {
      return heuristicLevel();
    }
  }

} /* namespace ipo */
