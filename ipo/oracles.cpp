#include "oracles.h"

#include <cassert>

#include "unique_rational_vectors.h"

using namespace soplex;

namespace ipo {

  Face::Face(std::size_t numVariables) :
      _synced(false), _worker(numVariables)
  {

  }

  Face::Face(std::size_t numVariables, const soplex::LPRowRational& inequality) :
      _synced(false), _worker(numVariables)
  {
    add(inequality);
  }

  Face::Face(std::size_t numVariables, const soplex::LPRowSetRational& inequalities) :
      _synced(false), _worker(numVariables)
  {
    add(inequalities);
  }

  Face::~Face()
  {

  }

  void Face::add(const soplex::LPRowRational& inequality)
  {
    if (inequality.lhs() > -infinity && inequality.rhs() < infinity && inequality.lhs() < inequality.rhs())
    {
      throw std::runtime_error("Cannot add a proper ranged row as face-defining inequality.");
    }
    _inequalities.add(inequality);
    _synced = false;
  }

  void Face::add(const soplex::LPRowSetRational& inequalities)
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

  void Face::ensureSync()
  {
    if (_synced)
      return;

    _worker.clear();
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
          _worker[row.index(p)] += row.value(p);
      }
      else
      {
        assert(rhs >= -infinity);
        _rhs -= lhs;
        const SVectorRational& row = _inequalities.rowVector(i);
        for (int p = row.size() - 1; p >= 0; --p)
          _worker[row.index(p)] -= row.value(p);
      }
    }

    _normal.clear();
    _normal = _worker;
    _largestAbsCoefficient = 0;
    for (int p = _normal.size() - 1; p >= 0; --p)
    {
      const Rational& x = _normal.value(p);
      if (x > 0 && x > _largestAbsCoefficient)
        _largestAbsCoefficient = x;
      if (x < 0 && x < -_largestAbsCoefficient)
        _largestAbsCoefficient = -x;
    }
  }

  Projection::Projection()
  {

  }

  Projection::Projection(const Space& space, const std::vector<std::size_t>& variableSubset)
  {
    _map.reserve(variableSubset.size());
    _shift.reserve(variableSubset.size());
    for (std::size_t v = 0; v < variableSubset.size(); ++v)
      addVariable(space, variableSubset[v]);
  }

  Projection::~Projection()
  {

  }

  void Projection::addVariable(const std::string& name, const soplex::SVectorRational& variableMap,
      const soplex::Rational& shift)
  {
    _names.push_back(name);
    _map.push_back(DSVectorRational(variableMap));
    _shift.push_back(shift);
  }

  void Projection::addVariable(const Space& space, std::size_t originalVariable,
      const soplex::Rational& shift)
  {
    _names.push_back(space[originalVariable]);
    _map.push_back(DSVectorRational(1));
    _map.back().add(originalVariable, Rational(1));
    _shift.push_back(shift);
  }

  void Projection::projectPoint(const DVectorRational& point, DSVectorRational& image) const
  {
    for (std::size_t i = 0; i < _map.size(); ++i)
    {
      Rational x = _shift[i] + _map[i] * point;
      if (x != 0)
        image.add(i, x);
    }
  }

  void OptimizationResult::reset(std::size_t numVariables)
  {
    _numVariables = numVariables;
    bestIndex = std::numeric_limits<std::size_t>::max();
    bestValue = -infinity;
    points.clear();
    directions.clear();
    objectives.clear();
  }

  DSVectorRational& OptimizationResult::newPoint()
  {
    points.push_back(new DSVectorRational);
    return *points.back();
  }

  void OptimizationResult::setFeasible(const soplex::VectorRational& objective)
  {
    objectives.resize(points.size());
    bestValue = -soplex::infinity;
    bestIndex = std::numeric_limits<std::size_t>::max();
    for (std::size_t p = 0; p < points.size(); ++p)
    {
      objectives[p] = *points[p] * objective;
      if (objectives[p] > bestValue)
      {
        bestIndex = p;
        bestValue = objectives[p];
      }
    }
  }

  void OptimizationResult::setBest(const Rational& value, std::size_t index)
  {
    bestIndex = index;
    bestValue = value;
  }

  void OptimizationResult::setInfeasible()
  {
    bestIndex = std::numeric_limits<std::size_t>::max();
    bestValue = -infinity;
  }

  void OptimizationResult::setUnbounded()
  {
    bestIndex = std::numeric_limits<std::size_t>::max();
    bestValue = infinity;
  }

  void OptimizationResult::filterDuplicates()
  {
    filterDuplicates(points);
    filterDuplicates(directions);
  }

  void OptimizationResult::filterDuplicates(std::vector<soplex::DSVectorRational*>& vectors)
  {
    /// Move vectors to unique which does the filtering.

    UniqueRationalVectors unique(_numVariables);
    for (std::size_t i = 0; i < vectors.size(); ++i)
      unique.insertFree(vectors[i]);

    /// Move everything from unique to vectors.

    vectors.resize(unique.size());
    for (std::size_t i = 0; i < unique.size(); ++i)
      vectors[i] = unique.get(i);
    unique.extractAll();
  }

#ifdef IPO_DEBUG
  void OptimizationResult::checkConsistent() const
  {
    if (points.size() != objectives.size())
      throw std::runtime_error("Numbers of points and objectives do not match!");
    if (!directions.empty() && bestValue < soplex::infinity)
      throw std::runtime_error("Not unbounded, but unbounded ray exists!");
    if (!points.empty() && (bestValue == soplex::infinity || bestValue == -soplex::infinity))
      throw std::runtime_error("Infeasible or unbounded, but point exists!");
  }

  bool OptimizationResult::hasDuplicates() const
  {
    const std::vector<DSVectorRational*>* vectors = points.empty() ? &directions : &points;
    if (vectors->size() < 2)
      return false;

    DVectorRational worker;
    worker.reDim(_numVariables, false);
    for (std::size_t s1 = 0; s1 < vectors->size(); ++s1)
    {
      const SVectorRational& v1 = *vectors->at(s1);
      worker.clear();
      worker.assign(v1);
      for (std::size_t s2 = s1 + 1; s2 < vectors->size(); ++s2)
      {
        const SVectorRational& v2 = *vectors->at(s2);
        if (v2.size() != v1.size())
          continue;
        bool dup = true;
        for (int p = v2.size() - 1; p >= 0; --p)
        {
          if (worker[v2.index(p)] != v2.value(p))
          {
            dup = false;
            break;
          }
        }
        if (dup)
          return true;
      }
    }
    return false;
  }
#endif

  OptimizationOracleBase::OptimizationOracleBase(const std::string& name, const Space& space) :
      _name(name), _space(space)
  {
    
  }

  OptimizationOracleBase::~OptimizationOracleBase()
  {

  }
  
  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::VectorRational& objective,
      bool forceOptimal)
  {
    result.reset(space().dimension());
    return run(result, objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::VectorReal& objective,
      bool forceOptimal)
  {
    result.reset(space().dimension());
    _objective.reDim(space().dimension(), false);
    _objective = objective;
    return run(result, _objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::SVectorRational& objective,
      bool forceOptimal)
  {
    result.reset(space().dimension());
    _objective.reDim(space().dimension(), false);
    _objective.clear();
    _objective.assign(objective);
    return run(result, _objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::SVectorReal& objective,
      bool forceOptimal)
  {
    result.reset(space().dimension());
    _objective.reDim(space().dimension(), false);
    _objective = objective;
    return run(result, _objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::run(OptimizationResult& result, const soplex::VectorRational& objective,
      const soplex::Rational* improveValue, bool forceOptimal)
  {
    throw std::runtime_error("You should implement OptimizationOracle::run() for a base class!");
  }

  ChainedOptimizationOracle::ChainedOptimizationOracle(FaceOptimizationOracleBase* first,
      FaceOptimizationOracleBase* second) : 
      FaceOptimizationOracleBase(first->name() + "+" + second->name(), first->space()),       
      _first(first), _second(second)
  {
    if (first->space() != second->space())
      throw std::runtime_error("ChainedOptimizationOracle: Spaces do not match.");
  }

  ChainedOptimizationOracle::~ChainedOptimizationOracle()
  {

  }

  void ChainedOptimizationOracle::run(OptimizationResult& result, const soplex::VectorRational& objective,
      const soplex::Rational* improveValue, bool forceOptimal)
  {
    if (improveValue)
      throw std::runtime_error("UNIMPLEMENTED.");

    _first->maximize(result, objective, forceOptimal);

//    std::cerr << "First oracles yields " << result.bestValue << std::endl;
//    std::cerr << " [forceOptimal = " << forceOptimal << ", result.optimal = " << result.optimal << "] " << std::flush;

    if (!forceOptimal || result.optimal)
      return;

    for (std::size_t i = 0; i < result.points.size(); ++i)
      delete result.points[i];
    result.points.clear();
    for (std::size_t i = 0; i < result.directions.size(); ++i)
      delete result.directions[i];
    result.directions.clear();

    _second->maximize(result, objective, forceOptimal);
  }

  void ChainedOptimizationOracle::faceEnabled(Face* face)
  {
    _first->setFace(face);
    _second->setFace(face);
  }

  void ChainedOptimizationOracle::faceDisabled(Face* face)
  {
    _first->setFace(NULL);
    _second->setFace(NULL);
  }

  FaceOptimizationOracleBase::FaceOptimizationOracleBase(const std::string& name,
    const Space& space) : OptimizationOracleBase(name, space), _face(NULL)
  {

  }

  FaceOptimizationOracleBase::~FaceOptimizationOracleBase()
  {

  }

  Face* FaceOptimizationOracleBase::setFace(Face* face)
  {
    if (face == _face)
      return face;

    Face* result = _face;
    if (_face != NULL)
      faceDisabled(_face);

    _face = face;
    if (_face != NULL)
      faceEnabled(_face);
    return result;
  }

//   ProjectedOptimizationOracle::ProjectedOptimizationOracle(const std::string& name, const Projection& projection,
//       OptimizationOracleBase* oracle) :
//       OptimizationOracleBase(name), _projection(projection), _oracle(oracle)
//   {
//     initialize(projection.names());
//     _liftedObjective.reDim(oracle->numVariables());
//   }
// 
//   ProjectedOptimizationOracle::~ProjectedOptimizationOracle()
//   {
// 
//   }
// 
//   void ProjectedOptimizationOracle::run(OptimizationResult& result, const VectorRational& objective,
//       const Rational* improveValue, bool forceOptimal)
//   {
//     _liftedObjective.clear();
//     for (std::size_t v = 0; v < numVariables(); ++v)
//     {
//       if (objective[v] == 0)
//         continue;
//       const SVectorRational& vector = _projection.map(v);
//       for (int p = vector.size() - 1; p >= 0; --p)
//         _liftedObjective[vector.index(p)] += vector.value(p) * objective[v];
//     }
// 
// //    std::cout << "Extension oracle has " << _oracle->numVariables() << " variables." << std::endl;
// //    std::cout << "Original objective is " << objective << std::endl;
// //    std::cout << "Lifted objective is " << _liftedObjective << std::endl;
// 
//     _oracle->maximize(_result, _liftedObjective, forceOptimal);
// 
//     if (_result.isFeasible())
//     {
//       result.reset(numVariables());
//       DVectorRational dense;
//       dense.reDim(_oracle->numVariables());
//       for (std::size_t i = 0; i < _result.points.size(); ++i)
//       {
//         dense.clear();
// //        _oracle->printVector(std::cout, _result.points[i]);
//         dense.assign(*_result.points[i]);
// //        std::cout << " -> ";
//         _projection.projectPoint(dense, result.newPoint());
// //        printVector(std::cout, result.points.back());
// //        std::cout << std::endl;
//         delete _result.points[i];
//       }
//       result.filterDuplicates();
//       result.setFeasible(objective);
//     }
//     else
//     {
//       throw std::runtime_error("NOT IMPLEMENTED.");
//     }
//   }

} /* namespace ipo */
