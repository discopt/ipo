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

  Projection::Projection(const OptimizationOracleBase* oracle, const std::vector<std::size_t>& variableSubset)
  {
    _map.reserve(variableSubset.size());
    _shift.reserve(variableSubset.size());
    for (std::size_t v = 0; v < variableSubset.size(); ++v)
      addVariable(oracle, variableSubset[v]);
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

  void Projection::addVariable(const OptimizationOracleBase* oracle, std::size_t originalVariable,
      const soplex::Rational& shift)
  {
    _names.push_back(oracle->variableName(originalVariable));
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
    rays.clear();
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
    filterDuplicates(rays);
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

#ifndef NDEBUG
  void OptimizationResult::checkConsistent() const
  {
    if (points.size() != objectives.size())
      throw std::runtime_error("Numbers of points and objectives do not match!");
    if (!rays.empty() && bestValue < soplex::infinity)
      throw std::runtime_error("Not unbounded, but unbounded ray exists!");
    if (!points.empty() && (bestValue == soplex::infinity || bestValue == -soplex::infinity))
      throw std::runtime_error("Infeasible or unbounded, but point exists!");
  }

  bool OptimizationResult::hasDuplicates() const
  {
    const std::vector<DSVectorRational*>* vectors = points.empty() ? &rays : &points;
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

  OptimizationOracleBase::OptimizationOracleBase(const std::string& name) :
      _name(name)
  {
#ifndef NDEBUG
    _initialized = false;
#endif
  }

  OptimizationOracleBase::~OptimizationOracleBase()
  {

  }

  void OptimizationOracleBase::initialize(const std::vector<std::string>& variableNames)
  {
#ifndef NDEBUG
    _initialized = true;
#endif
    _variableNames = variableNames;
  }

  void OptimizationOracleBase::initialize(const OptimizationOracleBase* oracle)
  {
#ifndef NDEBUG
    _initialized = true;
#endif
    _variableNames.resize(oracle->numVariables());
    for (std::size_t v = 0; v < _variableNames.size(); ++v)
      _variableNames[v] = oracle->variableName(v);
  }

  void OptimizationOracleBase::printRow(std::ostream& stream, const LPRowRational& row) const
  {
    const Rational& lhs = row.lhs();
    const Rational& rhs = row.rhs();
    printRow(stream, lhs > -soplex::infinity / 2 ? &lhs : NULL, rhs < soplex::infinity / 2 ? &rhs : NULL,
        row.rowVector());
  }

  void OptimizationOracleBase::printRow(std::ostream& stream, const LPRowSetRational& rows, std::size_t index) const
  {
    assert(index < rows.num());
    const Rational& lhs = rows.lhs(index);
    const Rational& rhs = rows.rhs(index);
    printRow(stream, lhs > -soplex::infinity / 2 ? &lhs : NULL, rhs < soplex::infinity / 2 ? &rhs : NULL,
        rows.rowVector(index));
  }

  void OptimizationOracleBase::printRows(std::ostream& stream, const LPRowSetRational& rows) const
  {
    for (int i = 0; i < rows.num(); ++i)
    {
      const Rational& lhs = rows.lhs(i);
      const Rational& rhs = rows.rhs(i);
      printRow(stream, lhs > -soplex::infinity / 2 ? &lhs : NULL, rhs < soplex::infinity / 2 ? &rhs : NULL,
          rows.rowVector(i));
      stream << "\n";
    }
  }

  void OptimizationOracleBase::printRow(std::ostream& stream, const Rational* lhs, const Rational* rhs,
      const SVectorRational& vector) const
  {
    bool equation = lhs != NULL && rhs != NULL && *lhs == *rhs;
    if (lhs != NULL)
      stream << *lhs << (equation ? " == " : " <= ");
    bool first = true;
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
      std::size_t v = vector.index(i);
      const Rational& x = vector.value(i);
      if (x < 0)
        stream << (first ? "-" : " - ");
      else if (!first)
        stream << " + ";
      if (x != 1 && x != -1)
      {
        if (x > 0)
          stream << x << '*';
        else
          stream << (-x) << '*';
      }
      stream << variableName(v);
      first = false;
    }
    if (first)
      stream << '0';
    if (rhs != NULL)
      stream << (equation ? " == " : " <= ") << *rhs;
  }

  void OptimizationOracleBase::printVector(std::ostream& stream, const soplex::SVectorRational* vector) const
  {
    bool delimit = false;
    for (std::size_t i = 0; i < vector->size(); ++i)
    {
      std::size_t v = vector->index(i);
      const soplex::Rational& x = vector->value(i);
      stream << (delimit ? ", " : "(") << variableName(v) << "=" << x;
      delimit = true;
    }
    if (delimit)
      stream << ")";
    else
      stream << "<zero-vector>";
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::VectorRational& objective,
      bool forceOptimal)
  {
    result.reset(numVariables());
    return run(result, objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::VectorReal& objective,
      bool forceOptimal)
  {
    result.reset(numVariables());
    _objective.reDim(numVariables(), false);
    _objective = objective;
    return run(result, _objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::SVectorRational& objective,
      bool forceOptimal)
  {
    result.reset(numVariables());
    _objective.reDim(numVariables(), false);
    _objective.clear();
    _objective.assign(objective);
    return run(result, _objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::maximize(OptimizationResult& result, const soplex::SVectorReal& objective,
      bool forceOptimal)
  {
    result.reset(numVariables());
    _objective.reDim(numVariables(), false);
    _objective = objective;
    return run(result, _objective, NULL, forceOptimal);
  }

  void OptimizationOracleBase::improve(OptimizationResult& result, const soplex::VectorRational& objective,
      const soplex::Rational& value, bool forceOptimal)
  {
    result.reset(numVariables());
    return run(result, objective, &value, forceOptimal);
  }

  void OptimizationOracleBase::improve(OptimizationResult& result, const soplex::VectorReal& objective,
      const soplex::Rational& value, bool forceOptimal)
  {
    result.reset(numVariables());
    _objective.reDim(numVariables(), false);
    _objective = objective;
    return run(result, _objective, &value, forceOptimal);
  }

  void OptimizationOracleBase::improve(OptimizationResult& result, const soplex::SVectorRational& objective,
      const soplex::Rational& value, bool forceOptimal)
  {
    result.reset(numVariables());
    _objective.reDim(numVariables(), false);
    _objective.clear();
    _objective.assign(objective);
    return run(result, _objective, &value, forceOptimal);
  }

  void OptimizationOracleBase::improve(OptimizationResult& result, const soplex::SVectorReal& objective,
      const soplex::Rational& value, bool forceOptimal)
  {
    result.reset(numVariables());
    _objective.reDim(numVariables(), false);
    _objective.clear();
    _objective = objective;
    return run(result, _objective, &value, forceOptimal);
  }

  void OptimizationOracleBase::run(OptimizationResult& result, const soplex::VectorRational& objective,
      const soplex::Rational* improveValue, bool forceOptimal)
  {
    throw std::runtime_error("You should implement OptimizationOracle::run() for a base class!");
  }

  ChainedOptimizationOracle::ChainedOptimizationOracle(OptimizationOracleBase* first, OptimizationOracleBase* second) :
      OptimizationOracleBase(first->name() + "+" + second->name()), _first(first), _second(second)
  {
    if (first->numVariables() != second->numVariables())
      throw std::runtime_error("ChainedOptimizationOracle: Ambient dimensions do not match.");

    std::size_t n = first->numVariables();

    for (std::size_t v = 0; v < n; ++v)
    {
      if (first->variableName(v) != second->variableName(v))
      {
        throw std::runtime_error(
            "ChainedOptimizationOracle: Variables \"" + first->variableName(v) + "\" and \"" + second->variableName(v)
                + "\" do not match.");
      }
    }

    initialize(first);
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
    for (std::size_t i = 0; i < result.rays.size(); ++i)
      delete result.rays[i];
    result.rays.clear();

    _second->maximize(result, objective, forceOptimal);
  }

  FaceOptimizationOracleBase::FaceOptimizationOracleBase(const std::string& name) :
      OptimizationOracleBase(name), _face(NULL)
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

  ProjectedOptimizationOracle::ProjectedOptimizationOracle(const std::string& name, const Projection& projection,
      OptimizationOracleBase* oracle) :
      OptimizationOracleBase(name), _projection(projection), _oracle(oracle)
  {
    initialize(projection.names());
    _liftedObjective.reDim(oracle->numVariables());
  }

  ProjectedOptimizationOracle::~ProjectedOptimizationOracle()
  {

  }

  void ProjectedOptimizationOracle::run(OptimizationResult& result, const VectorRational& objective,
      const Rational* improveValue, bool forceOptimal)
  {
    _liftedObjective.clear();
    for (std::size_t v = 0; v < numVariables(); ++v)
    {
      if (objective[v] == 0)
        continue;
      const SVectorRational& vector = _projection.map(v);
      for (int p = vector.size() - 1; p >= 0; --p)
        _liftedObjective[vector.index(p)] += vector.value(p) * objective[v];
    }

//    std::cout << "Extension oracle has " << _oracle->numVariables() << " variables." << std::endl;
//    std::cout << "Original objective is " << objective << std::endl;
//    std::cout << "Lifted objective is " << _liftedObjective << std::endl;

    _oracle->maximize(_result, _liftedObjective, forceOptimal);

    if (_result.isFeasible())
    {
      result.reset(numVariables());
      DVectorRational dense;
      dense.reDim(_oracle->numVariables());
      for (std::size_t i = 0; i < _result.points.size(); ++i)
      {
        dense.clear();
//        _oracle->printVector(std::cout, _result.points[i]);
        dense.assign(*_result.points[i]);
//        std::cout << " -> ";
        _projection.projectPoint(dense, result.newPoint());
//        printVector(std::cout, result.points.back());
//        std::cout << std::endl;
        delete _result.points[i];
      }
      result.filterDuplicates();
      result.setFeasible(objective);
    }
    else
    {
      throw std::runtime_error("NOT IMPLEMENTED.");
    }
  }

} /* namespace ipo */
