#include <ipo/oracles_mip.hpp>

namespace ipo
{

  RationalMIPExtender::RationalMIPExtender(const std::vector<bool>& integrality,
    const std::vector<std::pair<double, double>>& bounds)
    : _integrality(integrality), _currentFace(std::make_shared<ipo::Constraint<rational>>(
    ipo::alwaysSatisfiedConstraint<rational>()))
  {
    _spx.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _spx.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    _spx.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    _spx.setRealParam(soplex::SoPlex::OPTTOL, 0.0);
    _spx.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);

    _indices = new int[integrality.size()];
    _coefficients = new mpq_t[integrality.size()];
    _originalLowerBounds = new mpq_t[integrality.size()];
    _originalUpperBounds = new mpq_t[integrality.size()];
    for (std::size_t i = 0; i < integrality.size(); ++i)
    {
      mpq_init(_coefficients[i]);
      mpq_init(_originalLowerBounds[i]);
      mpq_set_d(_originalLowerBounds[i],
        isMinusInfinity(bounds[i].first) ? -soplex::infinity : bounds[i].first);
      mpq_init(_originalUpperBounds[i]);
      mpq_set_d(_originalUpperBounds[i],
        isPlusInfinity(bounds[i].second) ? soplex::infinity : bounds[i].second);
    }

    _spx.addColsRational(_coefficients, _originalLowerBounds, nullptr, nullptr, nullptr, nullptr,
      _integrality.size(), 0, _originalUpperBounds);
  }

  RationalMIPExtender::~RationalMIPExtender()
  {
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_clear(_coefficients[i]);
      mpq_clear(_originalLowerBounds[i]);
      mpq_clear(_originalUpperBounds[i]);
    }
    delete[] _coefficients;
    delete[] _originalLowerBounds;
    delete[] _originalUpperBounds;
  }

  void RationalMIPExtender::addConstraint(const Constraint<rational>& constraint)
  {
    mpq_t rationalLhs, rationalRhs;
    mpq_init(rationalLhs);
    mpq_init(rationalRhs);

    std::size_t i = 0;
    for (const auto& iter : constraint.vector())
    {
      _indices[i] = iter.first;
      mpq_set(_coefficients[i], iter.second.get_mpq_t());
      ++i;
    }

    mpq_set(rationalLhs, constraint.lhs().get_mpq_t());
    mpq_set(rationalRhs, constraint.rhs().get_mpq_t());

    _spx.addRowRational(&rationalLhs, _coefficients, _indices, constraint.vector().size(), &rationalRhs);

    mpq_clear(rationalLhs);
    mpq_clear(rationalRhs);
  }

  void RationalMIPExtender::addConstraint(const Constraint<double>& constraint)
  {
    mpq_t rationalLhs, rationalRhs;
    mpq_init(rationalLhs);
    mpq_init(rationalRhs);

    std::size_t i = 0;
    for (const auto& iter : constraint.vector())
    {
      _indices[i] = iter.first;
      mpq_set_d(_coefficients[i], iter.second);
      ++i;
    }

    mpq_set_d(rationalLhs, isMinusInfinity(constraint.lhs()) ? -soplex::infinity : constraint.lhs());
    mpq_set_d(rationalRhs, isPlusInfinity(constraint.rhs()) ? soplex::infinity : constraint.rhs());

    _spx.addRowRational(&rationalLhs, _coefficients, _indices, constraint.vector().size(), &rationalRhs);

    mpq_clear(rationalLhs);
    mpq_clear(rationalRhs);
  }

  void RationalMIPExtender::setFace(std::shared_ptr<Constraint<rational>> face)
  {
    if (face != _currentFace)
    {
      if (!_currentFace->isAlwaysSatisfied())
        _spx.removeRowRational(_spx.numRowsRational() - 1);
      if (!face->isAlwaysSatisfied())
        addConstraint(*face);

      _currentFace = face;
    }
  }

  OptimizationOracle<rational>::Result RationalMIPExtender::maximizeDouble(
    std::shared_ptr<OptimizationOracle<double>> approximateOracle,
    const double* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    // Set objective for LP.
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_set_d(_coefficients[0], objectiveVector[i]);
      _spx.changeObjRational(i, _coefficients[0]);
    }
    
    // Query
    OptimizationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumSolutions = query.maxNumSolutions;
    approximateQuery.minObjectiveValue = query.minObjectiveValue.approximation();
    approximateQuery.timeLimit = query.timeLimit;

    return solve(approximateOracle->maximize(objectiveVector, approximateQuery));
  }

  OptimizationOracle<rational>::Result RationalMIPExtender::maximize(
    std::shared_ptr<OptimizationOracle<double>> approximateOracle,
    const rational* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    std::vector<double> temp(_integrality.size());

    // Set objective for LP.
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_set(_coefficients[0], objectiveVector[i].get_mpq_t());
      _spx.changeObjRational(i, _coefficients[0]);
      temp[i] = objectiveVector[i].approximation();
    }

    // Query

    OptimizationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumSolutions = query.maxNumSolutions;
    approximateQuery.minObjectiveValue = query.minObjectiveValue.approximation();
    approximateQuery.timeLimit = query.timeLimit;

    return solve(approximateOracle->maximize(&temp[0], approximateQuery));
  }

  OptimizationOracle<rational>::Result RationalMIPExtender::solve(
    const OptimizationOracle<double>::Result& approximateResult)
  {
    OptimizationOracle<rational>::Result result;

    if (!approximateResult.rays.empty())
    {
      // Restore all integral columns to their original bounds.
      for (std::size_t v = 0; v < _integrality.size(); ++v)
      {
        if (_integrality[v])
          _spx.changeBoundsRational(v, _originalLowerBounds[v], _originalUpperBounds[v]);
      }

      // Solve the LP.
      soplex::SPxSolver::Status status = _spx.solve();

      if (status == soplex::SPxSolver::UNBOUNDED)
      {
        _spx.getPrimalRayRational(_coefficients, _integrality.size());
        sparse_vector<rational> ray;
        for (std::size_t i = 0; i < _integrality.size(); ++i)
        {
          if (_coefficients[i] != 0)
            ray.push_back(i, rational(_coefficients[i]));
        }
        result.rays.clear();
        result.rays.push_back(OptimizationOracle<rational>::Result::Ray(std::move(ray)));
      }
      else
      {
        std::stringstream ss;
        ss << "Error in RationalMIPExtender::solve: Unbounded case yields SoPlex status " << status
          << '.'; 
        throw std::runtime_error(ss.str());
      }
    }

    std::vector<OptimizationOracle<rational>::Result::Point> points;
    for (const auto& point : approximateResult.points)
    {
      // First, fix all integral columns to zero.
      for (std::size_t v = 0; v < _integrality.size(); ++v)
      {
        if (_integrality[v])
        {
          mpq_set_si(_coefficients[0], 0, 1);
          _spx.changeBoundsRational(v, _coefficients[0], _coefficients[0]);
        }
      }

      // Now go through the current solution vector and fix the integral variables.
      for (const auto& iter : point.vector)
      {
        if (_integrality[iter.first])
        {
          mpq_set_d(_coefficients[0], round(iter.second));
          _spx.changeBoundsRational(iter.first, _coefficients[0], _coefficients[0]);
        }
      }

      // Solve the LP.
      soplex::SPxSolver::Status status = _spx.solve();

      if (status == soplex::SPxSolver::OPTIMAL)
      {
        _spx.getPrimalRational(_coefficients, _integrality.size());
        sparse_vector<rational> point;
        for (std::size_t i = 0; i < _integrality.size(); ++i)
        {
          if (_coefficients[i] != 0)
            point.push_back(i, rational(mpq_class(_coefficients[i])));
        }
        points.push_back(OptimizationOracle<rational>::Result::Point(std::move(point),
          rational(*_spx.objValueRational().getMpqPtr())));
      }
      else if (status == soplex::SPxSolver::UNBOUNDED)
      {
        // If we have a ray already, then we ignore this point.
        if (!result.rays.empty())
          continue;

        // Otherwise, create a ray.
        _spx.getPrimalRayRational(_coefficients, _integrality.size());
        sparse_vector<rational> ray;
        for (std::size_t i = 0; i < _integrality.size(); ++i)
        {
          if (_coefficients[i] != 0)
            ray.push_back(i, rational(_coefficients[i]));
        }
        result.rays.push_back(OptimizationOracle<rational>::Result::Ray(std::move(ray)));
      }
      else
      {
        std::stringstream ss;
        ss << "Error in MakeRationalSolver::solve: bounded case yields SoPlex status " << status
          << '.'; 
        throw std::runtime_error(ss.str());
      }
    }
    result.points.swap(points);
    
    if (!result.rays.empty())
    {
      result.dualBound = plusInfinity();
    }
    else if (result.points.empty())
      result.dualBound = minusInfinity();
    else
    {
      result.dualBound = approximateResult.dualBound;
      for (const auto& point : result.points)
      {
        if (point.objectiveValue > result.dualBound)
          result.dualBound = point.objectiveValue;
      }
    }

    return result;
  }

  RationalMIPExtendedOptimizationOracle::RationalMIPExtendedOptimizationOracle(
    std::shared_ptr<RationalMIPExtender> extender,
    std::shared_ptr<OptimizationOracle<double>> approximateOracle,
    std::shared_ptr<Constraint<rational>> face)
    : OptimizationOracle<rational>("Rational " + approximateOracle->name()), _extender(extender),
    _approximateOracle(approximateOracle), _face(face)
  {
    assert(_extender);
    assert(_approximateOracle);
    assert(_face);

    _space = approximateOracle->space();
  }

  OptimizationOracle<rational>::Result RationalMIPExtendedOptimizationOracle::maximizeDouble(
    const double* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    _extender->setFace(_face);
    return _extender->maximizeDouble(_approximateOracle, objectiveVector, query);
  }

  OptimizationOracle<rational>::Result RationalMIPExtendedOptimizationOracle::maximize(
    const rational* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    _extender->setFace(_face);
    return _extender->maximize(_approximateOracle, objectiveVector, query);
  }

  RationalMIPExtendedSeparationOracle::RationalMIPExtendedSeparationOracle(
    std::shared_ptr<SeparationOracle<double>> approximateOracle,
    std::shared_ptr<Constraint<rational>> face)
    : SeparationOracle<rational>("Rational " + approximateOracle->name()),
    _approximateOracle(approximateOracle), _face(face)
  {
    _space = approximateOracle->space();
  }

  SeparationOracle<rational>::Result RationalMIPExtendedSeparationOracle::getInitial(
    const SeparationOracle<rational>::Query& query)
  {
    SeparationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumInequalities = query.maxNumInequalities;
    approximateQuery.timeLimit = query.timeLimit;

    auto approximateResult = _approximateOracle->getInitial(approximateQuery);
    SeparationOracle<rational>::Result result;
    result.constraints.reserve(approximateResult.numConstraints());
    result.hitTimeLimit = approximateResult.hitTimeLimit;
    for (const auto& constraint : approximateResult.constraints)
      result.constraints.push_back( constraintToRational(constraint) );

    return result;
  }

  SeparationOracle<rational>::Result RationalMIPExtendedSeparationOracle::separateDouble(
    const double* vector, bool isPoint, const SeparationOracle<rational>::Query& query)
  {
    SeparationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumInequalities = query.maxNumInequalities;
    approximateQuery.timeLimit = query.timeLimit;

    auto approximateResult = _approximateOracle->separateDouble(vector, isPoint, approximateQuery);
    SeparationOracle<rational>::Result result;
    result.constraints.reserve(approximateResult.numConstraints());
    result.hitTimeLimit = approximateResult.hitTimeLimit;
    for (const auto& constraint : approximateResult.constraints)
      result.constraints.push_back( constraintToRational(constraint) );

    return result;
  }

  SeparationOracle<rational>::Result RationalMIPExtendedSeparationOracle::separate(
    const rational* vector, bool isPoint, const SeparationOracle<rational>::Query& query)
  {
    double* temp = new double[space()->dimension()];
    for (std::size_t v = 0; v < space()->dimension(); ++v)
      temp[v] = vector[v].approximation();
    auto result = separateDouble(temp, isPoint, query);
    delete[] temp;
    return result;
  }

}
