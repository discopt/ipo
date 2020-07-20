#include <ipo/oracles_mip.hpp>

// #define IPO_DEBUG_ORACLES_MIP // Uncomment to debug this file.

#if defined (IPO_WITH_GMP)

#include "reconstruct.hpp"

namespace ipo
{

  RationalMIPExtender::RationalMIPExtender(const std::vector<bool>& integrality,
    const std::vector<std::pair<double, double>>& bounds)
    : _integrality(integrality), _currentFace(std::make_shared<ipo::Constraint<rational>>(
    ipo::alwaysSatisfiedConstraint<rational>()))
  {
    _hasContinuous = false;
    for (bool isIntegral : _integrality)
      _hasContinuous = _hasContinuous || !isIntegral;

//     if (!_hasContinuous)
//       throw std::runtime_error("Pure IP!");

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
      if (bounds[i].first == -std::numeric_limits<double>::infinity())
        mpq_set_d(_originalLowerBounds[i], -soplex::infinity);
      else
        mpq_reconstruct(_originalLowerBounds[i], bounds[i].first);
      mpq_init(_originalUpperBounds[i]);
      if (bounds[i].second == std::numeric_limits<double>::infinity())
        mpq_set_d(_originalUpperBounds[i], soplex::infinity);
      else
        mpq_reconstruct(_originalUpperBounds[i], bounds[i].second);
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
    delete[] _indices;
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
      mpq_class x = reconstruct(iter.second);
      mpq_set(_coefficients[i], x.get_mpq_t());
      ++i;
    }

    if (constraint.type() == LESS_OR_EQUAL)
      mpq_set_d(rationalLhs, -soplex::infinity);
    else
      mpq_reconstruct(rationalLhs, constraint.lhs());
    if (constraint.type() == GREATER_OR_EQUAL)
      mpq_set_d(rationalRhs, soplex::infinity);
    else
      mpq_reconstruct(rationalRhs, constraint.rhs());

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

  OptimizationOracle<rational>::Response RationalMIPExtender::maximizeDouble(
    std::shared_ptr<OptimizationOracle<double>> approximateOracle,
    const double* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    std::vector<rational> exactObjectiveVector(_integrality.size());
    for (std::size_t i = 0; i < _integrality.size(); ++i)
      exactObjectiveVector[i] = objectiveVector[i];

    return maximize(approximateOracle, &exactObjectiveVector[0], query);
  }

  void RationalMIPExtender::setZeroObjective()
  {
    mpq_set_d(_coefficients[0], 0.0);
    for (std::size_t i = 0; i < _integrality.size(); ++i)
      _spx.changeObjRational(i, _coefficients[0]);
  }

  void RationalMIPExtender::prepareRay()
  {
    // Restore all integral columns to their original bounds.
    for (std::size_t v = 0; v < _integrality.size(); ++v)
    {
      if (_integrality[v])
        _spx.changeBoundsRational(v, _originalLowerBounds[v], _originalUpperBounds[v]);
    }
  }

  void RationalMIPExtender::extractRay(OptimizationOracle<rational>::Response& result)
  {
    _spx.getPrimalRayRational(_coefficients, _integrality.size());
    auto ray = std::make_shared<sparse_vector<rational>>();
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      if (_coefficients[i] != 0)
        ray->push_back(i, rational(_coefficients[i]));
    }
    result.rays.push_back(OptimizationOracle<rational>::Response::Ray(ray));
  }

  void RationalMIPExtender::preparePoint(
    const OptimizationOracle<double>::Response::Point& approximatePoint)
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
    for (const auto& iter : *approximatePoint.vector)
    {
      if (_integrality[iter.first])
      {
        mpq_set_d(_coefficients[0], round(iter.second));
        _spx.changeBoundsRational(iter.first, _coefficients[0], _coefficients[0]);
      }
    }
  }
  
  void RationalMIPExtender::extractPoint(OptimizationOracle<rational>::Response& result,
    const rational* objectiveVector)
  {
    _spx.getPrimalRational(_coefficients, _integrality.size());
    auto point = std::make_shared<sparse_vector<rational>>();
    rational objectiveValue = 0;
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      if (mpq_sgn(_coefficients[i]) != 0)
      {
        rational x(mpq_class(_coefficients[i]));
        point->push_back(i, x);
        if (objectiveVector)
          objectiveValue += objectiveVector[i] * x;
      }
    }
    if (!objectiveVector)
      objectiveValue = rational(*_spx.objValueRational().getMpqPtr());
    result.points.push_back(OptimizationOracle<rational>::Response::Point(point, objectiveValue));
  }

  OptimizationOracle<rational>::Response RationalMIPExtender::maximize(
    std::shared_ptr<OptimizationOracle<double>> approximateOracle,
    const rational* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    std::vector<double> approximateObjectiveVector(_integrality.size());
    for (std::size_t v = 0; v < _integrality.size(); ++v)
      approximateObjectiveVector[v] = (double)objectiveVector[v];

    OptimizationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumSolutions = query.maxNumSolutions;
    approximateQuery.minObjectiveValue = query.minObjectiveValue.approximation();
    approximateQuery.timeLimit = query.timeLimit;
    
    OptimizationOracle<double>::Response approximateResponse = approximateOracle->maximize(
      &approximateObjectiveVector[0], approximateQuery);

#if defined(IPO_DEBUG_ORACLES_MIP)
    std::cout << "RationalMIPExtender::maximize. Approx. response: " << approximateResponse
      << std::endl;
#endif /* IPO_DEBUG_ORACLES_MIP */

    OptimizationOracle<rational>::Response response;
    response.hitTimeLimit = approximateResponse.hitTimeLimit;
    response.hasDualBound = approximateResponse.hasDualBound;
    response.dualBound = approximateResponse.dualBound;

    // Set LP objective.
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_set(_coefficients[0], objectiveVector[i].get_mpq_t());
      _spx.changeObjRational(i, _coefficients[0]);
    }

    if (approximateResponse.outcome == OPTIMIZATION_INFEASIBLE)
      response.outcome = OPTIMIZATION_INFEASIBLE;
    else if (approximateResponse.outcome == OPTIMIZATION_UNBOUNDED)
    {
      prepareRay();
      soplex::SPxSolver::Status status = _spx.solve();

      // Convert the ray by computing one for the LP relaxation.
      if (status == soplex::SPxSolver::UNBOUNDED)
      {
        extractRay(response);
        response.outcome = OPTIMIZATION_UNBOUNDED;
      }
      else
      {
        std::stringstream ss;
        ss << "Error in RationalMIPExtender::maximize. Unbounded approximate oracle with ray LP status "
          << status << '.';
        throw std::runtime_error(ss.str());
      }

      // Also convert all points but without objective function.
      if (!approximateResponse.points.empty())
      {
        setZeroObjective();

        for (const auto& point : approximateResponse.points)
        {
          preparePoint(point);
          soplex::SPxSolver::Status status = _spx.solve();

          if (status == soplex::SPxSolver::OPTIMAL)
            extractPoint(response, objectiveVector);
          else
          {
            std::stringstream ss;
            ss << "Error in RationalMIPExtender::maximize. Unbounded approximate oracle with point LP status "
              << status << '.';
            throw std::runtime_error(ss.str());
          }
        }
      }
    }
    else
    {
      assert(response.rays.empty());

      response.outcome = approximateResponse.outcome;
      response.primalBound = query.minObjectiveValue;
      for (const auto& point : approximateResponse.points)
      {
        preparePoint(point);
        soplex::SPxSolver::Status status = _spx.solve();

        if (status == soplex::SPxSolver::OPTIMAL)
        {
          extractPoint(response, nullptr);
          const rational& value = response.points.back().objectiveValue;
          if (response.hasDualBound && value > response.dualBound)
            response.dualBound = value;
          if (response.outcome == OPTIMIZATION_FEASIBLE && value > response.primalBound)
            response.primalBound = value;
        }
        else if (status == soplex::SPxSolver::UNBOUNDED)
        {
          if (!response.rays.empty())
            continue;

          extractRay(response);
          response.outcome = OPTIMIZATION_UNBOUNDED;
          response.primalBound = 0;
        }
        else if (status != soplex::SPxSolver::INFEASIBLE)
        {
          std::stringstream ss;
            ss << "Error in RationalMIPExtender::maximize. Optimal approximate oracle with point LP status "
                << status << '.';
            throw std::runtime_error(ss.str());
        }
      }

      if (response.outcome == OPTIMIZATION_FEASIBLE
        && approximateResponse.outcome == OPTIMIZATION_FEASIBLE
        && approximateResponse.primalBound == approximateQuery.minObjectiveValue)
      {
        response.primalBound = query.minObjectiveValue;
      }
    }

#if defined(IPO_DEBUG_ORACLES_MIP)
    std::cout << "RationalMIPExtender::maximize. Exact response: " << response << std::endl;
#endif /* IPO_DEBUG_ORACLES_MIP */

    return response;
  }

  RationalMIPExtendedOptimizationOracle::RationalMIPExtendedOptimizationOracle(
    RationalMIPExtender* extender,
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

  RationalMIPExtendedOptimizationOracle::~RationalMIPExtendedOptimizationOracle()
  {

  }

  OptimizationOracle<rational>::Response RationalMIPExtendedOptimizationOracle::maximizeDouble(
    const double* objectiveVector,
    const OptimizationOracle<rational>::Query& query)
  {
    _extender->setFace(_face);
    return _extender->maximizeDouble(_approximateOracle, objectiveVector, query);
  }

  OptimizationOracle<rational>::Response RationalMIPExtendedOptimizationOracle::maximize(
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

  RationalMIPExtendedSeparationOracle::~RationalMIPExtendedSeparationOracle()
  {

  }

  SeparationOracle<rational>::Response RationalMIPExtendedSeparationOracle::getInitial(
    const SeparationOracle<rational>::Query& query)
  {
    SeparationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumInequalities = query.maxNumInequalities;
    approximateQuery.timeLimit = query.timeLimit;

    auto approximateResponse = _approximateOracle->getInitial(approximateQuery);
    SeparationOracle<rational>::Response response;
    response.constraints.reserve(approximateResponse.numConstraints());
    response.hitTimeLimit = approximateResponse.hitTimeLimit;
    for (const auto& constraint : approximateResponse.constraints)
      response.constraints.push_back( constraintToRational(constraint) );

    return response;
  }

  SeparationOracle<rational>::Response RationalMIPExtendedSeparationOracle::separateDouble(
    const double* vector, bool isPoint, const SeparationOracle<rational>::Query& query)
  {
    SeparationOracle<double>::Query approximateQuery;
    approximateQuery.maxNumInequalities = query.maxNumInequalities;
    approximateQuery.timeLimit = query.timeLimit;

    auto approximateResponse = _approximateOracle->separateDouble(vector, isPoint, approximateQuery);
    SeparationOracle<rational>::Response response;
    response.constraints.reserve(approximateResponse.numConstraints());
    response.hitTimeLimit = approximateResponse.hitTimeLimit;
    for (const auto& constraint : approximateResponse.constraints)
      response.constraints.push_back( constraintToRational(constraint) );

    return response;
  }

  SeparationOracle<rational>::Response RationalMIPExtendedSeparationOracle::separate(
    const rational* vector, bool isPoint, const SeparationOracle<rational>::Query& query)
  {
    std::vector<double> tempVector(space()->dimension());
    for (std::size_t v = 0; v < space()->dimension(); ++v)
      tempVector[v] = vector[v].approximation();
    return separateDouble(&tempVector[0], isPoint, query);
  }

}

#endif /* IPO_WITH_GMP */
