#include <ipo/make_rational_soplex.hpp>

#include <iostream>

namespace ipo
{

  MakeRationalSolver::MakeRationalSolver(const std::vector<bool>& integrality,
      const std::vector<Constraint>& constraints)
  {
    _integrality = integrality;
    _spx.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MAXIMIZE);
    _spx.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
    _spx.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
    _spx.setRealParam(soplex::SoPlex::OPTTOL, 0.0);

    _indices = new int[integrality.size()];
    _coefficients = new mpq_t[integrality.size()];
    _originalLowerBounds = new mpq_t[integrality.size()];
    _originalUpperBounds = new mpq_t[integrality.size()];
    for (std::size_t i = 0; i < integrality.size(); ++i)
    {
      mpq_init(_coefficients[i]);
      mpq_set_d(_coefficients[i], 0.0);
      mpq_init(_originalLowerBounds[i]);
      mpq_set_d(_originalLowerBounds[i], -soplex::infinity);
      mpq_init(_originalUpperBounds[i]);
      mpq_set_d(_originalUpperBounds[i], soplex::infinity);
    }

    // Process bounds.
    for (const auto& constraint : constraints)
    {
      if (constraint.vector.size() == 1)
      {
        std::size_t coordinate = constraint.vector.coordinate(0);
        if (constraint.lhs.isFinite())
        {
          if (mpq_cmp(constraint.lhs.rational.get_mpq_t(), _originalLowerBounds[coordinate]) > 0)
            mpq_set(_originalLowerBounds[coordinate], constraint.lhs.rational.get_mpq_t());
        }
        if (constraint.rhs.isFinite())
        {
          if (mpq_cmp(constraint.rhs.rational.get_mpq_t(), _originalUpperBounds[coordinate]) < 0)
            mpq_set(_originalUpperBounds[coordinate], constraint.rhs.rational.get_mpq_t());
        }
      }
    }

    _spx.addColsRational(_coefficients, _originalLowerBounds, nullptr, nullptr, nullptr, nullptr,
      _integrality.size(), 0, _originalUpperBounds);

    // Process other constraints.
    mpq_t rationalLhs, rationalRhs;
    mpq_init(rationalLhs);
    mpq_init(rationalRhs);
    for (const auto& constraint : constraints)
    {
      if (constraint.vector.size() != 1)
      {
        for (std::size_t i = 0; i < constraint.vector.size(); ++i)
        {
          _indices[i] = constraint.vector.coordinate(i);
          mpq_set(_coefficients[i], constraint.vector.rational(i).get_mpq_t());
        }

        mpq_set(rationalLhs, constraint.lhs.rational.get_mpq_t());
        mpq_set(rationalRhs, constraint.rhs.rational.get_mpq_t());

        _spx.addRowRational(&rationalLhs, _coefficients, _indices, constraint.vector.size(), &rationalRhs);
      }
    }
    mpq_clear(rationalLhs);
    mpq_clear(rationalRhs);
  }

  MakeRationalSolver::~MakeRationalSolver()
  {
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_clear(_coefficients[i]);
      mpq_clear(_originalLowerBounds[i]);
      mpq_clear(_originalUpperBounds[i]);
    }
    delete[]_coefficients;
    delete[] _originalLowerBounds;
    delete[] _originalUpperBounds;
  }

  void MakeRationalSolver::addRow(const Constraint& constraint)
  {
    mpq_t rationalLhs, rationalRhs;
    mpq_init(rationalLhs);
    mpq_init(rationalRhs);

    for (std::size_t i = 0; i < constraint.vector.size(); ++i)
    {
      _indices[i] = constraint.vector.coordinate(i);
      mpq_set(_coefficients[i], constraint.vector.rational(i).get_mpq_t());
    }

    mpq_set(rationalLhs, constraint.lhs.rational.get_mpq_t());
    mpq_set(rationalRhs, constraint.rhs.rational.get_mpq_t());

    _spx.addRowRational(&rationalLhs, _coefficients, _indices, constraint.vector.size(), &rationalRhs);

    mpq_clear(rationalLhs);
    mpq_clear(rationalRhs);
  }

  void MakeRationalSolver::removeLastRow()
  {
    _spx.removeRowRational(_spx.numRowsRational() - 1);
  }

  void MakeRationalSolver::setObjective(const double* objectiveVector)
  {
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_set_d(_coefficients[0], objectiveVector[i]);
      _spx.changeObjRational(i, _coefficients[0]);
    }
  }

  void MakeRationalSolver::setObjective(const mpq_class* objectiveVector)
  {
    for (std::size_t i = 0; i < _integrality.size(); ++i)
    {
      mpq_set(_coefficients[0], objectiveVector[i].get_mpq_t());
      _spx.changeObjRational(i, _coefficients[0]);
    }
  }

  void MakeRationalSolver::solve(OptimizationOracle::Result& result)
  {
    if (!result.rays.empty())
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
        
        result.rays.clear();
        result.rays.push_back(OptimizationOracle::Result::Ray(
          Vector(_coefficients, _integrality.size(), true)));
      }
      else
      {
        std::stringstream ss;
        ss << "Error in MakeRationalSolver::solve: Unbounded case yields SoPlex status " << status
          << '.'; 
        throw std::runtime_error(ss.str());
      }
    }

    std::vector<OptimizationOracle::Result::Point> points;
    for (const auto& point : result.points)
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
      for (std::size_t i = 0; i < point.vector.size(); ++i)
      {
        std::size_t coordinate = point.vector.coordinate(i);
        if (_integrality[coordinate])
        {
          mpq_set_d(_coefficients[0], round(point.vector.real(i)));
          _spx.changeBoundsRational(coordinate, _coefficients[0], _coefficients[0]);
        }
      }

      // Solve the LP.
      soplex::SPxSolver::Status status = _spx.solve();

      if (status == soplex::SPxSolver::OPTIMAL)
      {
        _spx.getPrimalRational(_coefficients, _integrality.size());
        soplex::Rational objVal = _spx.objValueRational();
        points.push_back(OptimizationOracle::Result::Point(
          Vector(_coefficients, _integrality.size(), true), Value(objVal.getMpqPtr())));
      }
      else if (status == soplex::SPxSolver::UNBOUNDED)
      {
        // If we have a ray already, then we ignore this point.
        if (!result.rays.empty())
          continue;

        // Otherwise, create a ray.
         _spx.getPrimalRayRational(_coefficients, _integrality.size());

        result.rays.push_back(Vector(_coefficients, _integrality.size(), true));
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
      result.dualBound = plusInfinity();
    else if (result.points.empty())
      result.dualBound = minusInfinity();
    else
    {
      result.dualBound.rational = result.dualBound.real;
      for (const auto& point : result.points)
      {
        if (point.objectiveValue > result.dualBound.rational)
          result.dualBound = point.objectiveValue;
      }
    }
  }

} /* namespace ipo */
