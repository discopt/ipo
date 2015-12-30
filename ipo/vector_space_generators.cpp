#include "vector_space_generators.h"

#include <limits>
#include <stdexcept>

using namespace soplex;

namespace ipo {

  VectorSpaceGenerators::VectorSpaceGenerators() : _zeroRhs(false)
  {
    _spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    _spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    _spx.setRealParam(SoPlex::FEASTOL, 0.0);
    _spx.setBoolParam(SoPlex::RATFAC, true);
    _spx.setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF); // TODO: There seems to be a bug in the simplifier.
    _spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
  }

  VectorSpaceGenerators::~VectorSpaceGenerators()
  {

  }

  void VectorSpaceGenerators::reset(std::size_t ambientDimension)
  {
    _spx.clearLPRational();
    LPRowSetRational rowSet(ambientDimension);
    DSVectorRational vector;
    for (std::size_t r = 0; r < ambientDimension; ++r)
      rowSet.add(Rational(0), vector, Rational(0));
    _spx.addRowsRational(rowSet);
    _lazyVectors.clear();
    _zeroRhs = true;
  }

  std::size_t VectorSpaceGenerators::add(const SVectorRational& vector)
  {
    addLazy(vector);
    flushLazy();
    return _spx.numColsRational() - 1;
  }

  std::size_t VectorSpaceGenerators::add(const SVectorRational& vector, bool dependent)
  {
    if (dependent)
      _spx.addColRational(LPColRational(Rational(0), vector, Rational(0), Rational(0)));
    else
      _spx.addColRational(LPColRational(Rational(0), vector, infinity, -infinity));
    return _spx.numColsRational() - 1;
  }

  std::size_t VectorSpaceGenerators::addLazy(const SVectorRational& vector)
  {
    _spx.addColRational(LPColRational(Rational(0), vector, infinity, -infinity));
    _lazyVectors.push_back(_spx.numColsRational() - 1);
    return _spx.numColsRational() - 1;
  }

  void VectorSpaceGenerators::flushLazy()
  {
    if (!_zeroRhs)
    {
      DVectorRational rhs;
      rhs.reDim(_spx.numRowsRational(), true);
      _spx.changeRangeRational(rhs, rhs);
      _zeroRhs = true;
    }

    _solution.reDim(_spx.numColsRational(), false);
    for (std::size_t i = 0; i < _lazyVectors.size(); ++i)
    {
      std::size_t column = _lazyVectors[i];
      _spx.changeBoundsRational(column, Rational(1), Rational(1));

      SPxSolver::Status status = _spx.solve();
      if (status == SPxSolver::INFEASIBLE)
      {
        _spx.changeBoundsRational(column, -infinity, +infinity);
        continue;
      }

      if (status != SPxSolver::OPTIMAL)
        throw std::runtime_error("VectorSpaceGenerator::flushLazy failed to solve system.");

      _spx.getPrimalRational(_solution);
      std::size_t densestColumn = std::numeric_limits<std::size_t>::max();
      std::size_t densestNumNonzeros = 0;
      for (std::size_t c = 0; c < _spx.numColsRational(); ++c)
      {
        if (_solution[c] == 0)
          continue;

        std::size_t numNonzeros = _spx.colVectorRational(i).size();
        if (densestColumn == std::numeric_limits<std::size_t>::max() || numNonzeros > densestNumNonzeros)
        {
          densestColumn = i;
          densestNumNonzeros = numNonzeros;
        }
      }
      assert(densestColumn < std::numeric_limits<std::size_t>::max());

      _spx.changeBoundsRational(densestColumn, Rational(0), Rational(0));
      if (densestColumn != column)
      {
        _spx.changeBoundsRational(column, -soplex::infinity, +soplex::infinity);
        --i;
      }

    }
    _lazyVectors.clear();
  }

  bool VectorSpaceGenerators::isDependent(std::size_t index)
  {
    LPColRational col;
    _spx.getColRational(index, col);
    return col.lower() == 0;
  }

  bool VectorSpaceGenerators::isDependent(const SVectorRational& vector)
  {
    DVectorRational denseVector;
    denseVector.reDim(_spx.numRowsRational(), true);
    denseVector.assign(vector);
    return isDependent(denseVector);
  }

  bool VectorSpaceGenerators::isDependent(const VectorRational& vector)
  {
    _zeroRhs = false;
    _spx.changeRangeRational(vector, vector);
    SPxSolver::Status status = _spx.solve();
    if (status == SPxSolver::OPTIMAL)
      return true;
    else if (status != SPxSolver::INFEASIBLE)
    {
      std::cerr << "\nFailure Status: " << status << std::endl;
      _spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_HIGH);
      _spx.solve();
      _spx.writeFileRational("VectorSpaceGenerator.lp", NULL, NULL, NULL);
      throw std::runtime_error("VectorSpaceGenerators::isDependent could not solve system.");
    }
    return false;
  }

} /* namespace polycomb */
