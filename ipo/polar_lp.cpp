#include "polar_lp.h"

#include <algorithm>
#include <random>
#include <vector>

using namespace soplex;

//#define CUT_AGING_BASIC

namespace ipo {

  PolarLP::PolarLP(OracleBase* oracle, double initialPenalty, int maxAge) 
    : _oracle(oracle), _n(oracle->space().dimension()), _d(oracle->space().dimension() + 1), 
    _offsetLower(oracle->space().dimension() + 1), _offsetUpper(2 * (oracle->space().dimension() + 1)), _stabilizing(false),
    _maxAge(maxAge), _initialPenalty(initialPenalty), _stabPenalty(0), _lastMainObjective(0.0), _lastPenaltyCosts(0.0)
  {

    // TODO: have a temporary maxAge variable during each run which is increased if stuck.

    /// Main LP

//     std::cerr << "Creating new SoPlex instance mainLP." << std::endl;
    _mainLP = new SoPlex;
    _mainLP->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    _mainLP->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    _mainLP->setRealParam(SoPlex::FEASTOL, 0.0);
    _mainLP->setBoolParam(SoPlex::RATREC, true);
    _mainLP->setBoolParam(SoPlex::RATFAC, true);
    _mainLP->setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    _mainLP->setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    _mainLP->setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);

    LPColSetRational cols(_d);
    DSVectorRational vector;
    for (std::size_t c = 0; c < _d; ++c)
      cols.add(Rational(0), -infinity, vector, infinity);
    _mainLP->addColsRational(cols);

    /// Stabilization LP.

//     std::cerr << "Creating new SoPlex instance stabLP." << std::endl;
    _stabLP = new SoPlex;
    _stabLP->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_REAL);
    _stabLP->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_ONLYREAL);
    _stabLP->setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    _stabLP->setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    _stabLP->setIntParam(SoPlex::SIMPLIFIER, SoPlex::SIMPLIFIER_OFF);

    _offsetStabilizationRows = 2 * _d;
    LPColSetReal stabCols(3 * _d);
    DSVectorReal stabVector;
    for (std::size_t c = 0; c < _d; ++c)
      stabCols.add(0.0, -infinity, stabVector, infinity);
    for (std::size_t c = 0; c < 2 * _d; ++c)
      stabCols.add(0.0, 0.0, stabVector, infinity);
    _stabLP->addColsReal(stabCols);

    _stabColInfos.resize(_d);

    /// Create rows.

    LPRowSetReal stabRows(2 * _d);
    for (std::size_t c = 0; c < _d; ++c)
    {
      _stabColInfos[c].lowerRow = stabRows.num();
      stabVector.clear();
      stabVector.add(c, 1);
      stabVector.add(_offsetLower + c, 1.0);
      stabRows.add(0, stabVector, infinity);
      RowInfo ri = { 's', SparseVector(), -1 };
      _stabRowInfos.push_back(ri);
    }
    for (std::size_t c = 0; c < _d; ++c)
    {
      _stabColInfos[c].upperRow = stabRows.num();
      stabVector.clear();
      stabVector.add(c, 1.0);
      stabVector.add(_offsetUpper + c, -1.0);
      stabRows.add(-infinity, stabVector, 0);
      RowInfo ri = { 's', SparseVector(), -1 };
      _stabRowInfos.push_back(ri);
    }
    _stabLP->addRowsReal(stabRows);
  }

  PolarLP::~PolarLP()
  {
    delete _stabLP;
    delete _mainLP;
  }

  void PolarLP::setBounds(std::size_t column, const Rational& lower, const Rational& upper)
  {
    assert(column <= n());

    _mainLP->changeBoundsRational(column, lower, upper);
    _stabLP->changeBoundsReal(column, double(lower), double(upper));
  }

  std::size_t PolarLP::addConstraint(const Rational& lhs, const SVectorRational& row, const Rational& rhs)
  {
    for (int p = row.size() - 1; p >= 0; --p)
      assert(row.index(p) <= _n);

    _mainLP->addRowRational(LPRowRational(lhs, row, rhs));
    RowInfo ri = { 'C', SparseVector(), -1 };
    _rowInfos.push_back(ri);
    _constraintsToRows.push_back(_mainLP->numRowsRational() - 1);

    DSVectorReal realRow;
    for (int p = row.size() - 1; p >= 0; --p)
      realRow.add(row.index(p), double(row.value(p)));
    _stabLP->addRowReal(LPRowReal(double(lhs), realRow, double(rhs)));
    _stabRowInfos.push_back(ri);

    return _constraintsToRows.size() - 1;
  }

  void PolarLP::updateConstraint(std::size_t index, const LPRowRational& row)
  {
    for (int p = row.rowVector().size() - 1; p >= 0; --p)
      assert(row.rowVector().index(p) <= _n);

    assert(index < _constraintsToRows.size());
    std::size_t rowIndex = _constraintsToRows[index];
    assert(rowIndex < _mainLP->numRowsRational());

    _mainLP->changeRowRational(rowIndex, row);

    DSVectorReal realVector;
    for (int p = row.rowVector().size() - 1; p >= 0; --p)
      realVector.add(row.rowVector().index(p), double(row.rowVector().value(p)));
    _stabLP->changeRowReal(_offsetStabilizationRows + rowIndex,
        LPRowReal(double(row.lhs()), realVector, double(row.rhs())));
  }

  void PolarLP::updateConstraint(std::size_t index, const soplex::Rational& lhs, const
soplex::SVectorRational& normal,
      const soplex::Rational& rhs)
  {
    updateConstraint(index, LPRowRational(lhs, normal, rhs));
  }

  void PolarLP::updateConstraint(std::size_t index, const soplex::Rational& lhs, const
soplex::VectorRational& normal,
      const soplex::Rational& rhs)
  {
    DSVectorRational sparseNormal;
    for (std::size_t i = 0; i < normal.dim(); ++i)
    {
      if (normal[i] != 0)
        sparseNormal.add(i, normal[i]);
    }
    updateConstraint(index, LPRowRational(lhs, sparseNormal, rhs));
  }

  void PolarLP::updateObjective(const DenseVector& objective)
  {
    assert(objective.dim() == _d);

    for (std::size_t c = 0; c < _d; ++c)
    {
      if (_mainLP->objRational(c) == objective[c])
        continue;

      _mainLP->changeObjRational(c, objective[c]);
      _stabLP->changeObjReal(c, double(objective[c]));
    }
  }

  std::size_t PolarLP::addPointConstraint(const SparseVector& point)
  {
    DSVectorRational exactVector;
    for (std::size_t p = 0; p < point.size(); ++p)
      exactVector.add(point.index(p), point.value(p));
    exactVector.add(_n, Rational(-1));
    _mainLP->addRowRational(LPRowRational(-infinity, exactVector, Rational(0)));

    RowInfo ri = { 'P', point, -1 };
    _rowInfos.push_back(ri);
    _constraintsToRows.push_back(_mainLP->numRowsRational() - 1);

    DSVectorReal approximateVector;
    for (std::size_t p = 0; p < point.size(); ++p)
      approximateVector.add(point.index(p), point.approximation(p));
    approximateVector.add(_n, -1.0);
    _stabLP->addRowReal(LPRowReal(-infinity, approximateVector, 0.0));
    _stabRowInfos.push_back(ri);

    return _constraintsToRows.size() - 1;
  }

  std::size_t PolarLP::addRayConstraint(const SparseVector& ray)
  {
    DSVectorRational exactVector;
    for (std::size_t p = 0; p < ray.size(); ++p)
      exactVector.add(ray.index(p), ray.value(p));
    _mainLP->addRowRational(LPRowRational(-infinity, exactVector, Rational(0)));

    RowInfo ri = { 'R', ray, -1 };
    _rowInfos.push_back(ri);
    _constraintsToRows.push_back(_mainLP->numRowsRational() - 1);

    DSVectorReal approximateVector;
    for (std::size_t p = 0; p < ray.size(); ++p)
      approximateVector.add(ray.index(p), ray.approximation(p));
    _stabLP->addRowReal(LPRowReal(-infinity, approximateVector, 0.0));
    _stabRowInfos.push_back(ri);

    return _constraintsToRows.size() - 1;
  }

  void PolarLP::setBasis(const Basis& basis)
  {
    assert(basis.columnStatus.size() == _n + 1);
    assert(basis.constraintStatus.size() == _constraintsToRows.size());

    std::vector<SPxSolver::VarStatus> rowStatus(_mainLP->numRowsRational());
    std::vector<SPxSolver::VarStatus> colStatus(_mainLP->numColsRational());
    for (std::size_t c = 0; c <= _n; ++c)
      colStatus[c] = basis.columnStatus[c];
    for (std::size_t c = _n + 1; c < _mainLP->numColsRational(); ++c)
      colStatus[c] = SPxSolver::ON_LOWER;
    for (std::size_t i = 0; i < _constraintsToRows.size(); ++i)
      rowStatus[_constraintsToRows[i]] = basis.constraintStatus[i];
    for (std::size_t r = 0; r < _mainLP->numRowsRational(); ++r)
    {
      if (_rowInfos[r].type == 'C' || _rowInfos[r].type == 'P' || _rowInfos[r].type == 'R')
        continue;
      bool tight = false;
      if (_rowInfos[r].type == 'p')
      {
        for (std::size_t i = 0; i < basis.tightPoints.size(); ++i)
          tight |= (basis.tightPoints[i] == _rowInfos[r].vector);
      }
      else if (_rowInfos[r].type == 'r')
      {
        for (std::size_t i = 0; i < basis.tightRays.size(); ++i)
          tight |= (basis.tightRays[i] == _rowInfos[r].vector);
      }
      else
        assert(_rowInfos[r].type == 'b');
      rowStatus[r] = tight ? SPxSolver::ON_UPPER : SPxSolver::BASIC;
    }
    _mainLP->setBasis(&rowStatus[0], &colStatus[0]);
  }

  void PolarLP::getBasis(Basis& basis)
  {
    assert(_mainLP->hasBasis());

    std::vector<SPxSolver::VarStatus> rowStatus(_mainLP->numRowsRational());
    std::vector<SPxSolver::VarStatus> colStatus(_mainLP->numColsRational());

    _mainLP->getBasis(&rowStatus[0], &colStatus[0]);
    basis.columnStatus.resize(_n + 1);
    for (std::size_t c = 0; c <= _n; ++c)
      basis.columnStatus[c] = colStatus[c];
    basis.constraintStatus.resize(_constraintsToRows.size());
    for (std::size_t i = 0; i < _constraintsToRows.size(); ++i)
      basis.constraintStatus[i] = rowStatus[_constraintsToRows[i]];
    basis.tightPoints.clear();
    basis.tightRays.clear();
    for (std::size_t r = 0; r < _mainLP->numRowsRational(); ++r)
    {
      if (rowStatus[r] != SPxSolver::ON_UPPER)
        continue;
      if (_rowInfos[r].type == 'p' || _rowInfos[r].type == 'P')
        basis.tightPoints.push_back(_rowInfos[r].vector);
      else if (_rowInfos[r].type == 'r' || _rowInfos[r].type == 'R')
        basis.tightRays.push_back(_rowInfos[r].vector);
    }
  }

  Rational PolarLP::getObjectiveValue()
  {
    return _mainLP->objValueRational();
  }

  void PolarLP::getPrimalSolution(VectorRational& solution)
  {
    _currentPrimalSolution.reDim(_mainLP->numColsRational());
    if (!_mainLP->getPrimalRational(_currentPrimalSolution))
      throw std::runtime_error("PolarLP: No primal solution available.");
    for (std::size_t v = 0; v <= _n; ++v)
      solution[v] = _currentPrimalSolution[v];
  }

  void PolarLP::getPrimalSolution(DSVectorRational& solution)
  {
    _currentPrimalSolution.reDim(_mainLP->numColsRational());
    assert(_mainLP->getPrimalRational(_currentPrimalSolution));
    assert(solution.size() == 0);
    for (std::size_t v = 0; v <= _n; ++v)
    {
      if (_currentPrimalSolution[v] != 0)
        solution.add(v, _currentPrimalSolution[v]);
    }
  }
  
  SparseVector PolarLP::getPrimalSolution()
  {
    _currentPrimalSolution.reDim(_mainLP->numColsRational());
    assert(_mainLP->getPrimalRational(_currentPrimalSolution));

    std::size_t size = 0;
    for (std::size_t v = 0; v <= _n; ++v)
    {
      if (_currentPrimalSolution[v] != 0)
        ++size;
    }
    SparseVector result(size);
    for (std::size_t v = 0; v <= _n; ++v)
    {
      if (_currentPrimalSolution[v] != 0)
        result.add(v, _currentPrimalSolution[v]);
    }
  }

  void PolarLP::stabilizedPresolve()
  {
    const double APPROX_VIOLATION = 1.0e-5;

    /// Remove dynamic rows of exact LP.

    std::vector<int> rowPermutation;
    rowPermutation.resize(_mainLP->numRowsRational());
    for (std::size_t i = 0; i < _mainLP->numRowsRational(); ++i)
    {
      rowPermutation[i] = (_rowInfos[i].type == 'p' || _rowInfos[i].type == 'r') ? -1 : 0;
    }
    _mainLP->removeRowsRational(&rowPermutation[0]);
    _rowInfos.resize(_mainLP->numRowsRational());

    /// Reset stabilization procedures.

    int localMaxAge = _maxAge;
    int numLastCacheRounds = 0; // Number of subsequent rows without really new solutions.
    _stabilizing = true;
    _lastMainObjective = 0.0;
    _lastPenaltyCosts = 0.0;
    _stabPenalty = _initialPenalty;
    for (std::size_t c = 0; c < _d; ++c)
    {
      for (std::size_t h = 0; h < IPO_POLAR_STABIL_HIST_SIZE; ++h)
        _stabColInfos[c].history[h] = 0.0;
      _stabColInfos[c].trustLower = 0.0;
      _stabColInfos[c].trustUpper = 0.0;

      _stabLP->changeLhsReal(_stabColInfos[c].lowerRow, 0);
      _stabLP->changeRhsReal(_stabColInfos[c].upperRow, 0);
      _stabLP->changeObjReal(_offsetLower + c, -_stabPenalty);
      _stabLP->changeObjReal(_offsetUpper + c, -_stabPenalty);
    }

    std::size_t numRows;
    DVectorReal primalSolution;
    DVectorReal inequalityApproxDenseNormal;
    inequalityApproxDenseNormal.reDim(_n);
    double inequalityApproxRhs;
    DVectorRational inequalityRoundedDenseNormal;
    inequalityRoundedDenseNormal.reDim(_n);
    Rational inequalityRoundedRhs;
    DVectorReal dualSolution;
    primalSolution.reDim(3 * _d, false);
    std::vector<SPxSolver::VarStatus> rowStatus;
    std::vector<SPxSolver::VarStatus> columnStatus(_stabLP->numColsReal());
    bool firstRound = true; // Indicator telling us to remove all dynamic rows in the first round.
    double lastObjective = 0.0;
    while (true)
    {
      assert(_stabLP->numRowsReal() == _stabRowInfos.size());

      /// Perform cut aging

      numRows = _stabLP->numRowsReal();
      rowPermutation.resize(numRows);
      for (std::size_t i = 0; i < numRows; ++i)
      {
        rowPermutation[i] = 0;
        if (_stabRowInfos[i].type != 'p' && _stabRowInfos[i].type != 'r')
          continue;
        if (_stabRowInfos[i].age >= 0)
        {
#ifdef CUT_AGING_BASIC
          if (i < rowStatus.size() && rowStatus[i] == SPxSolver::ON_UPPER && !firstRound)
#else
          if (i < dualSolution.dim() && dualSolution[i] != 0)
#endif
          {
            _stabRowInfos[i].age = 0;
          }
          else
          {
            _stabRowInfos[i].age++;
            if (_stabRowInfos[i].age >= localMaxAge || firstRound)
            {
              rowPermutation[i] = -1;
            }
          }
        }
      }

      _stabLP->removeRowsReal(&rowPermutation[0]);
      for (std::size_t i = 0; i < numRows; ++i)
      {
        if (rowPermutation[i] < 0)
          continue;
        std::size_t newRow = rowPermutation[i];
        if (newRow >= 0 && newRow != i)
          _stabRowInfos[newRow] = _stabRowInfos[i];
      }
      _stabRowInfos.resize(_stabLP->numRowsReal());
      numRows = _stabLP->numRowsReal();

      /// Solve LP and evaluate result.

      onBeforeSolve(true);
      SPxSolver::Status status = _stabLP->solve();
      if (status != SPxSolver::OPTIMAL)
      {
        _stabLP->writeFileReal("stab.lp");

        std::stringstream ss;
        ss << "PolarLP: Stabilization LP is " << status <<
          ", although provably feasible and bounded.";
        throw std::runtime_error(ss.str());
      }

      rowStatus.resize(numRows);
      _stabLP->getBasis(&rowStatus[0], &columnStatus[0]);
      firstRound = false;
      if (_stabLP->objValueReal() < lastObjective)
        numLastCacheRounds = 0;
      lastObjective = _stabLP->objValueReal();

      if (!_stabLP->getPrimalReal(primalSolution))
        throw std::runtime_error("Stabilization: No primal solution available.");
      for (std::size_t v = 0; v < _n; ++v)
        inequalityApproxDenseNormal[v] = primalSolution[v];
      inequalityApproxRhs = primalSolution[_n];
      dualSolution.reDim(numRows, false);
      if (!_stabLP->getDualReal(dualSolution))
        throw std::runtime_error("Stabilization: No dual solution available.");

      // TODO: SOLUTION DUMP
//      std::cout << "\nSOLUTION:" << std::setw(5) << std::setprecision(3) << iteration << " " <<
// std::setw(6)
//          << primalSolution[10] << " " << std::setw(6) << primalSolution[11] << std::endl;
//      ++iteration;

      _lastMainObjective = 0;
      for (std::size_t c = 0; c < _d; ++c)
        _lastMainObjective += primalSolution[c] * _stabLP->objReal(c);
      _lastPenaltyCosts = _lastMainObjective - _stabLP->objValueReal();
      if (_lastPenaltyCosts < 0)
        _lastPenaltyCosts = 0;
      onAfterSolve(true);

      /// Stabilization.

      if (_stabLP->objValueReal() <= _stabLP->realParam(SoPlex::OPTTOL) && _stabPenalty > 1.0e-6)
      {
        _stabPenalty *= 0.5;
        for (std::size_t c = 0; c < _d; ++c)
        {
          _stabLP->changeObjReal(_offsetLower + c, -_stabPenalty);
          _stabLP->changeObjReal(_offsetUpper + c, -_stabPenalty);
        }
        onPenaltyDecrease();
        numLastCacheRounds = 0;
        continue;
      }

      /// Call oracle.

      for (std::size_t v = 0; v < _n; ++v)
      {
        if (fabs(inequalityApproxDenseNormal[v]) > 1.0e-7)
          inequalityRoundedDenseNormal[v] = inequalityApproxDenseNormal[v];
        else
          inequalityRoundedDenseNormal[v] = 0.0;
      }
      inequalityRoundedRhs = inequalityApproxRhs + APPROX_VIOLATION;

      onBeforeOracleCall();
      _oracle->maximize(_result, inequalityRoundedDenseNormal, ObjectiveBound(inequalityRoundedRhs,
        true));
      onAfterOracleCall(_result.isFeasible(), _result.points.size(), _result.directions.size(),
        false);

      if (_result.isUnbounded())
      {
        for (std::size_t i = 0; i < _result.directions.size(); ++i)
        {
          onBeforeAddRay();
          addRayRow(_result.directions[i].vector, true);
          onAfterAddRay();
        }
        continue;
      }
      else if (_result.isFeasible())
      {
        if (_result.points.front().objectiveValue > inequalityRoundedRhs)
        {
          for (std::size_t i = 0; i < _result.points.size(); ++i)
          {
            onBeforeAddPoint();
            addPointRow(_result.points[i].vector, true);
            onAfterAddPoint();
          }
          continue;
        }
      }
      else
      {
        throw std::runtime_error("Polar LP: Oracle claims infeasible.");
      }

      bool penalizing = _stabPenalty > 1.0e-6;

      /// If we are feasible, then quickly decrease penalty, if not already low.

      if (penalizing)
      {
        _stabPenalty *= 0.125;
        for (std::size_t c = 0; c < _d; ++c)
        {
          _stabLP->changeObjReal(_offsetLower + c, -_stabPenalty);
          _stabLP->changeObjReal(_offsetUpper + c, -_stabPenalty);
        }
        onPenaltyDecrease();
        continue;
      }
      else
      {
        /// If stabilization is complete, then add rows of certifying points and rays to exact LP.

        for (std::size_t r = 0; r < numRows; ++r)
        {
          if (r < rowStatus.size() && rowStatus[r] == SPxSolver::ON_UPPER)
          {
            if (_stabRowInfos[r].type == 'p')
              addPointRow(_stabRowInfos[r].vector, false);
            if (_stabRowInfos[r].type == 'r')
              addRayRow(_stabRowInfos[r].vector, false);
          }
        }

        break;
      }
    }

    _stabilizing = false;
  }

  void PolarLP::optimize(bool perturbeObjective)
  {
    _lastMainObjective = 0.0;
    _lastPenaltyCosts = 0.0;
    _stabPenalty = 0.0;

    int localMaxAge = _maxAge;
    DSVectorRational inequalityExactSparseNormal;
    DVectorRational inequalityExactDenseNormal;
    DVectorReal inequalityApproxDenseNormal;
    DVectorRational dualSolution;
    Rational inequalityRhs;
    std::vector<SPxSolver::VarStatus> rowStatus;
    std::vector<SPxSolver::VarStatus> columnStatus(_mainLP->numColsRational());
    std::vector<int> rowPermutation;

    // Extract objective and compute perturbed one.

    DVectorRational objectiveVector, perturbedVector;
    if (perturbeObjective)
    {
      const Rational PERTURBATION_DENOMINATOR = 1024*1024*1024;

      std::default_random_engine generator;
      std::uniform_int_distribution<int> distribution(1,1024);

      objectiveVector.reDim(numColumns());
      _mainLP->getObjRational(objectiveVector);
      perturbedVector = objectiveVector;
      LPRowRational row;
      for (std::size_t r = 0; r < numRows(); ++r)
      {
        _mainLP->getRowRational(r, row);
        int lhsLambda = (row.lhs() > -infinity) ? distribution(generator) : 0;
        int rhsLambda = (row.rhs() < infinity) ? distribution(generator) : 0;
        if (lhsLambda == rhsLambda)
          continue;
        Rational lambda = (rhsLambda - lhsLambda) / PERTURBATION_DENOMINATOR;
        const SVectorRational& rowVector = row.rowVector();
        for (int p = rowVector.size() - 1; p >= 0; --p)
        {
          perturbedVector[rowVector.index(p)] += lambda * rowVector.value(p);
        }
      }
      _mainLP->changeObjRational(perturbedVector);
    }

    // Main loop.

    while (true)
    {
      /// Perform cut aging

      std::size_t numRows = _mainLP->numRowsRational();
      rowStatus.resize(numRows);
      rowPermutation.resize(numRows);
      _mainLP->getBasis(&rowStatus[0], &columnStatus[0]);
      std::size_t lastStatic = 0;
      for (std::size_t i = 0; i < numRows; ++i)
      {
        rowPermutation[i] = 0;
        if (_rowInfos[i].type != 'p' && _rowInfos[i].type != 'r')
        {
          lastStatic = i;
          continue;
        }

        if (_rowInfos[i].age >= 0)
        {
#ifdef CUT_AGING_BASIC
          if (rowStatus[i] == SPxSolver::ON_UPPER)
#else
          if (i < dualSolution.dim() && dualSolution[i] != 0)
#endif
          {
            _rowInfos[i].age = 0;
          }
          else
          {
            _rowInfos[i].age++;
            if (_rowInfos[i].age >= localMaxAge)
              rowPermutation[i] = -1;
          }
        }
      }

      _mainLP->removeRowsRational(&rowPermutation[0]);
      for (std::size_t i = 0; i < numRows; ++i)
      {
        if (rowPermutation[i] < 0)
          continue;
        std::size_t newRow = rowPermutation[i];
        if (newRow >= 0)
          _rowInfos[newRow] = _rowInfos[i];
      }
      _rowInfos.resize(_mainLP->numRowsRational());
      numRows = _mainLP->numRowsRational();

      /// Solve LP and evaluate.

      onBeforeSolve(false);
      SPxSolver::Status status = _mainLP->solve();
      if (status != SPxSolver::OPTIMAL)
      {
        std::stringstream ss;
        ss << "Polar LP: Unexpected status " << status << ".";
        throw std::runtime_error(ss.str());
      }

      rowStatus.resize(_mainLP->numRowsRational());
      _mainLP->getBasis(&rowStatus[0], &columnStatus[0]);
//       int numBasic = 0;
//       for (int c = 0; c < _mainLP->numColsRational(); ++c)
//       {
//         if (columnStatus[c] == SPxSolver::BASIC)
//           numBasic++;
//         std::cout << "Col " << c << " is " << columnStatus[c] << std::endl;
//       }
//       for (int r = 0; r < _mainLP->numRowsRational(); ++r)
//       {
//         if (rowStatus[r] == SPxSolver::ON_UPPER)
//         {
//           std::cerr << "Row " << r << " is on upper." << std::endl;
//         }
//         else
//         {
//           std::cerr << "Row " << r << " is not on upper." << std::endl;
//         }
//       }
//       std::cout << "#basic cols: " << numBasic << std::endl;

      dualSolution.reDim(numRows, false);
      if (!_mainLP->getDualRational(dualSolution))
        throw std::runtime_error("Optimization: No dual solution available.");

      double obj = double(_mainLP->objValueRational());
      _lastMainObjective = obj;
      onAfterSolve(false);

      // Extract solution vector.

      _currentPrimalSolution.reDim(_mainLP->numColsRational());
      _mainLP->getPrimalRational(_currentPrimalSolution);
      inequalityExactSparseNormal.clear();
      inequalityExactDenseNormal.reDim(n());
      inequalityApproxDenseNormal.reDim(n());
      for (std::size_t v = 0; v < _n; ++v)
      {
        inequalityExactDenseNormal[v] = _currentPrimalSolution[v];
        inequalityApproxDenseNormal[v] = double(_currentPrimalSolution[v]);
        if (_currentPrimalSolution[v] != 0)
          inequalityExactSparseNormal.add(v, _currentPrimalSolution[v]);
      }
      inequalityRhs = _currentPrimalSolution[_n];

      onBeforeOracleCall();
      _oracle->maximize(_result, inequalityExactDenseNormal, ObjectiveBound(inequalityRhs, true));
      onAfterOracleCall(_result.isFeasible(), _result.points.size(),
        _result.directions.size(), false);

      if (_result.isUnbounded())
      {
        for (std::size_t i = 0; i < _result.directions.size(); ++i)
        {
          onBeforeAddRay();
          addRayRow(_result.directions[i].vector, false);
          onAfterAddRay();
        }
        continue;
      }
      else if (_result.isFeasible())
      {
        if (_result.points.front().objectiveValue > inequalityRhs)
        {
          for (std::size_t i = 0; i < _result.points.size(); ++i)
          {
            onBeforeAddPoint();
            addPointRow(_result.points[i].vector, false);
            onAfterAddPoint();
          }
          continue;
        }
      }
      else
      {
        throw std::runtime_error("Polar LP: Oracle claims infeasible.");
      }

      if (perturbeObjective)
      {
        _mainLP->changeObjRational(objectiveVector);
        perturbeObjective = false;
        SPxSolver::Status status = _mainLP->solve();
        if (status != SPxSolver::OPTIMAL)
          throw std::runtime_error("PolarLP: Main LP is not optimal.");

        rowStatus.resize(_mainLP->numRowsRational());
        _mainLP->getBasis(&rowStatus[0], &columnStatus[0]);
      }

      break;
    }
  }

  void PolarLP::reoptimizePerturbed()
  {
    const Rational PERTURBATION_DENOMINATOR = 1024*1024*1024;

    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution(1,1024);

//     Rational objectiveValue = _mainLP->objValueRational();
//     std::cout << "Objective = " << objectiveValue << std::endl;

    DVectorRational objectiveVector(numColumns());
    _mainLP->getObjRational(objectiveVector);

    DVectorRational perturbedVector(numColumns());
    perturbedVector = objectiveVector;
    LPRowRational row;
    for (std::size_t r = 0; r < numRows(); ++r)
    {
      _mainLP->getRowRational(r, row);
      int lhsLambda = (row.lhs() > -infinity) ? distribution(generator) : 0;
      int rhsLambda = (row.rhs() < infinity) ? distribution(generator) : 0;
      if (lhsLambda == rhsLambda)
        continue;
      Rational lambda = (rhsLambda - lhsLambda) / PERTURBATION_DENOMINATOR;
      const SVectorRational& rowVector = row.rowVector();
      for (int p = rowVector.size() - 1; p >= 0; --p)
      {
        perturbedVector[rowVector.index(p)] += lambda * rowVector.value(p);
      }
    }
//     for (std::size_t c = 0; c < numColumns(); ++c)
//     {
//       std::cout << "Obj #" << c << ": " << objectiveVector[c] << " -> " << perturbedVector[c] <<
// std::endl;
//     }

    _mainLP->changeObjRational(perturbedVector);

    SPxSolver::Status status = _mainLP->solve();
    if (status != SPxSolver::OPTIMAL)
    {
      std::stringstream ss;
      ss << "Polar LP (perturbed): Unexpected status " << status << ".";
      throw std::runtime_error(ss.str());
    }


//     std::vector<SPxSolver::VarStatus> rowStatus(_mainLP->numRowsRational());
//     std::vector<SPxSolver::VarStatus> columnStatus(_mainLP->numColsRational());
//     _mainLP->getBasis(&rowStatus[0], &columnStatus[0]);
//     int numBasic = 0;
//     for (int c = 0; c < _mainLP->numColsRational(); ++c)
//     {
//       if (columnStatus[c] == SPxSolver::BASIC)
//         numBasic++;
//       std::cout << "Perturbed: Col " << c << " is " << columnStatus[c] << std::endl;
//     }
//     std::cout << "Perturbed: #basic cols: " << numBasic << std::endl;

    _mainLP->changeObjRational(objectiveVector);
    status = _mainLP->solve();
    if (status != SPxSolver::OPTIMAL)
    {
      std::stringstream ss;
      ss << "Polar LP (unperturbed): Unexpected status " << status << ".";
      throw std::runtime_error(ss.str());
    }

    _currentPrimalSolution.reDim(_mainLP->numColsRational());
    _mainLP->getPrimalRational(_currentPrimalSolution);
    DSVectorRational sparsePrimal;
    sparsePrimal = _currentPrimalSolution;

    std::cout << _mainLP->objValueRational() << std::endl;
    std::cout << sparsePrimal << std::endl;
    std::cout << std::endl;

//     objectiveValue = _mainLP->objValueRational();
//     std::cout << "Objective = " << objectiveValue << std::endl;

/*
    _mainLP->getBasis(&rowStatus[0], &columnStatus[0]);
    numBasic = 0;
    for (int c = 0; c < _mainLP->numColsRational(); ++c)
    {
      if (columnStatus[c] == SPxSolver::BASIC)
        numBasic++;
      std::cout << "Reset: Col " << c << " is " << columnStatus[c] << std::endl;
    }
    std::cout << "Reset: #basic cols: " << numBasic << std::endl;*/
  }

  void PolarLP::onBeforeSolve(bool stabilizing)
  {

  }

  void PolarLP::onAfterSolve(bool stabilizing)
  {

  }

  void PolarLP::onPenaltyDecrease()
  {

  }

  void PolarLP::onBeforeCache()
  {

  }

  void PolarLP::onAfterCache(std::size_t numPoints, std::size_t numRays)
  {

  }

  void PolarLP::onBeforeOracleCall()
  {

  }

  void PolarLP::onAfterOracleCall(bool feasible, std::size_t numPoints, std::size_t numRays,
      bool lastIteration)
  {

  }

  void PolarLP::onBeforeAddPoint()
  {

  }

  void PolarLP::onAfterAddPoint()
  {

  }

  void PolarLP::onBeforeAddRay()
  {

  }

  void PolarLP::onAfterAddRay()
  {

  }

  void PolarLP::addPointRow(SparseVector& point, bool stabLP)
  {
    DSVectorRational vector;
    for (std::size_t p = 0; p < point.size(); ++p)
      vector.add(point.index(p), point.value(p));
    vector.add(_n, Rational(-1));
    RowInfo ri = { 'p', point, 0 };

    if (stabLP)
    {
      DSVectorReal realVector;
      for (int p = vector.size() - 1; p >= 0; --p)
        realVector.add(vector.index(p), double(vector.value(p)));
      _stabLP->addRowReal(LPRowReal(-infinity, realVector, 0.0));
      _stabRowInfos.push_back(ri);
    }
    else
    {
      _mainLP->addRowRational(LPRowRational(-infinity, vector, Rational(0)));
      _rowInfos.push_back(ri);
    }
  }

  void PolarLP::addRayRow(SparseVector& ray, bool stabLP)
  {
    RowInfo ri = { 'r', ray, 0 };
    if (stabLP)
    {
      DSVectorReal vector;
      for (std::size_t p = 0; p < ray.size(); ++p)
        vector.add(ray.index(p), ray.approximation(p));
      _stabLP->addRowReal(LPRowReal(minusInfinity, vector, 0.0));
      _stabRowInfos.push_back(ri);
    }
    else
    {
      DSVectorRational vector;
      for (std::size_t p = 0; p < ray.size(); ++p)
        vector.add(ray.index(p), ray.value(p));
      _mainLP->addRowRational(LPRowRational(minusInfinity, vector, Rational(0)));
      _rowInfos.push_back(ri);
    }
  }


} /* namespace ipo */
