#include "affine_hull.h"

#include <algorithm>

#include "reconstruct.h"
#include "vector_space_generators.h"
#include "vectors.h"

using namespace soplex;

//#define RANDOM_COBASIS_INDEX
//#define EXACT_CACHE_SEARCH


namespace ipo {

  class Factorization
  {
  public:
    Factorization(std::size_t numColumns) :
        _numRows(0), _columnVectors(numColumns), _columnBasics(numColumns, std::numeric_limits<std::size_t>::max())
    {

    }

    ~Factorization()
    {

    }

    inline std::size_t numRows() const
    {
      return _numRows;
    }

    inline std::size_t numColumns() const
    {
      return _columnVectors.size();
    }

    inline bool isBasic(std::size_t columnIndex) const
    {
      return _columnBasics[columnIndex] != std::numeric_limits<std::size_t>::max();
    }

    void solveRightFactor(DVectorRational& vector) const
    {
      for (std::size_t b = _numRows; b > 0; b--)
      {
        std::size_t pivotRow = b - 1;
        const DSVectorRational& column = _factorUcolumnVectors[pivotRow];
        assert(column.index(column.size() - 1) == pivotRow);
        vector[pivotRow] /= column.value(column.size() - 1);

        for (int p = column.size() - 2; p >= 0; --p)
        {
          vector[column.index(p)] -= column.value(p) * vector[pivotRow];
        }
      }
    }

    void solveRightApproximateFactor(DVectorReal& vector) const
    {
      for (std::size_t b = _numRows; b > 0; b--)
      {
        std::size_t pivotRow = b - 1;
        const DSVectorReal& column = _approximateFactorUcolumnVectors[pivotRow];
        assert(column.index(column.size() - 1) == pivotRow);
        vector[pivotRow] /= column.value(column.size() - 1);

        for (int p = column.size() - 2; p >= 0; --p)
        {
          vector[column.index(p)] -= column.value(p) * vector[pivotRow];
        }
      }
    }

    void solveLeftFactor(DVectorRational& vector) const
    {
      for (std::size_t pivotRow = 0; pivotRow < numRows(); ++pivotRow)
      {
        const DSVectorRational& column = _factorLcolumnVectors[pivotRow];
        assert(column.index(0) == pivotRow);
        vector[pivotRow] /= column.value(0);

        for (int p = column.size() - 1; p >= 1; --p)
        {
          vector[column.index(p)] -= column.value(p) * vector[pivotRow];
        }
      }
    }

    void solveLeftApproximateFactor(DVectorReal& vector) const
    {
      for (std::size_t pivotRow = 0; pivotRow < numRows(); ++pivotRow)
      {
        const DSVectorReal& column = _approximateFactorLcolumnVectors[pivotRow];
        assert(column.index(0) == pivotRow);
        vector[pivotRow] /= column.value(0);

        for (int p = column.size() - 1; p >= 1; --p)
        {
          vector[column.index(p)] -= column.value(p) * vector[pivotRow];
        }
      }
    }

    void computeKernelVector(std::size_t columnIndex, soplex::DVectorRational& result) const
    {
      soplex::DVectorRational rhs;
      rhs.reDim(numRows(), true);
      const SVectorRational& col = _columnVectors[columnIndex];
      for (int p = col.size() - 1; p >= 0; --p)
        rhs[col.index(p)] = -col.value(p);

      solveLeftFactor(rhs);
      solveRightFactor(rhs);

      std::size_t size = 1;
      for (std::size_t b = 0; b < _basicColumns.size(); ++b)
      {
        if (rhs[b] != 0)
          ++size;
      }
      result.reDim(numColumns(), false);
      result.clear();
      for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        result[_basicColumns[b]] = rhs[b];
      result[columnIndex] = Rational(1);
    }

    void computeApproximateKernelVector(std::size_t columnIndex, soplex::DVectorReal& result, double epsilon) const
    {
      DVectorReal rhs;
      rhs.reDim(numRows(), true);
      const SVectorRational& col = _columnVectors[columnIndex];
      for (int p = col.size() - 1; p >= 0; --p)
        rhs[col.index(p)] = -(double) col.value(p);

      solveLeftApproximateFactor(rhs);
      solveRightApproximateFactor(rhs);

      result.reDim(numColumns());
      result.clear();
      for (std::size_t b = 0; b < _basicColumns.size(); ++b)
      {
        if (fabs(rhs[b]) > epsilon)
          result[_basicColumns[b]] = rhs[b];
      }
      result[columnIndex] = 1.0;
    }

    void dump()
    {
      std::cerr << "ACM: Matrix is " << numRows() << "x" << numColumns() << "\n";
      std::cerr << "ACM: Basis:";
      for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        std::cerr << " " << _basicColumns[b];
      std::cerr << std::endl;
      for (std::size_t c = 0; c < numColumns(); ++c)
      {
        std::cerr << "ACM: Column " << c << ":";
        for (int p = 0; p < _columnVectors[c].size(); ++p)
        {
          std::cerr << "{" << _columnVectors[c].index(p) << ":" << _columnVectors[c].value(p) << "}";
        }
        std::cerr << std::endl;
      }
      std::cerr << "ACM: L:\n";
      for (std::size_t b = 0; b < _basicColumns.size(); ++b)
      {
        std::cerr << "ACM: Column " << b << " -> " << _basicColumns[b] << ":";
        for (int p = 0; p < _factorLcolumnVectors[b].size(); ++p)
        {
          std::cerr << "{" << _factorLcolumnVectors[b].index(p) << ":" << _factorLcolumnVectors[b].value(p) << "}";
        }
        std::cerr << std::endl;
      }
      std::cerr << "ACM: U:\n";
      for (std::size_t b = 0; b < _basicColumns.size(); ++b)
      {
        std::cerr << "ACM: Column " << b << " -> " << _basicColumns[b] << ":";
        for (int p = 0; p < _factorUcolumnVectors[b].size(); ++p)
        {
          std::cerr << "{" << _factorUcolumnVectors[b].index(p) << ":" << _factorUcolumnVectors[b].value(p) << "}";
        }
        std::cerr << std::endl;
      }

    }

    void addRow(const Vector& row, const Rational& last, std::size_t newBasicColumnIndex)
    {
      assert(_columnBasics[newBasicColumnIndex] == std::numeric_limits<std::size_t>::max());

      // Compute new column which is aligned right to the right factor.

      DVectorRational rhs;
      rhs.reDim(numRows(), true);
      rhs.assign(_columnVectors[newBasicColumnIndex]);
      solveLeftFactor(rhs);
      _factorUcolumnVectors.push_back(DSVectorRational());
      _factorUcolumnVectors.back() = rhs;
      _factorUcolumnVectors.back().sort();
      _basicColumns.push_back(newBasicColumnIndex);
      _columnBasics[newBasicColumnIndex] = numRows();

      _approximateFactorUcolumnVectors.push_back(DSVectorReal());
      for (int p = 0; p < _factorUcolumnVectors.back().size(); ++p)
      {
        _approximateFactorUcolumnVectors.back().add(_factorUcolumnVectors.back().index(p),
            (double) _factorUcolumnVectors.back().value(p));
      }

      // Copy new row to worker which is a row vector indexed the ordered basis.

      DVectorRational worker;
      worker.reDim(numRows() + 1, true);
      for (int p = row.size() - 1; p >= 0; --p)
      {
        if (_columnBasics[row.index(p)] != std::numeric_limits<std::size_t>::max())
          worker[_columnBasics[row.index(p)]] = row.value(p);
      }
      if (last != 0 && _columnBasics[numColumns() - 1] != std::numeric_limits<std::size_t>::max())
        worker[_columnBasics[numColumns() - 1]] = last;

      std::vector<int> positions(_basicColumns.size(), 0);
      for (std::size_t b = 0; b < _basicColumns.size() - 1; ++b)
      {
        if (worker[b] == 0)
          continue;

        std::size_t pivotRowIndex = b;
        std::size_t pivotColumnIndex = _basicColumns[b];
        assert(_factorUcolumnVectors[b].index(_factorUcolumnVectors[b].size() - 1) == pivotRowIndex);
        const Rational& pivotElement = _factorUcolumnVectors[b].value(_factorUcolumnVectors[b].size() - 1);
        assert(pivotElement != 0);
        const Rational lambda = -worker[b] / pivotElement;

        for (std::size_t i = b + 1; i < _basicColumns.size(); ++i)
        {
          const SVectorRational& currentRightColumn = _factorUcolumnVectors[i];
          while (positions[i] < currentRightColumn.size() && currentRightColumn.index(positions[i]) < pivotRowIndex)
            positions[i]++;

          // Pivot row entry may be zero:

          if (positions[i] == currentRightColumn.size() || currentRightColumn.index(positions[i]) > pivotRowIndex)
            continue;

          assert(currentRightColumn.index(positions[i]) == pivotRowIndex);

          worker[i] += lambda * currentRightColumn.value(positions[i]);
        }

        /// Put -lambda into left factor.
        _factorLcolumnVectors[b].add(_numRows, -lambda);
        _approximateFactorLcolumnVectors[b].add(_numRows, -(double) lambda);
      }
      if (worker[numRows()] == 0)
      {
        throw std::runtime_error("Invalid call to AffineComputationMatrix::addRow: Basis matrix is singular!");
      }

      // Diagonal entries.

      _factorLcolumnVectors.push_back(DSVectorRational());
      _factorLcolumnVectors[numRows()].add(numRows(), Rational(1));
      _factorUcolumnVectors[numRows()].add(numRows(), worker[numRows()]);

      _approximateFactorLcolumnVectors.push_back(DSVectorReal());
      _approximateFactorLcolumnVectors[numRows()].add(numRows(), 1.0);
      _approximateFactorUcolumnVectors[numRows()].add(numRows(), (double) worker[numRows()]);

      // Physically add the row.

      for (int p = row.size() - 1; p >= 0; --p)
      {
        _columnVectors[row.index(p)].add(numRows(), row.value(p));
      }
      if (last != 0)
        _columnVectors[numColumns() - 1].add(numRows(), last);
      ++_numRows;
    }

    void getLastRow(DVectorRational& row) const
    {
      row.reDim(numColumns(), true);
      if (numRows() == 0)
        return;
      int rowIndex = numRows() - 1;
      for (std::size_t c = 0; c < numColumns(); ++c)
      {
        int p = _columnVectors[c].size() - 1;
        if (p >= 0 && _columnVectors[c].index(p) == rowIndex)
        {
          row[c] = _columnVectors[c].value(p);
        }
      }
    }

  protected:
    std::size_t _numRows;
    std::vector<DSVectorRational> _columnVectors;
    std::vector<std::size_t> _basicColumns;
    std::vector<std::size_t> _columnBasics;
    std::vector<DSVectorRational> _factorLcolumnVectors;
    std::vector<DSVectorRational> _factorUcolumnVectors;
    std::vector<DSVectorReal> _approximateFactorLcolumnVectors;
    std::vector<DSVectorReal> _approximateFactorUcolumnVectors;
  };

  struct NonbasicColumn
  {
    bool exactValid;
    soplex::DVectorRational exactDirection;
    Rational exactRhs;
    bool approximateValid;
    soplex::DVectorReal approximateDirection;
    std::size_t sparsity;
    bool definesEquation;
    bool avoid;

    NonbasicColumn()
      : exactDirection(0)
    {
      exactValid = false;
      approximateValid = false;
      definesEquation = false;
      avoid = false;
      sparsity = 0;
    }

    bool requiredBefore(const NonbasicColumn& other) const
    {
      if (definesEquation)
        return false;
      if (other.definesEquation)
        return true;

      if (avoid)
        return false;
      if (other.avoid)
        return true;

      return false;
    }

    bool preferredBeforeSparsity(const NonbasicColumn& other) const
    {
      assert(!this->requiredBefore(other));
      assert(!other.requiredBefore(*this));

      return sparsity < other.sparsity;
    }
  };

  AffineHullState::AffineHullState()
  {

  }

  AffineHullState::~AffineHullState()
  {

  }

  class XAffineHull : public AffineHullState
  {
  public:
    XAffineHull(std::vector<AffineHullHandler*>& handlers, const std::shared_ptr<OracleBase>& oracle,
      InnerDescription& resultInnerDescription, AffineOuterDescription& resultOuterDescription)
      : _handlers(handlers), _oracle(oracle), _innerDescription(resultInnerDescription), _equations(resultOuterDescription),
      _infeasible(false), _n(oracle->space().dimension()), _factorization(_n + 1), _objectiveValuePositive(soplex::infinity),
      _objectiveValueNegative(soplex::infinity),
      _verifyIndex(std::numeric_limits<std::size_t>::max()),
      _verifySuccess(true), _lastApproximateDirectionSolves(0), _lastExactDirectionSolves(0),
      _lastExactDirectionNonzeros(0), _lastExactDirectionBitsize(std::numeric_limits<std::size_t>::max())
    {
      _innerDescription.points.clear();
      _innerDescription.rays.clear();
      _equations.clear();
      _directionPositiveVector.reDim(_n, false);
      _directionNegativeVector.reDim(_n, false);
      _columns.resize(_n + 1);
      _columns[_n].avoid = true;
    }

    virtual ~XAffineHull()
    {

    }

    virtual const Space& space() const
    {
      return _oracle->space();
    }

    virtual std::size_t numPoints() const
    {
      return _innerDescription.points.size();
    }

    virtual std::size_t numRays() const
    {
      return _innerDescription.rays.size();
    }

    virtual std::size_t numEquations() const
    {
      return _equations.size();
    }

    virtual std::size_t numCandidateEquations() const
    {
      return _candidateEquations.size();
    }

    virtual bool infeasible() const
    {
      return _infeasible;
    }

    virtual std::size_t verifyIndex() const
    {
      return _verifyIndex;
    }

    virtual bool verifySuccess() const
    {
      return _verifySuccess;
    }

    virtual std::size_t directionApproximateSolves() const
    {
      return _lastApproximateDirectionSolves;
    }

    virtual std::size_t directionExactSolves() const
    {
      return _lastExactDirectionSolves;
    }

    virtual std::size_t directionNonzeros()
    {
      if (_lastExactDirectionNonzeros == std::numeric_limits<std::size_t>::max())
      {
        _lastExactDirectionNonzeros = 0;
        for (std::size_t v = 0; v < _n; ++v)
        {
          if (_directionPositiveVector[v] != 0)
            ++_lastExactDirectionNonzeros;
        }
      }
      return _lastExactDirectionNonzeros;
    }

    virtual std::size_t directionBitsize()
    {
      if (_lastExactDirectionBitsize == std::numeric_limits<std::size_t>::max())
      {
        _lastExactDirectionBitsize = 0;
        for (std::size_t v = 0; v < _n; ++v)
        {
          if (_directionPositiveVector[v] != 0)
            _lastExactDirectionBitsize += _directionPositiveVector[v].sizeInBase(2);
        }
      }
      return _lastExactDirectionBitsize;
    }

    virtual const soplex::VectorRational& directionVector()
    {
      return _directionPositiveVector;
    }

    virtual HeuristicLevel oracleMaxHeuristicLevel() const
    {
      return _oracleMaxHeuristicLevel;
    }

    virtual HeuristicLevel oracleMinHeuristicLevel() const
    {
      return _oracleMinHeuristicLevel;
    }

    virtual HeuristicLevel oracleResultHeuristicLevel() const
    {
      return _oracleResultHeuristicLevel;
    }

    virtual std::size_t oracleNumPoints() const
    {
      return _result.points.size();
    }

    virtual std::size_t oracleNumRays() const
    {
      return _result.rays.size();
    }

    virtual const InnerDescription& innerDescription() const
    {
      return _innerDescription;
    }

    virtual const AffineOuterDescription& outerDescription() const
    {
      return _equations;
    }

    virtual const std::vector<LinearConstraint>& candidateEquations() const
    {
      return _candidateEquations;
    }

  protected:

    void notify(AffineHullHandler::Event event)
    {
      for (std::size_t i = 0; i < _handlers.size(); ++i)
        _handlers[i]->notify(event, *this);
    }

    void resetEquations()
    {
      _equationSpace.reset(_n);
      for (std::size_t i = 0; i < _equations.size(); ++i)
        _equationSpace.addLazy(_equations[i].normal());
      _equationSpace.flushLazy();
    }

    void initializeEquations(const std::vector<LinearConstraint>& givenEquations)
    {
      assert(_equations.empty());

      _equationSpace.reset(_n);
      for (std::size_t i = 0; i < givenEquations.size(); ++i)
        _equationSpace.addLazy(givenEquations[i].normal());
      _equationSpace.flushLazy();

      for (std::size_t i = 0; i < givenEquations.size(); ++i)
      {
        if (!_equationSpace.isDependent(i))
          _equations.push_back(givenEquations[i]);
      }

      resetEquations();
    }

    std::size_t updateApproximateDirections()
    {
      std::size_t count = 0;
      for (std::size_t c = 0; c <= _n; ++c)
      {
        DSVectorReal vector;
        if (!_factorization.isBasic(c) && !_columns[c].approximateValid)
        {
          vector.clear();
          _factorization.computeApproximateKernelVector(c, _columns[c].approximateDirection, paramApproximateDirectionEpsilon);
          _columns[c].approximateValid = true;
          ++count;
        }
      }
      return count;
    }

    void computeAllowedDirectionColumn()
    {
      _allowedDirectionColumns.clear();
      for (std::size_t c = 0; c <= _n; ++c)
      {
        if (_factorization.isBasic(c))
          continue;

        if (_allowedDirectionColumns.empty())
        {
          _allowedDirectionColumns.push_back(c);
          continue;
        }
        if (_columns[c].requiredBefore(_columns[_allowedDirectionColumns.front()]))
        {
          _allowedDirectionColumns.clear();
          _allowedDirectionColumns.push_back(c);
        }
        else if (!_columns[_allowedDirectionColumns.front()].requiredBefore(_columns[c]))
        {
          _allowedDirectionColumns.push_back(c);
        }
      }

      for (std::size_t i = 0; i < _allowedDirectionColumns.size(); ++i)
      {
        assert(!_columns[_allowedDirectionColumns[i]].definesEquation);
      }
    }

    std::size_t selectDirectionColumnSparseApproximate()
    {
      std::size_t column = _allowedDirectionColumns[0];
      for (std::size_t i = 1; i < _allowedDirectionColumns.size(); ++i)
      {
        if (_columns[_allowedDirectionColumns[i]].sparsity < _columns[column].sparsity)
          column = _allowedDirectionColumns[i];
      }
      return column;
    }

    std::size_t selectDirectionColumnRandom()
    {
      return _allowedDirectionColumns[rand() % _allowedDirectionColumns.size()];
    }

    void updateExactDirection(std::size_t column)
    {
      Timer timer;
      timer.start();

      NonbasicColumn& col = _columns[column];
      _factorization.computeKernelVector(column, col.exactDirection);

      if (timer.time() > limitExactDirectionTime)
      {
        throw std::runtime_error("Affine hull computation reached time limit per exact direction computation!");
      }

      col.exactRhs = col.exactDirection[_n];
      col.exactDirection.reDim(_n);
      col.exactValid = true;

      // Also update approximate direction.

      col.approximateDirection = col.exactDirection;
      col.approximateValid = true;
    }

    bool checkDirectionDepends(std::size_t column)
    {
      if (_equationSpace.isDependent(_columns[column].exactDirection))
      {
        _columns[column].definesEquation = true;
        return true;
      }
      return false;
    }

    void addPoint(Vector& point, std::size_t column, bool invalidate)
    {
      _innerDescription.points.push_back(point);

      notify(AffineHullHandler::POINT_BEGIN);

      _basicColumns.push_back(column);
      _factorization.addRow(point, Rational(-1), column);

      if (invalidate)
      {
        DVectorRational lastRow;
        _factorization.getLastRow(lastRow);
        lastRow.reDim(_n);
        for (std::size_t c = 0; c <= _n; ++c)
        {
          if (_factorization.isBasic(c))
            continue;

          _columns[c].approximateValid = false;
          if (_columns[c].exactValid)
          {
            Rational product = _columns[c].exactDirection * lastRow;
            product -= _columns[c].exactRhs;
            if (product != 0)
            {
              _columns[c].exactValid = false;
              _columns[c].definesEquation = false;
            }
          }
        }
      }

      notify(AffineHullHandler::POINT_END);
    }

    void addRay(Vector& ray, std::size_t column, bool invalidate)
    {
      _innerDescription.rays.push_back(ray);

      notify(AffineHullHandler::RAY_BEGIN);

      _basicColumns.push_back(column);
      _factorization.addRow(ray, Rational(0), column);

      if (invalidate)
      {
        DVectorRational lastRow;
        _factorization.getLastRow(lastRow);
        lastRow.reDim(_n);
        for (std::size_t c = 0; c <= _n; ++c)
        {
          if (_factorization.isBasic(c))
            continue;

          _columns[c].approximateValid = false;
          if (_columns[c].exactValid)
          {
            Rational product = _columns[c].exactDirection * lastRow;
            if (product != 0)
            {
              _columns[c].exactValid = false;
              _columns[c].definesEquation = false;
            }
          }
        }
      }

      notify(AffineHullHandler::RAY_END);
    }

    void oracleQuery(AffineHullHandler::Event event, const soplex::DVectorRational& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(), HeuristicLevel minHeuristicLevel = 0,
      HeuristicLevel maxHeuristicLevel = std::numeric_limits<HeuristicLevel>::max())
    {
//       std::cerr << "Querying oracle for objective vector " << objective << " with heuristics in [" << _oracleMinHeuristicLevel << "," 
//         << _oracleMaxHeuristicLevel << "]" << std::endl;
      
      _oracleMaxHeuristicLevel = maxHeuristicLevel;
      _oracleMinHeuristicLevel = minHeuristicLevel;

      notify(event);

      _oracle->maximize(_result, objective, objectiveBound, minHeuristicLevel, maxHeuristicLevel);
      _oracleResultHeuristicLevel = _result.heuristicLevel();

      _result.sortPoints();
//       std::cerr << "Maximum value = " << _result.objectiveValue() << std::endl;

      if (event == AffineHullHandler::ORACLE_VERIFY_BEGIN)
        _verifySuccess = _result.isFeasible() && _result.points.front().objectiveValue == objectiveBound.value;

      notify(AffineHullHandler::Event(int(event) + 1));
    }

    void findLastPoint()
    {
      _directionPositiveVector.clear();

      oracleQuery(AffineHullHandler::ORACLE_ZERO_BEGIN, _directionPositiveVector);

      if (_result.isUnbounded())
      {
        throw std::runtime_error("Oracle claims unbounded for zero objective vector.");
      }
      else if (_result.isFeasible())
        addPoint(_result.points.front().vector, _n, false);
      else if (!_innerDescription.rays.empty())
      {
        throw std::runtime_error("Oracle claims infeasible after returning points or rays.");
      }
    }

    void mainLoop(HeuristicLevel lastCheapHeuristic, HeuristicLevel lastHeuristic)
    {
      while (lowerBound() < candidateUpperBound())
      {
        notify(AffineHullHandler::LOOP);

        // Zero objective case.

        if (_innerDescription.points.empty() && _innerDescription.rays.size() == candidateUpperBound())
        {
          findLastPoint();
          return;
        }

        if (paramApproximateDirections)
        {
          // Approximate directions.

          notify(AffineHullHandler::DIRECTIONS_APPROXIMATE_BEGIN);
          _lastApproximateDirectionSolves = updateApproximateDirections();
          notify(AffineHullHandler::DIRECTIONS_APPROXIMATE_END);
        }

        // Exact directions.

        notify(AffineHullHandler::DIRECTIONS_EXACT_BEGIN);
        _lastExactDirectionSolves = 0;
        _allowedDirectionColumns.clear();
        do
        {
          computeAllowedDirectionColumn();
          assert(!_allowedDirectionColumns.empty());

          if (paramApproximateDirections)
            _directionColumn = selectDirectionColumnSparseApproximate();
          else
            _directionColumn = selectDirectionColumnRandom();

          assert(_directionColumn <= _n);

          updateExactDirection(_directionColumn);
          ++_lastExactDirectionSolves;
          _directionPositiveVector = _columns[_directionColumn].exactDirection; // TODO: Avoid copying: a pointer should suffice!
          _lastExactDirectionNonzeros = std::numeric_limits<std::size_t>::max();
          _lastExactDirectionBitsize = std::numeric_limits<std::size_t>::max();
          notify(AffineHullHandler::DIRECTIONS_EXACT_COMPUTED);
        }
        while (checkDirectionDepends(_directionColumn));
        _directionNegativeVector = -_directionPositiveVector;

        if (_innerDescription.points.empty())
        {
          _objectiveValuePositive = soplex::infinity;
          _objectiveValueNegative = soplex::infinity;
        }
        else
        {
          _objectiveValuePositive = *_innerDescription.points.begin() * _directionPositiveVector;
          _objectiveValueNegative = -_objectiveValuePositive;
        }

        notify(AffineHullHandler::DIRECTIONS_EXACT_END);

        // Loop for maximization / minimization.

        HeuristicLevel lastMaximizeHeuristicLevel = _oracle->heuristicLevel() + 1;
        HeuristicLevel lastMinimizeHeuristicLevel = _oracle->heuristicLevel() + 1;
        Rational* pObjectiveValue =  NULL;
        bool addEquation = true;
        while (lastMaximizeHeuristicLevel > lastHeuristic || lastMinimizeHeuristicLevel > lastHeuristic)
        {
          // Maximize or minimize? Choose depending on last heuristic level (higher has preference).

          if (lastMaximizeHeuristicLevel >= lastMinimizeHeuristicLevel)
          {
            pObjectiveValue = &_objectiveValuePositive;
            HeuristicLevel min = lastMaximizeHeuristicLevel - 1;
            if (min < lastCheapHeuristic)
              min = lastHeuristic;
            HeuristicLevel max = lastMaximizeHeuristicLevel - 1;
            oracleQuery(AffineHullHandler::ORACLE_MAXIMIZE_BEGIN, _directionPositiveVector,
              ObjectiveBound(_objectiveValuePositive, true), min, max);
            lastMaximizeHeuristicLevel = _result.heuristicLevel();
          }
          else
          {
            pObjectiveValue = &_objectiveValueNegative;
            HeuristicLevel min = lastMinimizeHeuristicLevel - 1;
            if (min < lastCheapHeuristic)
              min = lastHeuristic;
            HeuristicLevel max = lastMinimizeHeuristicLevel - 1;
            oracleQuery(AffineHullHandler::ORACLE_MINIMIZE_BEGIN, _directionNegativeVector,
              ObjectiveBound(_objectiveValueNegative, true), min, max);
            lastMinimizeHeuristicLevel = _result.heuristicLevel();
          }

          // Evaluate result.

          if (_result.isInfeasible())
          {
            assert(_innerDescription.points.empty() && _innerDescription.rays.empty());
            _infeasible = true;
            _equations.clear();
            _equations.push_back(LinearConstraint('=', zeroVector(), Rational(1)));
            _candidateEquations.clear();
            notify(AffineHullHandler::EQUATION_FINAL);
            return;
          }
          else if (_result.isUnbounded())
          {
            addRay(_result.rays.front().vector, _directionColumn, true);
            addEquation = false;
            break;
          }
          else
          {
            assert(_result.isFeasible());

            for (std::size_t i = 0; i < _result.points.size(); ++i)
            {
              Vector& vector = _result.points[i].vector;
              const Rational& value = _result.points[i].objectiveValue;
              if (value == *pObjectiveValue)
                continue;

              if (_innerDescription.points.empty())
              {
                addPoint(vector, _n, true);
                _objectiveValuePositive = value;
                _objectiveValueNegative = -_objectiveValuePositive;
              }
              else
              {
                addPoint(vector, _directionColumn, true);
                addEquation = false;
                break;
              }
            }
            if (!addEquation)
              break;
          }
        }

        // If not enough was found, add candidate equation.

        if (addEquation)
        {
          LinearConstraint equation = LinearConstraint('=', denseToVector(_directionPositiveVector), _objectiveValuePositive);
          _equationSpace.add(equation.normal());
          _columns[_directionColumn].definesEquation = true;

          if (lastMinimizeHeuristicLevel > 0 || lastMaximizeHeuristicLevel > 0)
          {
            _candidateEquations.push_back(equation);
            notify(AffineHullHandler::EQUATION_CANDIDATE);
          }
          else
          {
            _equations.push_back(equation);
            notify(AffineHullHandler::EQUATION_FINAL);
          }
        }
      }
    }

  public:

    void run(const std::vector<LinearConstraint>& givenEquations)
    {
      _oracleMaxHeuristicLevel = _oracle->heuristicLevel();
      notify(AffineHullHandler::BEGIN);
      initializeEquations(givenEquations);
      notify(AffineHullHandler::EQUATIONS_INITIALIZED);

      mainLoop(paramLastCheapHeuristic, paramLastModerateHeuristic);

      if (lowerBound() < upperBound() && !_candidateEquations.empty())
      {
        // Verify k equations with only k+1 oracle calls. If one fails, continue with exact loop.

        notify(AffineHullHandler::VERIFY_BEGIN);

        _directionNegativeVector.clear();
        _objectiveValueNegative = 0;
        for (std::size_t i = 0; _verifySuccess && i < _candidateEquations.size(); ++i)
        {
          vectorToDense(_candidateEquations[i].normal(), _directionPositiveVector);
          _objectiveValuePositive = _candidateEquations[i].rhs();
          _directionNegativeVector -= _directionPositiveVector;
          _objectiveValueNegative -= _objectiveValuePositive;

          _verifyIndex = i;
          oracleQuery(AffineHullHandler::ORACLE_VERIFY_BEGIN, _directionPositiveVector,
            ObjectiveBound(_objectiveValuePositive, true));
        }

        // Last oracle call.

        if (_verifySuccess)
        {
          _verifyIndex = _candidateEquations.size();
          oracleQuery(AffineHullHandler::ORACLE_VERIFY_BEGIN, _directionNegativeVector,
            ObjectiveBound(_objectiveValueNegative, true));
        }
        notify(AffineHullHandler::VERIFY_END);

        // If successful, make equations regular ones.


        if (_verifySuccess)
        {
          std::copy(_candidateEquations.begin(), _candidateEquations.end(), std::back_inserter(_equations));
          _candidateEquations.clear();
        }
      }

      // Verification failed, so start main loop again with exact oracle.

      if (lowerBound() < upperBound())
      {
        for (std::size_t c = 0; c <= _n; ++c)
          _columns[c].definesEquation = false;

        resetEquations();
        _candidateEquations.clear();

        mainLoop(std::numeric_limits<HeuristicLevel>::max(), 0);
      }

      notify(AffineHullHandler::END);
    }

  public:
    HeuristicLevel paramLastCheapHeuristic;
    HeuristicLevel paramLastModerateHeuristic;
    bool paramApproximateDirections;
    double paramApproximateDirectionEpsilon;
    double limitExactDirectionTime;

  protected:
    std::vector<AffineHullHandler*>& _handlers;
    std::shared_ptr<OracleBase> _oracle;
    InnerDescription& _innerDescription;
    std::vector<LinearConstraint>& _equations;
    std::vector<LinearConstraint> _candidateEquations;
    bool _infeasible;
    std::size_t _n;
    std::vector<std::size_t> _basicColumns;

    OracleResult _result;
    VectorSpaceGenerators _equationSpace;
    std::size_t _directionColumn;
    std::vector<std::size_t> _allowedDirectionColumns;
    soplex::DVectorRational _directionPositiveVector;
    soplex::DVectorRational _directionNegativeVector;
    Factorization _factorization;
    std::vector<NonbasicColumn> _columns;
    Rational _objectiveValuePositive;
    Rational _objectiveValueNegative;

    std::size_t _verifyIndex;
    bool _verifySuccess;
    std::size_t _lastApproximateDirectionSolves;
    std::size_t _lastExactDirectionSolves;
    std::size_t _lastExactDirectionNonzeros;
    std::size_t _lastExactDirectionBitsize;
    HeuristicLevel _oracleMaxHeuristicLevel;
    HeuristicLevel _oracleMinHeuristicLevel;
    HeuristicLevel _oracleResultHeuristicLevel;
  };

  AffineHullHandler::AffineHullHandler()
  {

  }

  AffineHullHandler::~AffineHullHandler()
  {

  }

  void affineHull(const std::shared_ptr<OracleBase>& oracle, InnerDescription& resultInnerDescription,
    std::vector<LinearConstraint>& resultOuterDescription, std::vector<AffineHullHandler*>& handlers,
    HeuristicLevel lastCheapHeuristic, HeuristicLevel lastModerateHeuristic,
    const std::vector<LinearConstraint>& givenEquations, bool approximateDirections, double exactDirectionTimeLimit)
  {
    XAffineHull algorithm(handlers, oracle, resultInnerDescription, resultOuterDescription);
    algorithm.paramLastCheapHeuristic = lastCheapHeuristic;
    algorithm.paramLastModerateHeuristic = lastModerateHeuristic;
    algorithm.paramApproximateDirections = approximateDirections;
    algorithm.paramApproximateDirectionEpsilon = 1.0e-7;
    algorithm.limitExactDirectionTime = exactDirectionTimeLimit;
    algorithm.run(givenEquations);
  }

  DebugAffineHullHandler::DebugAffineHullHandler(std::ostream& stream, bool printPointsAndRays, bool printEquations,
    bool printDirections)
    : _stream(stream), _printPointsAndRays(printPointsAndRays), _printEquations(printEquations), _printDirections(printDirections)
  {

  }

  DebugAffineHullHandler::~DebugAffineHullHandler()
  {

  }

  void DebugAffineHullHandler::notify(AffineHullHandler::Event event, AffineHullState& state)
  {
    switch (event)
    {
      case BEGIN:
        _stream << "AH: Computing affine hull in ambient space of dimension " << state.space().dimension() << ".\n";
      break;
      case LOOP:
        _stream << "AH: " << state.lowerBound() << " <= dimension ";
        if (state.candidateUpperBound() < state.upperBound())
          _stream << "<?= " << state.candidateUpperBound();
        _stream << "<= " << state.upperBound() << ".\n";
      break;
      case EQUATIONS_INITIALIZED:
        _stream << "AH: Initialized " << state.numEquations() << " given independent equations.\n";
      break;
      case DIRECTIONS_APPROXIMATE_BEGIN:
        _stream << "AH: Computing approximate directions.\n";
      break;
      case DIRECTIONS_APPROXIMATE_END:
        _stream << "AH: Computed " << state.directionApproximateSolves() << " approximate directions.\n";
      break;
      case DIRECTIONS_EXACT_BEGIN:
        _stream << "AH: Computing exact directions.\n";
      break;
      case DIRECTIONS_EXACT_COMPUTED:
        _stream << "AH: Computed an exact direction.\n";
      break;
      case DIRECTIONS_EXACT_END:
        _stream << "AH: Computed " << state.directionExactSolves() << " exact directions.\n";
        _stream << "AH: Current direction vector has " << state.directionNonzeros() << " nonzeros and an encoding length of "
          << state.directionBitsize() << " bits.\n";
        if (_printDirections)
        {
          _stream << "AH: ";
          state.space().printLinearForm(_stream, denseToVector(state.directionVector(), false));
          _stream << "\n";
        }
      break;
      case ORACLE_ZERO_BEGIN:
        _stream << "AH: Maximizing zero objective";
        if (state.oracleMaxHeuristicLevel() == std::numeric_limits<HeuristicLevel>::max())
          _stream << " (heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
        else
          _stream << " (" << state.oracleMaxHeuristicLevel() << " >= heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
      break;
      case ORACLE_MAXIMIZE_BEGIN:
        _stream << "AH: Maximizing direction objective";
        if (state.oracleMaxHeuristicLevel() == std::numeric_limits<HeuristicLevel>::max())
          _stream << " (heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
        else
          _stream << " (" << state.oracleMaxHeuristicLevel() << " >= heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
      break;
      case ORACLE_MINIMIZE_BEGIN:
        _stream << "AH: Minimizing direction objective";
        if (state.oracleMaxHeuristicLevel() == std::numeric_limits<HeuristicLevel>::max())
          _stream << " (heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
        else
          _stream << " (" << state.oracleMaxHeuristicLevel() << " >= heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
      break;
      case ORACLE_VERIFY_BEGIN:
        _stream << "AH: Maximizing verification objective";
        if (state.oracleMaxHeuristicLevel() == std::numeric_limits<HeuristicLevel>::max())
          _stream << " (heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
        else
          _stream << " (" << state.oracleMaxHeuristicLevel() << " >= heurLevel >= " << state.oracleMinHeuristicLevel() << ")\n";
      break;
      case ORACLE_ZERO_END:
      case ORACLE_MAXIMIZE_END:
      case ORACLE_MINIMIZE_END:
      case ORACLE_VERIFY_END:
        if (state.oracleNumPoints() > 0)
          _stream << "AH: Oracle returned " << state.oracleNumPoints() << " points";
        else if (state.oracleNumRays() > 0)
          _stream << "AH: Oracle returned " << state.oracleNumRays() << " rays";
        else
          _stream << "AH: Oracle claimed infeasible";
        _stream << ", heurLevel = " << state.oracleResultHeuristicLevel() << "\n";
      break;
      case POINT_BEGIN:
        _stream << "AH: Adding a point.\n";
        if (_printPointsAndRays)
        {
          _stream << "AH: ";
          state.space().printVector(_stream, state.innerDescription().points.back());
          _stream << "\n";
        }
      break;
      case POINT_END:
        _stream << "AH: Added a point.\n";
      break;
      case RAY_BEGIN:
        _stream << "AH: Adding a ray.\n";
        if (_printPointsAndRays)
        {
          _stream << "AH: ";
          state.space().printVector(_stream, state.innerDescription().rays.back());
          _stream << "\n";
        }
      break;
      case RAY_END:
        _stream << "AH: Added a ray.\n";
      break;
      case EQUATION_CANDIDATE:
        _stream << "AH: Added a candidate equation.\n";
        if (_printEquations)
        {
          _stream << "AH: ";
          state.space().printLinearConstraint(_stream, state.candidateEquations().back());
          _stream << "\n";
        }
      break;
      case EQUATION_FINAL:
        _stream << "AH: Added an equation.\n";
        if (_printEquations)
        {
          _stream << "AH: ";
          state.space().printLinearConstraint(_stream, state.outerDescription().back());
          _stream << "\n";
        }
      break;
      case VERIFY_BEGIN:
        _stream << "AH: Starting to verify " << state.numCandidateEquations() << " candidate equations.\n";
      break;
      case VERIFY_END:
        if (state.verifySuccess())
          _stream << "AH: Verified " << state.numCandidateEquations() << " candidate equations.\n";
        else
          _stream << "AH: Failed to verify " << state.numCandidateEquations() << " candidate equations.\n";
      break;
      case END:
        _stream << "AH: Dimension is equal to " << state.lowerBound() << "." << std::endl;
      break;
      default:
        _stream << "AH: Unhandled event " << event << ".\n" << std::endl;
      break;
    }
    _stream << std::flush;
  }

  StatisticsAffineHullHandler::StatisticsAffineHullHandler()
  {
    reset();
  }

  StatisticsAffineHullHandler::~StatisticsAffineHullHandler()
  {

  }

  void StatisticsAffineHullHandler::reset()
  {
    _timer.reset();
    _numOracleQueries = 0;
    _numHeuristicLevelAnswers.clear();
    _numDirectionApproximateSolves = 0;
    _numDirectionExactSolves = 0;
    _sumDirectionNonzeros = 0;
    _maxDirectionBitsize = 0;
    _numFactorizations = 0;
    _timeApproximateDirections = 0.0;
    _timeExactDirections = 0.0;
    _timeFactorizations = 0.0;
    _timeOracles = 0.0;
    _timeMainLoop = 0.0;
    _timeVerification = 0.0;
    _timeLastExactDirectionBegin = 0.0;
    _timeLastMainLoop = 0.0;
    _timeLastVerification = 0.0;
    _timeLastEvent = 0.0;
#ifdef IPO_DEBUG
    _lastEvent = END;
#endif
  }

  void StatisticsAffineHullHandler::notify(AffineHullHandler::Event event, AffineHullState& state)
  {
    if (event == BEGIN)
    {
      _timer.reset();
      _timer.start();
      _numHeuristicLevelAnswers.resize(state.oracleMaxHeuristicLevel() + 1, 0);
    }
    double time = _timer.time();
    switch (event)
    {
      case BEGIN:
        _timeLastMainLoop = time;
      break;
      case DIRECTIONS_APPROXIMATE_END:
        _numDirectionApproximateSolves += state.directionApproximateSolves();
#ifdef IPO_DEBUG
        assert(_lastEvent == DIRECTIONS_APPROXIMATE_BEGIN);
#endif
        _timeApproximateDirections += time - _timeLastEvent;
      break;
      case DIRECTIONS_EXACT_BEGIN:
        _timeLastExactDirectionBegin = time;
      break;
      case DIRECTIONS_EXACT_COMPUTED:
        _numDirectionExactSolves++;
        _sumDirectionNonzeros += state.directionNonzeros();
        _maxDirectionBitsize = std::max(_maxDirectionBitsize, state.directionBitsize());
      break;
      case DIRECTIONS_EXACT_END:
        _timeExactDirections += time - _timeLastExactDirectionBegin;
      break;
      case ORACLE_ZERO_END:
      case ORACLE_MAXIMIZE_END:
      case ORACLE_MINIMIZE_END:
      case ORACLE_VERIFY_END:
#ifdef IPO_DEBUG
        assert(_lastEvent == ORACLE_ZERO_BEGIN || _lastEvent == ORACLE_MAXIMIZE_BEGIN || _lastEvent == ORACLE_MINIMIZE_BEGIN
          || _lastEvent == ORACLE_VERIFY_BEGIN);
#endif
        _numOracleQueries++;
        _timeOracles += time - _timeLastEvent;
        assert(state.oracleResultHeuristicLevel() < _numHeuristicLevelAnswers.size());
        _numHeuristicLevelAnswers[state.oracleResultHeuristicLevel()]++;
      break;
      case POINT_END:
      case RAY_END:
        _numFactorizations++;
#ifdef IPO_DEBUG
        assert(_lastEvent == POINT_BEGIN || _lastEvent == RAY_BEGIN);
#endif
        _timeFactorizations += time - _timeLastEvent;
      break;
      case VERIFY_BEGIN:
        _timeMainLoop += time - _timeLastMainLoop;
        _timeLastVerification = time;
      break;
      case VERIFY_END:
        _timeVerification += time - _timeLastVerification;
        _timeLastMainLoop = time;
      break;
      case END:
        _timeMainLoop += time - _timeLastMainLoop;
      break;
    }
#ifdef IPO_DEBUG
    _lastEvent = event;
#endif
    _timeLastEvent = time;
  }




  namespace AffineHull {

    class Factorization
    {
    public:
      Factorization(std::size_t numColumns) :
          _numRows(0), _columnVectors(numColumns), _columnBasics(numColumns, std::numeric_limits<std::size_t>::max())
      {

      }

      ~Factorization()
      {

      }

      inline std::size_t numRows() const
      {
        return _numRows;
      }

      inline std::size_t numColumns() const
      {
        return _columnVectors.size();
      }

      inline bool isBasic(std::size_t columnIndex) const
      {
        return _columnBasics[columnIndex] != std::numeric_limits<std::size_t>::max();
      }

      void solveRightFactor(DVectorRational& vector) const
      {
        for (std::size_t b = _numRows; b > 0; b--)
        {
          std::size_t pivotRow = b - 1;
          const DSVectorRational& column = _factorUcolumnVectors[pivotRow];
          assert(column.index(column.size() - 1) == pivotRow);
          vector[pivotRow] /= column.value(column.size() - 1);

          for (int p = column.size() - 2; p >= 0; --p)
          {
            vector[column.index(p)] -= column.value(p) * vector[pivotRow];
          }
        }
      }

      void solveRightApproximateFactor(DVectorReal& vector) const
      {
        for (std::size_t b = _numRows; b > 0; b--)
        {
          std::size_t pivotRow = b - 1;
          const DSVectorReal& column = _approximateFactorUcolumnVectors[pivotRow];
          assert(column.index(column.size() - 1) == pivotRow);
          vector[pivotRow] /= column.value(column.size() - 1);

          for (int p = column.size() - 2; p >= 0; --p)
          {
            vector[column.index(p)] -= column.value(p) * vector[pivotRow];
          }
        }
      }

      void solveLeftFactor(DVectorRational& vector) const
      {
        for (std::size_t pivotRow = 0; pivotRow < numRows(); ++pivotRow)
        {
          const DSVectorRational& column = _factorLcolumnVectors[pivotRow];
          assert(column.index(0) == pivotRow);
          vector[pivotRow] /= column.value(0);

          for (int p = column.size() - 1; p >= 1; --p)
          {
            vector[column.index(p)] -= column.value(p) * vector[pivotRow];
          }
        }
      }

      void solveLeftApproximateFactor(DVectorReal& vector) const
      {
        for (std::size_t pivotRow = 0; pivotRow < numRows(); ++pivotRow)
        {
          const DSVectorReal& column = _approximateFactorLcolumnVectors[pivotRow];
          assert(column.index(0) == pivotRow);
          vector[pivotRow] /= column.value(0);

          for (int p = column.size() - 1; p >= 1; --p)
          {
            vector[column.index(p)] -= column.value(p) * vector[pivotRow];
          }
        }
      }

      void computeKernelVector(std::size_t columnIndex, soplex::DVectorRational& result) const
      {
        soplex::DVectorRational rhs;
        rhs.reDim(numRows(), true);
        const SVectorRational& col = _columnVectors[columnIndex];
        for (int p = col.size() - 1; p >= 0; --p)
          rhs[col.index(p)] = -col.value(p);

        solveLeftFactor(rhs);
        solveRightFactor(rhs);

        std::size_t size = 1;
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        {
          if (rhs[b] != 0)
            ++size;
        }
        result.reDim(numColumns(), false);
        result.clear();
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
          result[_basicColumns[b]] = rhs[b];
        result[columnIndex] = Rational(1);
      }

      void computeApproximateKernelVector(std::size_t columnIndex, soplex::DVectorReal& result, double epsilon) const
      {
        DVectorReal rhs;
        rhs.reDim(numRows(), true);
        const SVectorRational& col = _columnVectors[columnIndex];
        for (int p = col.size() - 1; p >= 0; --p)
          rhs[col.index(p)] = -(double) col.value(p);

        solveLeftApproximateFactor(rhs);
        solveRightApproximateFactor(rhs);

        result.reDim(numColumns());
        result.clear();
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        {
          if (fabs(rhs[b]) > epsilon)
            result[_basicColumns[b]] = rhs[b];
        }
        result[columnIndex] = 1.0;
      }

      void dump()
      {
        std::cerr << "ACM: Matrix is " << numRows() << "x" << numColumns() << "\n";
        std::cerr << "ACM: Basis:";
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
          std::cerr << " " << _basicColumns[b];
        std::cerr << std::endl;
        for (std::size_t c = 0; c < numColumns(); ++c)
        {
          std::cerr << "ACM: Column " << c << ":";
          for (int p = 0; p < _columnVectors[c].size(); ++p)
          {
            std::cerr << "{" << _columnVectors[c].index(p) << ":" << _columnVectors[c].value(p) << "}";
          }
          std::cerr << std::endl;
        }
        std::cerr << "ACM: L:\n";
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        {
          std::cerr << "ACM: Column " << b << " -> " << _basicColumns[b] << ":";
          for (int p = 0; p < _factorLcolumnVectors[b].size(); ++p)
          {
            std::cerr << "{" << _factorLcolumnVectors[b].index(p) << ":" << _factorLcolumnVectors[b].value(p) << "}";
          }
          std::cerr << std::endl;
        }
        std::cerr << "ACM: U:\n";
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        {
          std::cerr << "ACM: Column " << b << " -> " << _basicColumns[b] << ":";
          for (int p = 0; p < _factorUcolumnVectors[b].size(); ++p)
          {
            std::cerr << "{" << _factorUcolumnVectors[b].index(p) << ":" << _factorUcolumnVectors[b].value(p) << "}";
          }
          std::cerr << std::endl;
        }

      }

      void addRow(const Vector& row, const Rational& last, std::size_t newBasicColumnIndex)
      {
        assert(_columnBasics[newBasicColumnIndex] == std::numeric_limits<std::size_t>::max());

        /// Compute new column which is aligned right to the right factor.

        DVectorRational rhs;
        rhs.reDim(numRows(), true);
        rhs.assign(_columnVectors[newBasicColumnIndex]);
        solveLeftFactor(rhs);
        _factorUcolumnVectors.push_back(DSVectorRational());
        _factorUcolumnVectors.back() = rhs;
        _factorUcolumnVectors.back().sort();
        _basicColumns.push_back(newBasicColumnIndex);
        _columnBasics[newBasicColumnIndex] = numRows();

        _approximateFactorUcolumnVectors.push_back(DSVectorReal());
        for (int p = 0; p < _factorUcolumnVectors.back().size(); ++p)
        {
          _approximateFactorUcolumnVectors.back().add(_factorUcolumnVectors.back().index(p),
              (double) _factorUcolumnVectors.back().value(p));
        }

        /// Copy new row to worker which is a row vector indexed the ordered basis.

        DVectorRational worker;
        worker.reDim(numRows() + 1, true);
        for (int p = row.size() - 1; p >= 0; --p)
        {
          if (_columnBasics[row.index(p)] != std::numeric_limits<std::size_t>::max())
            worker[_columnBasics[row.index(p)]] = row.value(p);
        }
        if (last != 0 && _columnBasics[numColumns() - 1] != std::numeric_limits<std::size_t>::max())
          worker[_columnBasics[numColumns() - 1]] = last;

        std::vector<int> positions(_basicColumns.size(), 0);
        for (std::size_t b = 0; b < _basicColumns.size() - 1; ++b)
        {
          if (worker[b] == 0)
            continue;

          std::size_t pivotRowIndex = b;
          std::size_t pivotColumnIndex = _basicColumns[b];
          assert(_factorUcolumnVectors[b].index(_factorUcolumnVectors[b].size() - 1) == pivotRowIndex);
          const Rational& pivotElement = _factorUcolumnVectors[b].value(_factorUcolumnVectors[b].size() - 1);
          assert(pivotElement != 0);
          const Rational lambda = -worker[b] / pivotElement;

          for (std::size_t i = b + 1; i < _basicColumns.size(); ++i)
          {
            const SVectorRational& currentRightColumn = _factorUcolumnVectors[i];
            while (positions[i] < currentRightColumn.size() && currentRightColumn.index(positions[i]) < pivotRowIndex)
              positions[i]++;

            /// Pivot row entry may be zero:
            if (positions[i] == currentRightColumn.size() || currentRightColumn.index(positions[i]) > pivotRowIndex)
              continue;

            assert(currentRightColumn.index(positions[i]) == pivotRowIndex);

            worker[i] += lambda * currentRightColumn.value(positions[i]);
          }

          /// Put -lambda into left factor.
          _factorLcolumnVectors[b].add(_numRows, -lambda);
          _approximateFactorLcolumnVectors[b].add(_numRows, -(double) lambda);
        }
        if (worker[numRows()] == 0)
        {
          throw std::runtime_error("Invalid call to AffineComputationMatrix::addRow: Basis matrix is singular!");
        }

        /// Diagonal entries.
        _factorLcolumnVectors.push_back(DSVectorRational());
        _factorLcolumnVectors[numRows()].add(numRows(), Rational(1));
        _factorUcolumnVectors[numRows()].add(numRows(), worker[numRows()]);

        _approximateFactorLcolumnVectors.push_back(DSVectorReal());
        _approximateFactorLcolumnVectors[numRows()].add(numRows(), 1.0);
        _approximateFactorUcolumnVectors[numRows()].add(numRows(), (double) worker[numRows()]);

        /// Physically add the row.
        for (int p = row.size() - 1; p >= 0; --p)
        {
          _columnVectors[row.index(p)].add(numRows(), row.value(p));
        }
        if (last != 0)
          _columnVectors[numColumns() - 1].add(numRows(), last);
        ++_numRows;
      }

      void getLastRow(DVectorRational& row) const
      {
        row.reDim(numColumns(), true);
        if (numRows() == 0)
          return;
        int rowIndex = numRows() - 1;
        for (std::size_t c = 0; c < numColumns(); ++c)
        {
          int p = _columnVectors[c].size() - 1;
          if (p >= 0 && _columnVectors[c].index(p) == rowIndex)
          {
            row[c] = _columnVectors[c].value(p);
          }
        }
      }

    protected:
      std::size_t _numRows;
      std::vector<DSVectorRational> _columnVectors;
      std::vector<std::size_t> _basicColumns;
      std::vector<std::size_t> _columnBasics;
      std::vector<DSVectorRational> _factorLcolumnVectors;
      std::vector<DSVectorRational> _factorUcolumnVectors;
      std::vector<DSVectorReal> _approximateFactorLcolumnVectors;
      std::vector<DSVectorReal> _approximateFactorUcolumnVectors;
    };

    /// Main implementation

    struct NonbasicColumn
    {
      bool exactValid;
      soplex::DVectorRational exactDirection;
      Rational exactRhs;
      bool approxValid;
      soplex::DVectorReal approxDirection;
      std::size_t sparsity;
      bool definesEquation;
      bool avoid;

      NonbasicColumn()
        : exactDirection(0)
      {
        exactValid = false;
        approxValid = false;
        definesEquation = false;
        avoid = false;
        sparsity = 0;
      }

      bool operator<(const NonbasicColumn& other) const
      {
        if (definesEquation)
          return false;
        if (other.definesEquation)
          return true;

        if (avoid)
          return false;
        if (other.avoid)
          return true;

        return sparsity < other.sparsity;
      }
    };

    class Implementation
    {
    public:
      Implementation()
        : _oracle(NULL), _output(NULL), _factorization(NULL)
      {
        _directionColumn = std::numeric_limits<std::size_t>::max();
        _commonValue = infinity;
        _numCacheQueries = 0;
        _numCacheHits = 0;
        _numHeuristicCalls = 0;
        _numOracleCalls = 0;
        _numApproximateDirectionSolves = 0;
        _numExactDirectionSolves = 0;
        _lastDirectionBitsize = 0;
        _maxDirectionBitsize = 0;
        _infeasible = false;
        _espace = new VectorSpaceGenerators();
      }

      virtual ~Implementation()
      {
        if (_factorization != NULL)
          delete _factorization;
        delete _espace;
      }

      int dimensionLowerBound() const
      {
        return int(_spanningPoints.size() + _spanningRays.size()) - 1;
      }

      int dimensionUnsafeUpperBound() const
      {
        if (_infeasible)
          return -1;
        else
          return int(n()) - int(_irredundantEquations.size() + _potentialEquations.size());
      }

      int dimensionSafeUpperBound() const
      {
        if (_infeasible)
          return -1;
        else
          return int(n()) - int(_irredundantEquations.size());
      }

      void initializeEquations()
      {
        assert(_potentialEquations.empty());
        assert(_redundantEquations.empty());

        _espace->reset(n());
        for (std::size_t i = 0; i < _irredundantEquations.size(); ++i)
          _espace->addLazy(_irredundantEquations[i].normal());
        _espace->flushLazy();

        std::vector<LinearConstraint> newIrredundantEquations;
        std::vector<LinearConstraint> newRedundantEquations;
        for (std::size_t i = 0; i < _irredundantEquations.size(); ++i)
        {
          if (_espace->isDependent(i))
            newIrredundantEquations.push_back(_irredundantEquations[i]);
          else
            newRedundantEquations.push_back(_irredundantEquations[i]);
        }
        _irredundantEquations = newIrredundantEquations;
        _redundantEquations = newRedundantEquations;
      }

//       void removePotentialEquations()
//       {
// #ifdef IPO_DEBUG
//         // Check that we first have the irredundant and then the potential equations.
//
//         int lastIrredundant = -1;
//         for (std::size_t i = 0; i < _irredundantEquations.size(); ++i)
//         {
//           if (_irredundantEquations[i] > lastIrredundant)
//             lastIrredundant = _irredundantEquations[i];
//         }
//         int firstPotential = _equations->num();
//         for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
//         {
//           if (_potentialEquations[i] < firstPotential)
//             firstPotential = _potentialEquations[i];
//         }
//         assert(firstPotential > lastIrredundant);
// #endif
//
//         /// Remove equations.
//
//         std::vector<int> remove;
//         for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
//           remove.push_back(_potentialEquations[i]);
//         _equations->remove(&remove[0], remove.size());
//
//         /// Recompute dependencies.
//         initializeEquations();
//       }

      void removeRedundantEquations()
      {
        std::size_t oldSize = _redundantEquations.size();
        _redundantEquations.clear();

        _output->onRemovedRedundantEquations(oldSize);
      }

      void updateApproximateDirections()
      {
        for (std::size_t c = 0; c <= n(); ++c)
        {
          DSVectorReal vector;
          if (!_factorization->isBasic(c) && !_columns[c].approxValid)
          {
            ++_numApproximateDirectionSolves;
            vector.clear();
            _factorization->computeApproximateKernelVector(c, _columns[c].approxDirection, 1.0e-7);
            _columns[c].approxValid = true;
          }
        }
      }

      void oracleMaximize(std::size_t minHeuristic)
      {
        _oracle->maximize(_result, _directionVector, _objectiveBound,
          std::numeric_limits<std::size_t>::max(), minHeuristic);

        if (_result.isInfeasible())
        {
          _infeasible = true;
          std::copy(_irredundantEquations.begin(), _irredundantEquations.end(), std::back_inserter(_redundantEquations));
          _irredundantEquations.resize(1);
          _irredundantEquations[0] = LinearConstraint('=', zeroVector(), Rational(1));
          _potentialEquations.clear();
          _output->onEquation();
        }
      }

      void addPoint(Vector& point, std::size_t column, bool invalidate)
      {
        _spanningPoints.push_back(point);
        _basicColumns.push_back(column);
        _factorization->addRow(point, Rational(-1), column);

        if (invalidate)
        {
          DVectorRational lastRow;
          _factorization->getLastRow(lastRow);
          for (std::size_t c = 0; c <= n(); ++c)
          {
            if (_factorization->isBasic(c))
              continue;

            _columns[c].approxValid = false;
            if (_columns[c].exactValid)
            {
              Rational product = _columns[c].exactDirection * lastRow;
              product -= _columns[c].exactRhs;
              if (product != 0)
              {
                _columns[c].exactValid = false;
                _columns[c].definesEquation = false;
              }
            }
          }
        }
      }

      void addRay(Vector& ray, std::size_t column, bool invalidate)
      {
        _spanningRays.push_back(ray);
        _basicColumns.push_back(column);
        _factorization->addRow(ray, Rational(0), column);

        if (invalidate)
        {
          DVectorRational lastRow;
          _factorization->getLastRow(lastRow);
          for (std::size_t c = 0; c <= n(); ++c)
          {
            if (_factorization->isBasic(c))
              continue;

            _columns[c].approxValid = false;
            if (_columns[c].exactValid)
            {
              Rational product = _columns[c].exactDirection * lastRow;
              if (product != 0)
              {
                _columns[c].exactValid = false;
                _columns[c].definesEquation = false;
              }
            }
          }
        }
      }

      void findLastPoint()
      {
        _directionVector.clear();
        _directionVector.reDim(n());

        _output->onBeforeOracleZero();

        _oracle->maximize(_result, _directionVector);

        _output->onAfterOracleZero(_result.points.size());

        if (_result.isUnbounded())
        {
          throw std::runtime_error("Oracle claims unbounded for zero objective vector.");
        }
        else if (_result.isInfeasible())
        {
          if (_spanningRays.empty())
            return;
          else
            throw std::runtime_error("Oracle claims infeasible after returning points or rays.");
        }
        else
        {
          assert(_result.isFeasible());

          _output->onBeforePoint(false);
          addPoint(_result.points.front().vector, n(), false);
          _output->onAfterPoint(false);

          return;
        }
      }

      void selectDirectionColumn()
      {
        _directionColumn = std::numeric_limits<std::size_t>::max();

#ifdef RANDOM_COBASIS_INDEX
        std::vector<std::size_t> candidates;
        bool avoid = true;

        for (std::size_t c = 0; c <= n(); ++c)
        {
          if (_factorization->isBasic(c))
            continue;
          if (_columns[c].definesEquation)
            continue;

          candidates.push_back(c);
          if (!_columns[c].avoid)
            avoid = false;
        }

        if (!avoid)
        {
          for (std::size_t i = 0; i < candidates.size(); ++i)
          {
            if (_columns[candidates[i]].avoid)
            {
              candidates.erase(candidates.begin() + i);
              i--;
            }
          }
        }

        std::random_shuffle(candidates.begin(), candidates.end());
        _directionColumn = candidates[0];
#else
        for (std::size_t c = 0; c <= n(); ++c)
        {
          if (_factorization->isBasic(c))
            continue;

          if (_directionColumn == std::numeric_limits<std::size_t>::max()
            || _columns[c] < _columns[_directionColumn])
          {
            _directionColumn = c;
          }
        }
#endif

        if (_columns[_directionColumn].definesEquation)
        {
          throw std::runtime_error(
            "AffineHull: Selected column defines an equation - there seems to be a bug!");
        }
      }

      void updateExactDirection(std::size_t column)
      {
        ++_numExactDirectionSolves;
        NonbasicColumn& col = _columns[column];
        _factorization->computeKernelVector(column, col.exactDirection);

        col.exactRhs = col.exactDirection[n()];
        col.exactValid = true;

        /// Also update approximate direction.

        col.approxDirection = col.exactDirection;
        col.approxValid = true;
      }

      bool checkDirectionDepends()
      {
        if (_espace->isDependent(_columns[_directionColumn].exactDirection))
        {
          _columns[_directionColumn].definesEquation = true;
          return true;
        }
        return false;
      }

      int getUpperBound(int minHeuristic)
      {
        int bound = int(n()) - int(_irredundantEquations.size());
        if (minHeuristic > 0)
          bound -= int(_potentialEquations.size());
        return bound;
      }

      void mainLoop(std::size_t minHeuristic)
      {
        while (dimensionLowerBound() < (minHeuristic > 0 ? dimensionUnsafeUpperBound() :
          dimensionSafeUpperBound()))
        {
          /// Approximate directions.

          _output->onBeforeApproximateDirections();
          std::size_t oldNumApproximateDirectionSolves = _numApproximateDirectionSolves;
          updateApproximateDirections();
          _output->onAfterApproximateDirections(
            _numApproximateDirectionSolves - oldNumApproximateDirectionSolves);

          /// Zero objective case.

          if (_spanningPoints.empty()
            && _spanningRays.size() == getUpperBound(minHeuristic))
          {
            return findLastPoint();
          }

          /// Exact directions.

          std::size_t oldNumExactDirectionSolves = _numExactDirectionSolves;
          _output->onBeforeExactDirections();
          do
          {
            selectDirectionColumn();
            assert(_directionColumn <= n());
            updateExactDirection(_directionColumn);
          }
          while (checkDirectionDepends());

          // We obtain the direction vector from a nonbasic column.
          _lastDirectionBitsize = 0;
          for (std::size_t v = 0; v < n(); ++v)
          {
            _directionVector[v] = _columns[_directionColumn].exactDirection[v];
            _approximateDirectionVector[v] = double(_directionVector[v]);
            if (_directionVector[v] != 0)
              _lastDirectionBitsize += _directionVector[v].sizeInBase(2);
          }

          if (_spanningPoints.empty())
            _commonValue = soplex::infinity;
          else
            _commonValue = *_spanningPoints.begin() * _directionVector;

          // Measure bitsize.

          _maxDirectionBitsize = std::max(_maxDirectionBitsize, _lastDirectionBitsize);
          _output->onAfterExactDirections(_numExactDirectionSolves - oldNumExactDirectionSolves);

          /*
           * TODO: To avoid oracle calls of low heuristic levels we might maximize/minimize
           * specific levels taking turns.
           */

          /// Maximize direction.

          _objectiveBound.value = _commonValue;
          _objectiveBound.strict = true;
          _output->onBeforeOracleMaximize();
          oracleMaximize(minHeuristic);
          _output->onAfterOracleMaximize(_result.points.size(), _result.rays.size());

          if (_result.isInfeasible())
          {
            if (!_spanningPoints.empty() || !_spanningRays.empty())
            {
              throw std::runtime_error(
                "AffineHull: Oracle claims infeasible though we know some points/rays.");
            }
            break;
          }
          else if (_result.isUnbounded())
          {
            _output->onBeforeRay();
            addRay(_result.rays.front().vector, _directionColumn, true);
            _output->onAfterRay();
            continue;
          }
          else
          {
            assert(_result.isFeasible());

            if (_spanningPoints.empty())
            {
              // Case |S| = 0: We might add two points.

              const Rational& firstObjectiveValue = _result.points.front().objectiveValue;
              Vector* firstPoint = &_result.points.front().vector;
              Vector* secondPoint = NULL;
              for (std::size_t i = 1; i < _result.points.size(); ++i)
              {
                if (_result.points[i].objectiveValue != firstObjectiveValue)
                {
                  secondPoint = &_result.points[i].vector;
                  break;
                }
              }

              // We always add the maximizer.

              _output->onBeforePoint(secondPoint > 0);
              addPoint(*firstPoint, n(), true);
              if (secondPoint != NULL)
              {
                // If there exists a point with different objective value, we add it as well.

                 std::size_t differingVariable = differingIndex(*firstPoint, *secondPoint);
                 assert(differingVariable != std::numeric_limits<std::size_t>::max());
                 addPoint(*secondPoint, differingVariable, true);
              }
              _output->onAfterPoint(secondPoint != NULL);
              if (secondPoint != NULL)
                continue;
              else
                _commonValue = _result.points.front().objectiveValue;
            }
            else
            {
              // We add any point that has objective different from common value.

              Vector* differentObjectiveVector = NULL;
              for (std::size_t i = 0; i < _result.points.size(); ++i)
              {
                if (_result.points[i].objectiveValue != _commonValue)
                {
                  differentObjectiveVector = &_result.points[i].vector;
                  break;
                }
              }
              if (differentObjectiveVector != NULL)
              {
                _output->onBeforePoint(false);
                addPoint(*differentObjectiveVector, _directionColumn, true);
                _output->onAfterPoint(false);
                continue;
              }
            }
          }
          std::size_t maximizationHeuristicLevel = _result.heuristicLevel();

          /// Minimize

          _directionVector *= -1;
          _commonValue *= -1;

          _objectiveBound.value = _commonValue;
          _output->onBeforeOracleMinimize();
          oracleMaximize(minHeuristic);
          _output->onAfterOracleMinimize(_result.points.size(), _result.rays.size());

          if (_result.isInfeasible())
          {
            if (!_spanningPoints.empty() || !_spanningRays.empty())
            {
              throw std::runtime_error(
                "AffineHull: Oracle claims infeasible though we know some points/rays.");
            }
            break;
          }
          else if (_result.isUnbounded())
          {
            _output->onBeforeRay();
            addRay(_result.rays.front().vector, _directionColumn, true);
            _output->onAfterRay();
            continue;
          }
          else
          {
            Vector* differentObjectiveVector = NULL;
            for (std::size_t i = 0; i < _result.points.size(); ++i)
            {
              if (_result.points[i].objectiveValue != _commonValue)
              {
                differentObjectiveVector = &_result.points[i].vector;
                break;
              }
            }

            if (differentObjectiveVector != NULL)
            {
              _output->onBeforePoint(false);
              addPoint(*differentObjectiveVector, _directionColumn, true);
              _output->onAfterPoint(false);
              continue;
            }
          }

          std::size_t minimizationHeuristicLevel = _result.heuristicLevel();

          _directionVector *= -1;
          _commonValue *= -1;

          // We have (probably) found an equation.

          LinearConstraint equation = LinearConstraint('=', denseToVector(_directionVector), _commonValue);
          _espace->add(equation.normal());
          _columns[_directionColumn].definesEquation = true;

          if (minimizationHeuristicLevel > 0 || maximizationHeuristicLevel > 0)
          {
            _potentialEquations.push_back(equation);
            _output->onPotentialEquation();
          }
          else
          {
            _irredundantEquations.push_back(equation);
            _output->onEquation();
          }
        }
      }

      int run(std::vector<LinearConstraint>& equations, const std::shared_ptr<OracleBase>& oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification, bool removeRedundantEqns)
      {
        /// Free data from previous run.

        if (_factorization != NULL)
          delete _factorization;

        /// Reset data.

        _irredundantEquations = equations;
        _oracle = oracle;
        _output = &output;
        _spanningPoints.clear();
        _spanningRays.clear();
        _basicColumns.clear();
        _irredundantEquations.clear();
        _potentialEquations.clear();
        _espace->reset(n());
        _factorization = new Factorization(n() + 1);
        _columns.resize(n() + 1);
        _columns[n()].avoid = true;
        _directionColumn = std::numeric_limits<std::size_t>::max();
        _directionVector.reDim(n(), false);
        _approximateDirectionVector.reDim(n(), false);
        _commonValue = infinity;
        _numCacheQueries = 0;
        _numCacheHits = 0;
        _numHeuristicCalls = 0;
        _numOracleCalls = 0;
        _numApproximateDirectionSolves = 0;
        _numExactDirectionSolves = 0;
        _lastDirectionBitsize = 0;
        _maxDirectionBitsize = 0;
        _infeasible = false;

        /// Start.

        _output->_implementation = this;
        _output->onStart();

        if (_irredundantEquations.size() > 0)
        {
          initializeEquations();
          _output->onAddedInitialEquations(_irredundantEquations.size());
        }

        mainLoop(minHeuristicBeforeVerification);

        if (!_infeasible && dimensionSafeUpperBound() > dimensionLowerBound())
        {
          /// Verify k equations with only k+1 oracle calls. If one fails, continue with exact loop.

          _output->onBeforeVerifyDelayed(_potentialEquations.size());

          Rational rhsSum = 0;
          DVectorRational objectiveSum;
          objectiveSum.reDim(n(), true);
          bool failure = false;
          for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
          {
            vectorToDense(_potentialEquations[i].normal(), _directionVector);
            const Rational& rhs = _potentialEquations[i].rhs();
            objectiveSum += _directionVector;
            rhsSum += rhs;

            _directionVector *= -1;

            _output->onBeforeOracleVerify(i);
            _objectiveBound.value = -rhs;
            _objectiveBound.strict = true;
            oracleMaximize(0);
            ++_numOracleCalls;
            _output->onAfterOracleVerify(_result.points.size(), _result.rays.size());

            if (_result.isFeasible())
            {
              if (_result.points.front().objectiveValue != _objectiveBound.value)
              {
                failure = true;
                break;
              }
            }
            else
            {
              failure = true;
              break;
            }
          }

          if (!failure)
          {
            _directionVector = objectiveSum;
            _objectiveBound.value = rhsSum;
            _objectiveBound.strict = true;
            _output->onBeforeOracleVerify(_potentialEquations.size());
            oracleMaximize(0);
            ++_numOracleCalls;
            _output->onAfterOracleVerify(_result.points.size(), _result.rays.size());

            failure = !_result.isFeasible() || _result.points.front().objectiveValue != rhsSum;
          }

          if (failure)
          {
            _output->onAfterVerifyDelayed(0);
          }
          else
          {
            std::size_t numVerified = _potentialEquations.size();
            std::copy(_potentialEquations.begin(), _potentialEquations.end(),
                std::back_inserter(_irredundantEquations));
            _potentialEquations.clear();
            _output->onAfterVerifyDelayed(numVerified);
          }
        }

        if (dimensionSafeUpperBound() > dimensionLowerBound())
        {
          /// Run loop allowing only exact oracle usage.

          for (std::size_t c = 0; c <= n(); ++c)
            _columns[c].definesEquation = false;

          if (!_potentialEquations.empty())
          {
            _potentialEquations.clear();

            /// Reinitialize edeps.

            _espace->reset(n()); // TODO: Could be improved
            for (std::size_t i = 0; i < _irredundantEquations.size(); ++i)
              _espace->addLazy(_irredundantEquations[i].normal());
            _espace->flushLazy();
          }

          mainLoop(0);
        }

        if (!_infeasible && removeRedundantEqns)
          removeRedundantEquations();

        /// End.

        _output->onEnd();
        _output->_implementation = NULL;

        equations = _irredundantEquations;

        return int(_spanningPoints.size() + _spanningRays.size()) - 1;
      }

      friend class InformationBase;
      friend class Result;

    protected:
      std::size_t n() const
      {
        return _oracle->space().dimension();
      }

    protected:
      std::vector<LinearConstraint> _irredundantEquations;
      std::vector<LinearConstraint> _redundantEquations;
      std::vector<LinearConstraint> _potentialEquations;
      std::shared_ptr<OracleBase> _oracle;
      OutputBase* _output;
      std::vector<Vector> _spanningPoints;
      std::vector<Vector> _spanningRays;
      std::vector<std::size_t> _basicColumns;
      VectorSpaceGenerators *_espace;
      Factorization* _factorization;
      std::vector<NonbasicColumn> _columns;
      OracleResult _result;
      std::size_t _directionColumn;
      ObjectiveBound _objectiveBound;
      soplex::DVectorRational _directionVector;
      soplex::DVectorReal _approximateDirectionVector;
      Rational _commonValue;
      std::size_t _numCacheQueries;
      std::size_t _numCacheHits;
      std::size_t _numHeuristicCalls;
      std::size_t _numOracleCalls;
      std::size_t _numApproximateDirectionSolves;
      std::size_t _numExactDirectionSolves;
      std::size_t _lastDirectionBitsize;
      std::size_t _maxDirectionBitsize;
      bool _infeasible;
    };

    Result::Result()
    {
      _implementation = new Implementation();
    }

    Result::~Result()
    {
      assert(hasImplementation());
      delete _implementation;
    }

    int Result::dimension() const
    {
      ensureImplementation();
      int lower = _implementation->dimensionLowerBound();
      int upper = _implementation->dimensionSafeUpperBound();
      if (lower != upper)
        throw std::runtime_error("AffineHull: Dimension not yet determined.");
      return lower;
    }

    int Result::run(std::vector<LinearConstraint>& equations,
      const std::shared_ptr<OracleBase>& oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification, bool removeRedundantEquations)
    {
      return _implementation->run(equations, oracle, output, minHeuristicBeforeVerification, removeRedundantEquations);
    }

    int run(std::vector<LinearConstraint>& equations, const std::shared_ptr<OracleBase>& oracle, OutputBase& output,
      std::size_t minHeuristicBeforeVerification, bool removeRedundantEquations)
    {
      Implementation implementation;
      return implementation.run(equations, oracle, output, minHeuristicBeforeVerification, removeRedundantEquations);
    }

    InformationBase::InformationBase() :
        _implementation(NULL)
    {

    }

    InformationBase::~InformationBase()
    {

    }

    const std::string& InformationBase::oracleName() const
    {
      ensureImplementation();
      return _implementation->_oracle->name();
    }

    std::size_t InformationBase::numVariables() const
    {
      ensureImplementation();
      return _implementation->n();
    }

    int InformationBase::dimensionLowerBound() const
    {
      ensureImplementation();
      return _implementation->dimensionLowerBound();
    }

    int InformationBase::dimensionUnsafeUpperBound() const
    {
      ensureImplementation();
      return _implementation->dimensionUnsafeUpperBound();
    }

    int InformationBase::dimensionSafeUpperBound() const
    {
      ensureImplementation();
      return _implementation->dimensionSafeUpperBound();
    }

    const Vector& InformationBase::spanningPoint(std::size_t i) const
    {
      ensureImplementation();
      return _implementation->_spanningPoints[i];
    }

    std::size_t InformationBase::numSpanningPoints() const
    {
      ensureImplementation();
      return _implementation->_spanningPoints.size();
    }

    const Vector& InformationBase::spanningRay(std::size_t i) const
    {
      ensureImplementation();
      return _implementation->_spanningRays[i];
    }

    std::size_t InformationBase::numSpanningRays() const
    {
      ensureImplementation();
      return _implementation->_spanningRays.size();
    }

    const std::vector<std::size_t>& InformationBase::basicColumns() const
    {
      ensureImplementation();
      return _implementation->_basicColumns;
    }

    const std::vector<LinearConstraint>& InformationBase::irredundantEquations() const
    {
      ensureImplementation();
      return _implementation->_irredundantEquations;
    }

    const std::vector<LinearConstraint>& InformationBase::redundantEquations() const
    {
      ensureImplementation();
      return _implementation->_redundantEquations;
    }

    const std::vector<LinearConstraint>& InformationBase::potentialEquations() const
    {
      ensureImplementation();
      return _implementation->_potentialEquations;
    }

    std::size_t InformationBase::numIrredundantEquations() const
    {
      return irredundantEquations().size();
    }

    std::size_t InformationBase::numRedundantEquations() const
    {
      return redundantEquations().size();
    }

    std::size_t InformationBase::numPotentialEquations() const
    {
      return potentialEquations().size();
    }

    std::size_t InformationBase::numCacheQueries() const
    {
      ensureImplementation();
      return _implementation->_numCacheQueries;
    }

    std::size_t InformationBase::numCacheHits() const
    {
      ensureImplementation();
      return _implementation->_numCacheHits;
    }

    std::size_t InformationBase::numHeuristicCalls() const
    {
      ensureImplementation();
      return _implementation->_numHeuristicCalls;
    }

    std::size_t InformationBase::numOracleCalls() const
    {
      ensureImplementation();
      return _implementation->_numOracleCalls;
    }

    std::size_t InformationBase::numApproximateDirectionSolves() const
    {
      ensureImplementation();
      return _implementation->_numApproximateDirectionSolves;
    }

    std::size_t InformationBase::numExactDirectionSolves() const
    {
      ensureImplementation();
      return _implementation->_numExactDirectionSolves;
    }

    std::size_t InformationBase::lastDirectionBitsize() const
    {
      ensureImplementation();
      return _implementation->_lastDirectionBitsize;
    }

    std::size_t InformationBase::maxDirectionBitsize() const
    {
      ensureImplementation();
      return _implementation->_maxDirectionBitsize;
    }

    bool InformationBase::hasImplementation() const
    {
      return _implementation != NULL;
    }

    void InformationBase::ensureImplementation() const
    {
      if (!hasImplementation())
        throw std::runtime_error("AffineHull: Output is not associated to an implementation.");
    }

    OutputBase::OutputBase()
    {

    }

    OutputBase::~OutputBase()
    {

    }

    bool OutputBase::isRunning() const
    {
      return hasImplementation();
    }

    void OutputBase::onStart()
    {

    }

    void OutputBase::onAddedInitialEquations(std::size_t numAllEquations)
    {

    }

    void OutputBase::onBeforeApproximateDirections()
    {

    }

    void OutputBase::onAfterApproximateDirections(std::size_t numComputed)
    {

    }

    void OutputBase::onBeforeExactDirections()
    {

    }

    void OutputBase::onAfterExactDirections(std::size_t numComputed)
    {

    }

    void OutputBase::onBeforeCache()
    {

    }

    void OutputBase::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {

    }

    void OutputBase::onBeforeOracleZero()
    {

    }

    void OutputBase::onAfterOracleZero(std::size_t numPoints)
    {

    }

    void OutputBase::onBeforeOracleMaximize()
    {

    }

    void OutputBase::onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays)
    {

    }

    void OutputBase::onBeforeOracleMinimize()
    {

    }

    void OutputBase::onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays)
    {

    }

    void OutputBase::onBeforeOracleVerify(std::size_t verifyIndex)
    {

    }

    void OutputBase::onAfterOracleVerify(std::size_t numPoints, std::size_t numRays)
    {

    }

    void OutputBase::onBeforePoint(bool twoPoints)
    {

    }

    void OutputBase::onAfterPoint(bool twoPoints)
    {

    }

    void OutputBase::onBeforeRay()
    {

    }

    void OutputBase::onAfterRay()
    {

    }

    void OutputBase::onPotentialEquation()
    {

    }

    void OutputBase::onEquation()
    {

    }

    void OutputBase::onBeforeVerifyImmediate()
    {

    }

    void OutputBase::onAfterVerifyImmediate(bool success)
    {

    }

    void OutputBase::onBeforeVerifyDelayed(std::size_t numVerifications)
    {

    }

    void OutputBase::onAfterVerifyDelayed(std::size_t numVerifications)
    {

    }

    void OutputBase::onRemovedRedundantEquations(std::size_t numRemoved)
    {

    }

    void OutputBase::onEnd()
    {

    }

    QuietOutput::QuietOutput()
    {

    }

    QuietOutput::~QuietOutput()
    {

    }

    ProgressOutput::ProgressOutput(std::size_t indent)
    {
      _indent.resize(indent, ' ');
      _numVerificationCalls = std::numeric_limits<std::size_t>::max();
      _timeApproxDirections = 0;
      _timeExactDirections = 0;
      _timeCache = 0;
      _timeHeuristics = 0;
      _timeOracles = 0;
      _timeFactorization = 0;
      _timeStarted = 0;
      _lastTime = 0;
      _lastHeuristic = false;
    }

    ProgressOutput::~ProgressOutput()
    {

    }

    double ProgressOutput::timeStamp()
    {
      double time = _timer.time();
      double elapsed = time - _lastTime;
      _lastTime = time;
      return elapsed;
    }

    void ProgressOutput::onProgress()
    {
      std::cout << _indent << "Points: " << numSpanningPoints() << ", Rays: " << numSpanningRays() << ",  "
          << dimensionLowerBound() << " <= dim <= ";
      if (dimensionUnsafeUpperBound() < dimensionSafeUpperBound())
        std::cout << dimensionUnsafeUpperBound() << " (heuristic)";
      else
        std::cout << dimensionSafeUpperBound();
      std::cout << ",  Heuristic calls: " << numHeuristicCalls() << ", Oracle calls: " << numOracleCalls() << ".\n"
          << std::flush;
    }

    void ProgressOutput::onStart()
    {
      _timer.start();
      _lastTime = _timer.time();
      _timeStarted = _lastTime;
    }

    void ProgressOutput::onAddedInitialEquations(std::size_t numAllEquations)
    {
      onProgress();
    }

    void ProgressOutput::onBeforeApproximateDirections()
    {
      timeStamp();
    }

    void ProgressOutput::onAfterApproximateDirections(std::size_t numComputed)
    {
      _timeApproxDirections += timeStamp();
    }

    void ProgressOutput::onBeforeExactDirections()
    {
      timeStamp();
    }

    void ProgressOutput::onAfterExactDirections(std::size_t numComputed)
    {
      _timeExactDirections += timeStamp();
    }

    void ProgressOutput::onBeforeCache()
    {
      timeStamp();
    }

    void ProgressOutput::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      _timeCache += timeStamp();
    }

    void ProgressOutput::onBeforeOracleZero()
    {
      timeStamp();
    }

    void ProgressOutput::onAfterOracleZero(std::size_t numPoints)
    {
//       if (_lastHeuristic)
//         _timeHeuristics += timeStamp();
//       else
//         _timeOracles += timeStamp();
    }

    void ProgressOutput::onBeforeOracleMaximize()
    {
//       _lastHeuristic = !forceOptimal;
      timeStamp();
    }

    void ProgressOutput::onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays)
    {
      if (_lastHeuristic)
        _timeHeuristics += timeStamp();
      else
        _timeOracles += timeStamp();
    }

    void ProgressOutput::onBeforeOracleMinimize()
    {
//       _lastHeuristic = !forceOptimal;
      timeStamp();
    }

    void ProgressOutput::onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays)
    {
      if (_lastHeuristic)
        _timeHeuristics += timeStamp();
      else
        _timeOracles += timeStamp();
    }

    void ProgressOutput::onBeforeOracleVerify(std::size_t verifyIndex)
    {
      std::cout << _indent << "(calling exact oracle " << (verifyIndex + 1) << "/" << _numVerificationCalls << "..."
          << std::flush;
    }

    void ProgressOutput::onAfterOracleVerify(std::size_t numPoints, std::size_t numRays)
    {
      _timeOracles += timeStamp();
      std::cout << " done)" << std::endl;
    }

    void ProgressOutput::onBeforePoint(bool twoPoints)
    {
      timeStamp();
    }

    void ProgressOutput::onAfterPoint(bool twoPoints)
    {
      onProgress();
      _timeFactorization += timeStamp();
    }

    void ProgressOutput::onBeforeRay()
    {
      timeStamp();
    }

    void ProgressOutput::onAfterRay()
    {
      onProgress();
      _timeFactorization += timeStamp();
    }

    void ProgressOutput::onPotentialEquation()
    {
      onProgress();
    }

    void ProgressOutput::onEquation()
    {
      onProgress();
    }

    void ProgressOutput::onBeforeVerifyImmediate()
    {
      _numVerificationCalls = 2;
    }

    void ProgressOutput::onBeforeVerifyDelayed(std::size_t numVerifications)
    {
      _numVerificationCalls = numVerifications + 1;
    }

    void ProgressOutput::onAfterVerifyDelayed(std::size_t numVerifications)
    {
      onProgress();
    }

    void ProgressOutput::onEnd()
    {
      std::cout << _indent << "Max. direction bitsize: " << maxDirectionBitsize() << "\n";
      std::cout << _indent << "Number of invocations:            Heuristics: " << std::setw(6) << numHeuristicCalls()
          << ", Oracles: " << std::setw(6) << numOracleCalls() << ", Cache: " << std::setw(6) << numCacheQueries()
          << ", Directions: " << std::setw(6) << numApproximateDirectionSolves() << " (approx.), " << std::setw(6)
          << numExactDirectionSolves() << " (exact), Factorizations: " << std::setw(6) << (dimensionLowerBound() + 1)
          << "\n";
      std::cout << _indent << "Timings (s): Overall:   " << std::setw(8) << (_timer.time() - _timeStarted)
          << ", Heuristics: " << std::setw(6) << _timeHeuristics << ", Oracles: " << std::setw(6) << _timeOracles
          << ", Cache: " << std::setw(6) << _timeCache << ", Directions: " << std::setw(6) << _timeApproxDirections
          << " (approx.), " << std::setw(6) << _timeExactDirections << " (exact), Factorizations: " << std::setw(6)
          << _timeFactorization << std::endl;
      _timer.stop();
    }

    DebugOutput::DebugOutput() :
        _lastTime(0.0)
    {

    }

    DebugOutput::~DebugOutput()
    {

    }

    void DebugOutput::printStatus()
    {
      std::cout << "n=" << numVariables() << ", |S|=" << numSpanningPoints() << ", |R|=" << numSpanningRays()
          << ", n-eqs=" << dimensionSafeUpperBound() << " (n-unsafe.eqs=" << dimensionUnsafeUpperBound() << "), #calls="
          << numHeuristicCalls() << "/" << numOracleCalls() << ", cache=" << numCacheHits() << "/" << numCacheQueries()
          << ", #dirs=" << numApproximateDirectionSolves() << "/" << numExactDirectionSolves() << ", bitsize="
          << lastDirectionBitsize() << "<=" << maxDirectionBitsize() << ": ";
    }

    void DebugOutput::timeStamp(const std::string& event, bool printNewline)
    {
      double time = _timer.time();
      double elapsed = time - _lastTime;
      _lastTime = time;
      std::map<std::string, double>::iterator iter = _times.find(event);
      if (iter == _times.end())
        _times[event] = elapsed;
      else
        iter->second += elapsed;

      if (printNewline)
        std::cout << std::endl;
    }

    void DebugOutput::onStart()
    {
      printStatus();
      _timer.start();
      _lastTime = _timer.time();
      timeStamp();
    }

    void DebugOutput::onAddedInitialEquations(std::size_t numAllEquations)
    {
      printStatus();
      std::cout << "Added " << numAllEquations << " initial equations.";
      timeStamp();
    }

    void DebugOutput::onBeforeApproximateDirections()
    {
      printStatus();
      std::cout << "Computing approximate directions.";
      timeStamp();
    }

    void DebugOutput::onAfterApproximateDirections(std::size_t numComputed)
    {
      printStatus();
      std::cout << "Computed " << numComputed << " approximate directions.";
      timeStamp("approximate directions");
    }

    void DebugOutput::onBeforeExactDirections()
    {
      printStatus();
      std::cout << "Computing exact directions.";
      timeStamp();
    }

    void DebugOutput::onAfterExactDirections(std::size_t numComputed)
    {
      printStatus();
      std::cout << "Computed " << numComputed << " exact directions.";
      timeStamp("exact directions");
    }

    void DebugOutput::onBeforeCache()
    {
      printStatus();
      std::cout << "Searching cache.";
      timeStamp();
    }

    void DebugOutput::onAfterCache(std::size_t numPoints, std::size_t numRays)
    {
      printStatus();
      if (numPoints > 0)
        std::cout << "Found " << numPoints << " points.";
      else if (numRays > 0)
        std::cout << "Found " << numRays << " rays.";
      else
        std::cout << "Found nothing.";
      timeStamp("cache");
    }

    void DebugOutput::onBeforeOracleZero()
    {
      printStatus();
      std::cout << "Calling oracle with zero objective.";
      timeStamp();
    }

    void DebugOutput::onAfterOracleZero(std::size_t numPoints)
    {
      printStatus();
      if (numPoints > 0)
        std::cout << "Oracle returned " << numPoints << " points.";
      else
        std::cout << "Oracle claims infeasible.";
      timeStamp("oracle-zero");
    }

    void DebugOutput::onBeforeOracleMaximize()
    {
      printStatus();
//       if (forceOptimal)
        std::cout << "Calling oracle for maximization.";
//       else
//         std::cout << "Calling oracle for unsafe maximization.";
      timeStamp();
    }

    void DebugOutput::onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays)
    {
      printStatus();
      if (numPoints > 0)
        std::cout << "Oracle returned " << numPoints << " points.";
      else if (numRays > 0)
        std::cout << "Oracle returned " << numRays << " rays.";
      else
        std::cout << "Oracle claims infeasible.";
      timeStamp("oracle-maximize");
    }

    void DebugOutput::onBeforeOracleMinimize()
    {
      printStatus();
//       if (forceOptimal)
        std::cout << "Calling oracle for minimization.";
//       else
//         std::cout << "Calling oracle for unsafe minimization.";
      timeStamp();
    }

    void DebugOutput::onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays)
    {
      printStatus();
      if (numPoints > 0)
        std::cout << "Oracle returned " << numPoints << " points.";
      else if (numRays > 0)
        std::cout << "Oracle returned " << numRays << " rays.";
      else
        std::cout << "Oracle claims infeasible.";
      timeStamp("oracle-minimize");
    }

    void DebugOutput::onBeforeOracleVerify(std::size_t verifyIndex)
    {
      printStatus();
      std::cout << "Calling oracle for verification.";
      timeStamp();
    }

    void DebugOutput::onAfterOracleVerify(std::size_t numPoints, std::size_t numRays)
    {
      printStatus();
      if (numPoints > 0)
        std::cout << "Oracle returned " << numPoints << " points.";
      else if (numRays > 0)
        std::cout << "Oracle returned " << numRays << " rays.";
      else
        std::cout << "Oracle claims infeasible.";
      timeStamp("oracle-verify");
    }

    void DebugOutput::onBeforePoint(bool twoPoints)
    {
      printStatus();
      if (twoPoints)
        std::cout << "Adding two points.";
      else
        std::cout << "Adding a point.";
      timeStamp();
    }

    void DebugOutput::onAfterPoint(bool twoPoints)
    {
      printStatus();
      if (twoPoints)
        std::cout << "Added two points.";
      else
        std::cout << "Added a point.";
      timeStamp("add-point");
    }

    void DebugOutput::onBeforeRay()
    {
      printStatus();
      std::cout << "Adding a ray.";
      timeStamp();
    }

    void DebugOutput::onAfterRay()
    {
      printStatus();
      std::cout << "Added a ray.";
      timeStamp("add-ray");
    }

    void DebugOutput::onPotentialEquation()
    {
      printStatus();
      std::cout << "Added a potential equation.";
      timeStamp();
    }

    void DebugOutput::onEquation()
    {
      printStatus();
      std::cout << "Added an equation.";
      timeStamp();
    }

    void DebugOutput::onBeforeVerifyImmediate()
    {
      printStatus();
      std::cout << "Verifying a potential equation.";
      timeStamp();
    }

    void DebugOutput::onAfterVerifyImmediate(bool success)
    {
      printStatus();
      if (success)
        std::cout << "Verified potential equation.";
      else
        std::cout << "Potential equation not verified. Switching to exact oracle.";
      timeStamp();
    }

    void DebugOutput::onBeforeVerifyDelayed(std::size_t numVerifications)
    {
      printStatus();
      std::cout << "Starting delayed verification for " << numVerifications << " potential equations.";
      timeStamp();
    }

    void DebugOutput::onAfterVerifyDelayed(std::size_t numVerifications)
    {
      printStatus();
      if (numVerifications > 0)
        std::cout << "Verified " << numVerifications << " potential equations.";
      else
        std::cout << "Verification failed. Switching to exact oracle.";
      timeStamp();
    }

    void DebugOutput::onRemovedRedundantEquations(std::size_t numRemoved)
    {
      printStatus();
      std::cout << "Removed " << numRemoved << " redundant equations.";
      timeStamp();
    }

    void DebugOutput::onEnd()
    {
      for (std::map<std::string, double>::const_iterator iter = _times.begin(); iter != _times.end(); ++iter)
      {
        if (!iter->first.empty())
          std::cout << "Time for " << iter->first << ": " << iter->second << std::endl;
      }
    }

  }

}

