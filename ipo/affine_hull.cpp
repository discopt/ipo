#include "affine_hull.h"

#include <algorithm>

#include "reconstruct.h"
#include "vector_space_generators.h"

using namespace soplex;

//#define RANDOM_COBASIS_INDEX
//#define EXACT_CACHE_SEARCH


namespace ipo {

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

      void computeKernelVector(std::size_t columnIndex, DSVectorRational& result) const
      {
        DVectorRational rhs;
        rhs.reDim(numRows(), true);
        const SVectorRational& col = _columnVectors[columnIndex];
        for (int p = col.size() - 1; p >= 0; --p)
          rhs[col.index(p)] = -col.value(p);

        solveLeftFactor(rhs);
        solveRightFactor(rhs);

        result.clear();
        result.add(columnIndex, Rational(1));
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        {
          if (rhs[b] != 0)
            result.add(_basicColumns[b], rhs[b]);
        }
      }

      void computeApproximateKernelVector(std::size_t columnIndex, DSVectorReal& result, double epsilon) const
      {
        DVectorReal rhs;
        rhs.reDim(numRows(), true);
        const SVectorRational& col = _columnVectors[columnIndex];
        for (int p = col.size() - 1; p >= 0; --p)
          rhs[col.index(p)] = -(double) col.value(p);

        solveLeftApproximateFactor(rhs);
        solveRightApproximateFactor(rhs);

        result.clear();
        result.add(columnIndex, 1.0);
        for (std::size_t b = 0; b < _basicColumns.size(); ++b)
        {
          if (fabs(rhs[b]) > epsilon)
            result.add(_basicColumns[b], rhs[b]);
        }
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

      void addRow(const SVectorRational& row, const Rational& last, std::size_t newBasicColumnIndex)
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
      DSVectorRational exactDirection;
      Rational exactRhs;
      bool approxValid;
      DSVectorReal approxDirection;
      bool definesEquation;
      bool avoid;

      NonbasicColumn()
      {
        exactValid = false;
        approxValid = false;
        definesEquation = false;
        avoid = false;
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

        return approxDirection.size() < other.approxDirection.size();
      }
    };

    class Implementation
    {
    public:
      Implementation()
        : _points(NULL), _directions(NULL), _equations(NULL), _oracle(NULL), _output(NULL),
        _factorization(NULL)
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
      }

      virtual ~Implementation()
      {
        if (_factorization != NULL)
          delete _factorization;
      }

      int dimensionLowerBound() const
      {
        return int(_spanningPoints.size() + _spanningDirections.size()) - 1;
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
        _espace.reset(n());
        for (int e = 0; e < _equations->num(); ++e)
          _espace.addLazy(_equations->rowVector(e));
        _espace.flushLazy();

        _irredundantEquations.clear();
        for (int e = 0; e < _equations->num(); ++e)
        {
          if (!_espace.isDependent(e))
            _irredundantEquations.push_back(e);
        }
      }

      void removePotentialEquations()
      {
#ifdef IPO_DEBUG
        int lastIrredundant = -1;
        for (std::size_t i = 0; i < _irredundantEquations.size(); ++i)
        {
          if (_irredundantEquations[i] > lastIrredundant)
            lastIrredundant = _irredundantEquations[i];
        }
        int firstPotential = _equations->num();
        for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
        {
          if (_potentialEquations[i] < firstPotential)
            firstPotential = _potentialEquations[i];
        }
        assert(firstPotential > lastIrredundant);
#endif

        /// Remove equations.

        std::vector<int> remove;
        for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
          remove.push_back(_potentialEquations[i]);
        _equations->remove(&remove[0], remove.size());

        /// Recompute dependencies.
        initializeEquations();
      }

      void removeRedundantEquations()
      {
        std::vector<int> redundantEquations;
        for (std::size_t e = 0; e < _equations->num(); ++e)
        {
          if (_espace.isDependent(e))
            redundantEquations.push_back(e);
        }
        _equations->remove(&redundantEquations[0], redundantEquations.size());
        _irredundantEquations.clear();
        for (std::size_t e = 0; e < _equations->num(); ++e)
          _irredundantEquations.push_back(e);
        _potentialEquations.clear();
        /// TODO: After this, all equations are irredundant and not potential anymore, so
// approximation status got lost.

        _output->onRemovedRedundantEquations(redundantEquations.size());
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
        _oracle->maximize(_result, _sparseDirectionVector, _objectiveBound,
          std::numeric_limits<std::size_t>::max(), minHeuristic);

        if (_result.isInfeasible())
        {
          _infeasible = true;
          DSVectorRational zeroVector;
          _equations->add(Rational(1), zeroVector, Rational(1));
          _irredundantEquations.clear();
          _irredundantEquations.push_back(_equations->num() - 1);
          _potentialEquations.clear();
          _output->onEquation();
        }
      }

      void addPoint(std::size_t pointIndex, std::size_t column, bool invalidate)
      {
        _spanningPoints.push_back(pointIndex);
        _basicColumns.push_back(column);
        _factorization->addRow(*_points->vector(pointIndex), Rational(-1), column);

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

      void addDirection(std::size_t directionIndex, std::size_t column, bool invalidate)
      {
        _spanningDirections.push_back(directionIndex);
        _basicColumns.push_back(column);
        _factorization->addRow(*_directions->vector(directionIndex), Rational(0), column);

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
        _sparseDirectionVector.clear();

        _output->onBeforeOracleZero();

        _oracle->maximize(_result, _sparseDirectionVector);

        _output->onAfterOracleZero(_result.points.size());

        if (_result.isUnbounded())
        {
          throw std::runtime_error("Oracle claims unbounded for zero objective vector.");
        }
        else if (_result.isInfeasible())
        {
          if (_spanningDirections.empty())
            return;
          else
            throw std::runtime_error("Oracle claims infeasible after returning points or rays.");
        }
        else
        {
          assert(_result.isFeasible());

          _result.addToContainers(*_points, *_directions);

          _output->onBeforePoint(false);
          addPoint(_result.points.front().index, n(), false);
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

          if (_directionColumn == std::numeric_limits<std::size_t>::max() || _columns[c] <
_columns[_directionColumn])
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

        col.exactDirection.sort();
        col.exactRhs = 0;
        if (col.exactDirection.size() > 0)
        {
          int pos = col.exactDirection.size() - 1;
          if (col.exactDirection.index(pos) == n())
          {
            col.exactRhs = col.exactDirection.value(pos);
            col.exactDirection.remove(pos);
          }
        }
        col.exactValid = true;

        /// Also update approximate direction.

        col.approxDirection = col.exactDirection;
        col.approxValid = true;
      }

      bool checkDirectionDepends()
      {
        if (_espace.isDependent(_columns[_directionColumn].exactDirection))
        {
          _columns[_directionColumn].definesEquation = true;
          return true;
        }
        return false;
      }

//       struct PointRayInfo
//       {
//         std::size_t sparsity;
//         std::size_t index;
//
//         PointRayInfo(std::size_t sp, std::size_t idx) :
//             sparsity(sp), index(idx)
//         {
//
//         }
//
//         bool operator<(const PointRayInfo& other) const
//         {
//           return sparsity < other.sparsity;
//         }
//       };
//
//       std::size_t searchCache(UniqueRationalVectorsBase* objects, bool points)
//       {
//         std::vector<PointRayInfo> info;
//         for (std::size_t i = objects->first(); i < objects->size(); i = objects->next(i))
//         {
//           info.push_back(PointRayInfo(objects->nonzeros(i), i));
//         }
//         std::sort(info.begin(), info.end());
//
//         std::size_t bestIndex = std::numeric_limits<std::size_t>::max();
//
//         /// Find the sparsest point/ray with good approx. scalar product.
//
//         double threshold = 1.e-3;
//         Rational localCommonValue = points ? _commonValue : Rational(0);
//         double approximateCommonValue = points ? double(localCommonValue) : 0.0;
//         std::size_t mostDifferent = std::numeric_limits<std::size_t>::max();
//         double mostDifference = 0.0;
//         for (std::size_t i = 0; i < info.size(); ++i)
//         {
//           std::size_t index = info[i].index;
// #ifndef EXACT_CACHE_SEARCH
//           double approxProduct = *objects->approximation(index) * _denseApproximateDirectionVector;
//           approxProduct = fabs(approxProduct - approximateCommonValue);
//           if (approxProduct > mostDifference)
//           {
//             mostDifferent = index;
//             mostDifference = approxProduct;
//           }
//
//           if (approxProduct > threshold)
// #endif
//           {
//             Rational product = *objects->vector(index) * _denseDirectionVector;
//             if (product != localCommonValue)
//               return index;
//             else
//               threshold *= 10;
//           }
//         }
//
//         /// If no attempt succeeds, take the one with the largest approx. discrepancy.
//
//         if (mostDifference > 0.0)
//         {
//           Rational product = *objects->vector(mostDifferent) * _denseDirectionVector;
//           if (product != localCommonValue)
//           {
//             return mostDifferent;
//           }
//         }
//
//         return std::numeric_limits<std::size_t>::max();
//       }
//
//       bool searchPoints()
//       {
//         std::size_t bestIndex = std::numeric_limits<std::size_t>::max();
//
//         if (_spanningPoints.empty())
//         {
//           std::size_t index = _points->first();
//           if (index < _points->size())
//             bestIndex = index;
//         }
//         else
//           bestIndex = searchCache(_points, true);
//
//         if (bestIndex < std::numeric_limits<std::size_t>::max())
//         {
//           _output->onAfterCache(1, 0);
//           _output->onBeforePoint(false);
//           addPoint(bestIndex, _spanningPoints.empty() ? n() : _directionColumn, true);
//           _output->onAfterPoint(false);
//           return true;
//         }
//         else
//           return false;
//       }
//
//       bool searchRays()
//       {
//         std::size_t bestIndex = searchCache(_directions, false);
//
//         if (bestIndex < std::numeric_limits<std::size_t>::max())
//         {
//           _output->onAfterCache(0, 1);
//           _output->onBeforeDirection();
//           addDirection(bestIndex, _directionColumn, true);
//           _output->onAfterDirection();
//           return true;
//         }
//         else
//           return false;
//       }

      std::size_t findPointsDifferingVariable(std::size_t firstPointIndex, std::size_t
secondPointIndex)
      {
        const DSVectorRational& firstPoint = *_points->vector(firstPointIndex);
        const DSVectorRational& secondPoint = *_points->vector(secondPointIndex);

        DVectorRational difference;
        difference.reDim(n(), true);
        difference.assign(firstPoint);
        for (int p = secondPoint.size() - 1; p >= 0; --p)
          difference[secondPoint.index(p)] -= secondPoint.value(p);
        for (std::size_t v = 0; v < n(); ++v)
        {
          if (_factorization->isBasic(v))
            continue;
          if (difference[v] != 0)
          {
            return v;
          }
        }

        return std::numeric_limits<std::size_t>::max();
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
            && _spanningDirections.size() == getUpperBound(minHeuristic))
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
          _sparseDirectionVector = _columns[_directionColumn].exactDirection;
          _denseDirectionVector.clear();
          _denseDirectionVector.assign(_sparseDirectionVector);
          for (std::size_t v = 0; v < n(); ++v)
            _denseApproximateDirectionVector[v] = double(_denseDirectionVector[v]);
          if (_spanningPoints.empty())
            _commonValue = infinity;
          else
          {
            _commonValue = *_points->vector(*_spanningPoints.begin()) * _denseDirectionVector;
          }

          /// Measure bitsize.

          _lastDirectionBitsize = 0;
          for (int p = _sparseDirectionVector.size() - 1; p >= 0; --p)
            _lastDirectionBitsize += _sparseDirectionVector.value(p).sizeInBase(2);
          _maxDirectionBitsize = std::max(_maxDirectionBitsize, _lastDirectionBitsize);
          _output->onAfterExactDirections(_numExactDirectionSolves - oldNumExactDirectionSolves);

          /*
           * TODO: To avoid oracle calls of low heuristic levels we might maximize/minimize
           * specific levels taking turns.
           */

//           std::cout << "Objective: ";
//           _oracle->space().printLinearForm(std::cout, &_sparseDirectionVector);
//           std::cout << std::endl;

          /// Maximize direction.

          _objectiveBound.value = _commonValue;
          _objectiveBound.strict = true;
          _output->onBeforeOracleMaximize();
          oracleMaximize(minHeuristic);
          _output->onAfterOracleMaximize(_result.points.size(), _result.directions.size());

//           std::cout << "max ";
//           if (_result.isFeasible())
//             std::cout << "= " << _result.points.front().objectiveValue << std::endl;
//           else if (_result.isUnbounded())
//             std::cout << " is unbounded." << std::endl;
//           else
//             std::cout << " is infeasible." << std::endl;

          _result.addToContainers(*_points, *_directions);

          if (_result.isInfeasible())
          {
            if (!_spanningPoints.empty() || !_spanningDirections.empty())
            {
              throw std::runtime_error(
                "AffineHull: Oracle claims infeasible though we know some points/directions.");
            }
            break;
          }
          else if (_result.isUnbounded())
          {
            _output->onBeforeDirection();
            addDirection(_result.directions.front().index, _directionColumn, true);
            _output->onAfterDirection();
            continue;
          }
          else
          {
            assert(_result.isFeasible());

            if (_spanningPoints.empty())
            {
              // Case |S| = 0: We might add two points.

              const Rational& firstObjectiveValue = _result.points.front().objectiveValue;
              std::size_t secondPoint = 0;
              for (std::size_t i = 1; i < _result.points.size(); ++i)
              {
                if (_result.points[i].objectiveValue != firstObjectiveValue)
                {
                  secondPoint = i;
                  break;
                }
              }
              std::size_t firstIndex = _result.points.front().index;

              // We always add the maximizer.

              _output->onBeforePoint(secondPoint > 0);
              addPoint(firstIndex, n(), true);
              if (secondPoint > 0)
              {
                // If there exist a point with different objective value, we add it as well.

                 std::size_t secondIndex = _result.points[secondPoint].index;
                 std::size_t differingVariable = findPointsDifferingVariable(firstIndex,
                   secondIndex);
                 assert(differingVariable != std::numeric_limits<std::size_t>::max());
                 addPoint(secondIndex, differingVariable, true);
              }
              _output->onAfterPoint(secondPoint > 0);
              if (secondPoint > 0)
                continue;
              else
                _commonValue = _result.points.front().objectiveValue;
            }
            else
            {
              // We add any point that has objective different from common value.

              std::size_t index = std::numeric_limits<std::size_t>::max();
              for (std::size_t i = 0; i < _result.points.size(); ++i)
              {
                if (_result.points[i].objectiveValue != _commonValue)
                {
                  index = _result.points[i].index;
                  break;
                }
              }
              if (index != std::numeric_limits<std::size_t>::max())
              {
                _output->onBeforePoint(false);
                addPoint(index, _directionColumn, true);
                _output->onAfterPoint(false);
                continue;
              }
            }
          }
          std::size_t maximizationHeuristicLevel = _result.heuristic();

          /// Minimize

          _sparseDirectionVector *= -1;
          _commonValue *= -1;

          _output->onBeforeOracleMinimize();
          oracleMaximize(minHeuristic);
          _output->onAfterOracleMinimize(_result.points.size(), _result.directions.size());

//           std::cout << "min ";
//           if (_result.isFeasible())
//             std::cout << "= " << -_result.points.front().objectiveValue << std::endl;
//           else if (_result.isUnbounded())
//             std::cout << " is unbounded." << std::endl;
//           else
//             std::cout << " is infeasible." << std::endl;

          _result.addToContainers(*_points, *_directions);

          if (_result.isInfeasible())
          {
            if (!_spanningPoints.empty() || !_spanningDirections.empty())
            {
              throw std::runtime_error(
                "AffineHull: Oracle claims infeasible though we know some points/directions.");
            }
            break;
          }
          else if (_result.isUnbounded())
          {
            _output->onBeforeDirection();
            addDirection(_result.points.front().index, _directionColumn, true);
            _output->onAfterDirection();
            continue;
          }
          else
          {
            std::size_t index = std::numeric_limits<std::size_t>::max();
            for (std::size_t i = 0; i < _result.points.size(); ++i)
            {
              if (_result.points[i].objectiveValue > _commonValue)
              {
                index = _result.points[i].index;
                break;
              }
            }
            if (index != std::numeric_limits<std::size_t>::max())
            {
              _output->onBeforePoint(false);
              addPoint(index, _directionColumn, true);
              _output->onAfterPoint(false);
              continue;
            }
          }

          std::size_t minimizationHeuristicLevel = _result.heuristic();

          _sparseDirectionVector *= -1;
          _commonValue *= -1;

          _espace.add(_sparseDirectionVector, false);
          _equations->add(_commonValue, _sparseDirectionVector, _commonValue);
          _columns[_directionColumn].definesEquation = true;

          if (minimizationHeuristicLevel > 0 || maximizationHeuristicLevel > 0)
          {
            _potentialEquations.push_back(_equations->num() - 1);
            _output->onPotentialEquation();
          }
          else
          {
            _irredundantEquations.push_back(_equations->num() - 1);
            _output->onEquation();
          }
        }
      }

      int run(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions,
        LPRowSetRational& equations, OracleBase* oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification, bool removeRedundantEqns)
      {
        /// Free data from previous run.

        if (_factorization != NULL)
          delete _factorization;

        /// Reset data.

        _points = &points;
        _directions = &directions;
        _equations = &equations;
        _oracle = oracle;
        _output = &output;
        _spanningPoints.clear();
        _spanningDirections.clear();
        _basicColumns.clear();
        _irredundantEquations.clear();
        _potentialEquations.clear();
        _espace.reset(n());
        _factorization = new Factorization(n() + 1);
        _columns.resize(n() + 1);
        _columns[n()].avoid = true;
        _directionColumn = std::numeric_limits<std::size_t>::max();
        _denseDirectionVector.reDim(n(), false);
        _denseApproximateDirectionVector.reDim(n(), false);
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

        if (_equations->num() > 0)
        {
          initializeEquations();
          _output->onAddedInitialEquations(_equations->num());
        }

        mainLoop(minHeuristicBeforeVerification);

//         if (!_infeasible && dimensionSafeUpperBound() > dimensionLowerBound())
//         {
//           /// Verify k equations with only k+1 oracle calls. If one fails, continue with exact
// loop.
//
//           _output->onBeforeVerifyDelayed(_potentialEquations.size());
//
//           DSVectorRational objective;
//           Rational rhsSum = 0;
//           DVectorRational objectiveSum;
//           objectiveSum.reDim(n(), true);
//           bool failure = false;
//           for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
//           {
//             objective = _equations->rowVector(_potentialEquations[i]);
//             objectiveSum += objective;
//
//             objective *= -1;
//             rhsSum += _equations->rhs(_potentialEquations[i]);
//
//             _output->onBeforeOracleVerify(i);
//             _oracle->maximize(_result, objective, true);
//             ++_numOracleCalls;
//             _output->onAfterOracleVerify(_result.points.size(), _result.directions.size());
//
//             if (_result.isFeasible())
//             {
//
//               for (std::size_t j = 0; j < _result.points.size(); ++j)
//                 _points->insertFree(_result.points[j]);
//
//               if (_result.bestValue != -_equations->rhs(_potentialEquations[i]))
//               {
//                 failure = true;
//                 break;
//               }
//             }
//             else
//             {
//               for (std::size_t j = 0; j < _result.directions.size(); ++j)
//                 _directions->insertFree(_result.directions[j]);
//               failure = true;
//               break;
//             }
//           }
//
//           if (!failure)
//           {
//             _output->onBeforeOracleVerify(_potentialEquations.size());
//             _oracle->maximize(_result, objectiveSum, true);
//             ++_numOracleCalls;
//             _output->onAfterOracleVerify(_result.points.size(), _result.directions.size());
//
//             for (std::size_t j = 0; j < _result.points.size(); ++j)
//               _points->insertFree(_result.points[j]);
//             for (std::size_t j = 0; j < _result.directions.size(); ++j)
//               _directions->insertFree(_result.directions[j]);
//
//             failure = !_result.isFeasible() || _result.bestValue != rhsSum;
//           }
//
//           if (failure)
//           {
//             _output->onAfterVerifyDelayed(0);
//           }
//           else
//           {
//             std::size_t numVerified = _potentialEquations.size();
//             std::copy(_potentialEquations.begin(), _potentialEquations.end(),
//                 std::back_inserter(_irredundantEquations));
//             _potentialEquations.clear();
//             _output->onAfterVerifyDelayed(numVerified);
//           }
//         }
//
//         if (optionOracleUsage() != ORACLE_NEVER && dimensionSafeUpperBound() >
// dimensionLowerBound())
//         {
//           /// Run loop allowing only exact oracle usage.
//
//           for (std::size_t c = 0; c <= n(); ++c)
//             _columns[c].definesEquation = false;
//
//           if (!_potentialEquations.empty())
//           {
//             std::vector<int> remove;
//             for (std::size_t i = 0; i < _potentialEquations.size(); ++i)
//               remove.push_back(_potentialEquations[i]);
//             _equations->remove(&remove[0], remove.size());
//             _potentialEquations.clear();
//
//             /// Reinitialize edeps.
//
//             _espace.reset(n()); // TODO: Could be improved
//             for (int e = 0; e < _equations->num(); ++e)
//               _espace.addLazy(_equations->rowVector(e));
//             _espace.flushLazy();
//           }
//
//           mainLoop(false, false);
//         }

        if (!_infeasible && removeRedundantEqns)
          removeRedundantEquations();

        /// End.

        _output->onEnd();
        _output->_implementation = NULL;

        return int(_spanningPoints.size() + _spanningDirections.size()) - 1;
      }

      friend class InformationBase;
      friend class Result;

    protected:
      std::size_t n() const
      {
        return _oracle->space().dimension();
      }

    protected:
      UniqueRationalVectorsBase* _points;
      UniqueRationalVectorsBase* _directions;
      LPRowSetRational* _equations;
      OracleBase* _oracle;
      OutputBase* _output;
      VectorSubset _spanningPoints;
      VectorSubset _spanningDirections;
      VectorSubset _basicColumns;
      VectorSubset _irredundantEquations;
      VectorSubset _potentialEquations;
      VectorSpaceGenerators _espace;
      Factorization* _factorization;
      std::vector<NonbasicColumn> _columns;
      OracleResult _result;
      std::size_t _directionColumn;
      ObjectiveBound _objectiveBound;
      DSVectorRational _sparseDirectionVector;
      DVectorRational _denseDirectionVector;
      DVectorReal _denseApproximateDirectionVector;
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

    int Result::run(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions,
        LPRowSetRational& equations, OracleBase* oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification, bool removeRedundantEquations)
    {
      return _implementation->run(points, directions, equations, oracle, output,
        minHeuristicBeforeVerification, removeRedundantEquations);
    }

    int run(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions,
      LPRowSetRational& equations, OracleBase* oracle, OutputBase& output,
      std::size_t minHeuristicBeforeVerification, bool removeRedundantEquations)
    {
      Implementation implementation;
      return implementation.run(points, directions, equations, oracle, output,
        minHeuristicBeforeVerification, removeRedundantEquations);
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

    const VectorSubset& InformationBase::spanningPoints() const
    {
      ensureImplementation();
      return _implementation->_spanningPoints;
    }

    std::size_t InformationBase::numSpanningPoints() const
    {
      return spanningPoints().size();
    }

    const VectorSubset& InformationBase::spanningDirections() const
    {
      ensureImplementation();
      return _implementation->_spanningDirections;
    }

    std::size_t InformationBase::numSpanningDirections() const
    {
      return spanningDirections().size();
    }

    const VectorSubset& InformationBase::basicColumns() const
    {
      ensureImplementation();
      return _implementation->_basicColumns;
    }

    const VectorSubset& InformationBase::irredundantEquations() const
    {
      ensureImplementation();
      return _implementation->_irredundantEquations;
    }

    std::size_t InformationBase::numIrredundantEquations() const
    {
      return irredundantEquations().size();
    }

    const VectorSubset& InformationBase::potentialEquations() const
    {
      return _implementation->_potentialEquations;
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

    void OutputBase::onBeforeDirection()
    {

    }

    void OutputBase::onAfterDirection()
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
      std::cout << _indent << "Points: " << numSpanningPoints() << ", Rays: " << numSpanningDirections() << ",  "
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

    void ProgressOutput::onBeforeDirection()
    {
      timeStamp();
    }

    void ProgressOutput::onAfterDirection()
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
      std::cout << "n=" << numVariables() << ", |S|=" << numSpanningPoints() << ", |R|=" << numSpanningDirections()
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

    void DebugOutput::onBeforeDirection()
    {
      printStatus();
      std::cout << "Adding a ray.";
      timeStamp();
    }

    void DebugOutput::onAfterDirection()
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

