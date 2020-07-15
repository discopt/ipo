#include <ipo/affine_hull.hpp>

#include "lu.hpp"
#include "redundancy.hpp"

#include <iostream>

namespace ipo
{
  AffineHullQuery::AffineHullQuery()
    : epsilonConstraints(1.0e-9), epsilonSafety(1.0e-6),
    timeLimit(std::numeric_limits<double>::infinity())
  {

  }
  
  template <typename T, typename IsZero>
  class AffineComplement
  {
  public:
    AffineComplement(std::size_t numVariables, const IsZero isZero = IsZero())
      : _numVariables(numVariables), _isZero(isZero), _lu(isZero),
      _columns(numVariables + 1), _lastColumn(0)
    {
      for (ColumnData& columnData : _columns)
      {
        columnData.basisIndex = std::numeric_limits<std::size_t>::max();
        columnData.definesEquation = false;
      }
    }

    inline std::size_t numVariables() const
    {
      return _numVariables;
    }

    inline std::size_t rank() const
    {
      return _basisIndexToColumn.size();
    }

    void add(const sparse_vector<T>& row, const T& last, std::size_t newBasicColumn)
    {
//       std::cout << "Adding row " << row << " with last = " << last << " for column " << newBasicColumn << std::endl;

      // Add column to basis.
      assert(newBasicColumn < _columns.size());
      assert(_columns[newBasicColumn].basisIndex == std::numeric_limits<std::size_t>::max());
      _columns[newBasicColumn].basisIndex = _basisIndexToColumn.size();
      _basisIndexToColumn.push_back(newBasicColumn);

      // Copy part of new row to basis matrix.
      std::vector<T> newRow(rank(), T(0));
      for (const auto& iter : row)
      {
        std::size_t basisIndex = _columns[iter.first].basisIndex;
        if (basisIndex < std::numeric_limits<std::size_t>::max())
          newRow[basisIndex] = iter.second;
      }
      if (last != 0 && _columns.back().basisIndex != std::numeric_limits<std::size_t>::max())
        newRow[_columns.back().basisIndex] = last;

      // Copy new basic column.
      std::vector<T> newColumn(rank() - 1, T(0));
      for (std::size_t i = 0; i < _columns[newBasicColumn].rows.size(); ++i)
        newColumn[_columns[newBasicColumn].rows[i]] = _columns[newBasicColumn].entries[i];

      // Update LU decomposition.
      _lu.extend(&newRow[0], &newColumn[0], newRow.back());

      // Add the row.
      for (const auto& iter : row)
      {
        ColumnData& colulmnData = _columns[iter.first];
        colulmnData.rows.push_back(_basisIndexToColumn.size() - 1);
        colulmnData.entries.push_back(T(0));
        colulmnData.entries.back() = iter.second;
      }

      // Invalidate equations.
      for (auto& columnData : _columns)
        columnData.definesEquation = false;
    }

    void computeKernelVector(std::size_t column, sparse_vector<T>& vector)
    {
      std::vector<T> rhs(rank(), T(0));
      const ColumnData& columnData = _columns[column];
      for (std::size_t i = 0; i < columnData.rows.size(); ++i)
        rhs[columnData.rows[i]] = -columnData.entries[i];

      _lu.solveLeft(&rhs[0]);
      _lu.solveUpper(&rhs[0]);

      vector.clear();
      for (std::size_t c = 0; c < numVariables(); ++c)
      {
        if (c == column)
          vector.push_back(c, T(1));
        else if (_columns[c].basisIndex < std::numeric_limits<std::size_t>::max()
          && !_isZero(rhs[_columns[c].basisIndex]))
        {
          vector.push_back(c, rhs[_columns[c].basisIndex]);
        }
      }
    }

    std::size_t selectColumn()
    {
      std::size_t bestColumn = std::numeric_limits<std::size_t>::max();
      std::size_t c = (_lastColumn + 1) % (_columns.size() - 1);
      while (c != _lastColumn)
      {
        if (_columns[c].basisIndex == std::numeric_limits<std::size_t>::max() 
          && !_columns[c].definesEquation && (bestColumn == std::numeric_limits<std::size_t>::max() 
            || _columns[c].entries.size() < _columns[bestColumn].entries.size()))
        {
          bestColumn = c;
        }
        c = (c+1) % (_columns.size() - 1);
      }
      if (bestColumn == std::numeric_limits<std::size_t>::max())
        bestColumn = _columns.size() - 1;

      _lastColumn = bestColumn;
      return bestColumn;
    }

    void markEquation(std::size_t column)
    {
      _columns[column].definesEquation = true;
    }

  protected:
    struct ColumnData
    {
      std::size_t basisIndex;
      std::vector<std::size_t> rows;
      std::vector<T> entries;
      bool definesEquation;
    };

    std::size_t _numVariables;
    IsZero _isZero;
    IncrementalLUFactorization<T, IsZero> _lu;
    std::vector<ColumnData> _columns;
    std::vector<std::size_t> _basisIndexToColumn;

    /**
     * \brief The last selected column.
     *
     * The last selected column. A kernel vector for a column may change if a new point or ray is
     * added. To save time for checking whether the selected kernel vector lies in the span of the
     * equations we prefer columns after the last returned one.
     */
    std::size_t _lastColumn;
  };

  static
  double elapsedTime(const std::chrono::time_point<std::chrono::system_clock>& started)
  {
    return std::chrono::duration<double>(std::chrono::duration<double>(
      std::chrono::system_clock::now() - started)).count();
  }

  template <typename T, typename IsZero>
  void affineHullImplementation(std::shared_ptr<Polyhedron<T, IsZero>> polyhedron,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultRays,
    std::vector<Constraint<T>>& resultEquations, const AffineHullQuery& query,
    const std::vector<Constraint<T>>& knownEquations, IsZero isZero, AffineHullResultCommon& result)
  {
    auto timeStarted = std::chrono::system_clock::now();
    auto timeComponent = std::chrono::system_clock::now();
    std::size_t n = polyhedron->space()->dimension();
    result.upperBound = n;

    auto redundancyCheck = EquationRedundancyCheck<T, IsZero>(n, isZero);
    for (auto equation : knownEquations)
    {
      if (redundancyCheck.add(equation) == EQUATION_INCONSISTENT)
      {
        assert(resultEquations.empty());
        resultEquations.push_back(neverSatisfiedConstraint<T>());
        result.upperBound = -1;
        return;
      }
    }

    std::cout << "Initial " << knownEquations.size() << " equations have rank "
      << redundancyCheck.rank() << "." << std::endl;
    result.upperBound = n - redundancyCheck.rank();

    auto affineComplement = AffineComplement<T, IsZero>(n, isZero);
    std::vector<T> objective(n);
    while (result.lowerBound < result.upperBound)
    {
      std::cout << "Iteration of affine hull computation: " << result.lowerBound
        << " <= dim <= " << result.upperBound << "." << std::endl;

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      // If we have dim(P) rays but no points yet, we solve the feasibility problem.

      if (resultPoints.empty() && int(resultRays.size()) == result.upperBound)
      {
        typename OptimizationOracle<T>::Query oracleQuery;
        for (std::size_t v = 0; v < n; ++v)
          objective[v] = 0;
        oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
        timeComponent = std::chrono::system_clock::now();
        typename OptimizationOracle<T>::Result oracleResult = polyhedron->maximize(
          &objective[0], oracleQuery);
        result.timeOracles += elapsedTime(timeComponent);

        if (oracleResult.isUnbounded())
          throw std::runtime_error("maximize(<zero vector>) returned unbounded result.");
        else if (oracleResult.isFeasible())
        {
          std::cout << "  -> adding last point" << std::endl;
          resultPoints.push_back(oracleResult.points.front().vector);
          result.dimension = result.lowerBound = result.upperBound;
          return;
        }
        else
          throw std::runtime_error("maximize(<zero vector>) returned infeasible after some rays.");
      }

      sparse_vector<T> kernelDirectionVector;
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      while (kernelDirectionVector.empty())
      {
        timeComponent = std::chrono::system_clock::now();
        kernelDirectionColumn = affineComplement.selectColumn();
        std::cout << "  Column " << kernelDirectionColumn << std::flush;

        affineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector);
        result.timeKernel += elapsedTime(timeComponent);
        assert(!kernelDirectionVector.empty());

        std::cout << " has kernel vector " << kernelDirectionVector << std::flush;

        if (elapsedTime(timeStarted) >= query.timeLimit)
        {
          result.dimension = AFFINEHULL_ERROR_TIMEOUT;
          return;
        }

        timeComponent = std::chrono::system_clock::now();
        EquationRedundancy redundancy = redundancyCheck.test(kernelDirectionVector);
        result.timeEquations += elapsedTime(timeComponent);
  
        if (redundancy == EQUATION_INDEPENDENT)
        {
          std::cout << " which is independent of the equations." << std::endl;
          break;
        }
        else
        {
          std::cout << " which is in the span of the equations." << std::endl;
          affineComplement.markEquation(kernelDirectionColumn);
          kernelDirectionVector.clear();
        }

        if (elapsedTime(timeStarted) >= query.timeLimit)
        {
          result.dimension = AFFINEHULL_ERROR_TIMEOUT;
          return;
        }
      }

      // Prepare for maximize kernelDirectionVector.

      std::cout << "  Maximization" << std::flush;
      for (std::size_t v = 0; v < n; ++v)
        objective[v] = 0;
      for (const auto& iter : kernelDirectionVector)
        objective[iter.first] = iter.second;

      typename OptimizationOracle<T>::Query oracleQuery;
      if (!resultPoints.empty())
      {
        oracleQuery.minObjectiveValue = *resultPoints.front() * kernelDirectionVector;
        std::cout << " with common objective value " << (double)oracleQuery.minObjectiveValue;
      }
      oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
      if (oracleQuery.timeLimit <= 0)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      timeComponent = std::chrono::system_clock::now();
      typename OptimizationOracle<T>::Result oracleResult = polyhedron->maximize(
        &objective[0], oracleQuery);
      result.timeOracles += elapsedTime(timeComponent);

      std::cout << " yields " << oracleResult.points.size() << " points and "
        << oracleResult.rays.size() << " rays with optimum in [" << oracleResult.primalBound << ","
        << oracleResult.dualBound << "]" << std::endl;
      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      if (oracleResult.isInfeasible())
      {
        resultEquations.clear();
        resultEquations.push_back(neverSatisfiedConstraint<T>());
        std::cout << "  -> infeasible." << std::endl;
      }
      else if (oracleResult.isUnbounded())
      {
        resultRays.push_back(oracleResult.rays.front().vector);
        timeComponent = std::chrono::system_clock::now();
        affineComplement.add(*oracleResult.rays.front().vector, 0, kernelDirectionColumn);
        result.timePointsRays += elapsedTime(timeComponent);
        std::cout << "  -> adding a ray." << std::endl;
        ++result.lowerBound;
      }
      else if (resultPoints.empty())
      {
        assert(oracleResult.isFeasible());

        // Add the maximizer.

        resultPoints.push_back(oracleResult.points.front().vector);
        timeComponent = std::chrono::system_clock::now();
        affineComplement.add(*resultPoints.back(), 1, n);
        result.timePointsRays += elapsedTime(timeComponent);
        ++result.lowerBound;

        // If another point has a different objective value, we also add that and continue.
        for (std::size_t p = 2; p < oracleResult.points.size(); ++p)
        {
          if (!isZero(oracleResult.points[p].objectiveValue - oracleResult.points.front().objectiveValue))
          {
            resultPoints.push_back(oracleResult.points[p].vector);
            timeComponent = std::chrono::system_clock::now();
            affineComplement.add(*resultPoints.back(), 1, kernelDirectionColumn);
            result.timePointsRays += elapsedTime(timeComponent);
            ++result.lowerBound;
            std::cout << "  -> adding two points with objective values "
              << (double)oracleResult.points.front().objectiveValue << " and " 
              << (double)oracleResult.points[p].objectiveValue << "." << std::endl;
            break;
          }
        }
        if (resultPoints.size() == 2)
          continue;

        oracleQuery.minObjectiveValue = oracleResult.points.front().objectiveValue;
        std::cout << "  -> adding a point with objective value "
          << (double)oracleResult.points.front().objectiveValue << std::endl;
      }
      else
      {
        assert(oracleResult.isFeasible());

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : oracleResult.points)
        {
          if (!isZero(point.objectiveValue - oracleQuery.minObjectiveValue))
          {
            resultPoints.push_back(point.vector);
            timeComponent = std::chrono::system_clock::now();
            affineComplement.add(*point.vector, 1, kernelDirectionColumn);
            result.timePointsRays += elapsedTime(timeComponent);
            ++result.lowerBound;
            std::cout << "  -> adding a point with objective value "
              << (double)oracleResult.points.front().objectiveValue << std::endl;
            success = true;
            break;
          }
        }
        if (success)
          continue;
      }

      // Prepare minimization.

      for (std::size_t v = 0; v < n; ++v)
        objective[v] *= -1;
      oracleQuery.minObjectiveValue *= -1;

      oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
      if (oracleQuery.timeLimit <= 0)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      timeComponent = std::chrono::system_clock::now();
      oracleResult = polyhedron->maximize(&objective[0], oracleQuery);
      result.timeOracles += elapsedTime(timeComponent);

      std::cout << "  Minimization yields " << oracleResult.points.size() << " points and "
        << oracleResult.rays.size() << " rays with optimum in [" << oracleResult.primalBound << ","
        << oracleResult.dualBound << "]" << std::endl;

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      if (oracleResult.isUnbounded())
      {
        resultRays.push_back(oracleResult.rays.front().vector);

        timeComponent = std::chrono::system_clock::now();
        affineComplement.add(*oracleResult.rays.front().vector, 0, kernelDirectionColumn);
        result.timePointsRays += elapsedTime(timeComponent);

        ++result.lowerBound;
        std::cout << "  -> adding a ray." << std::endl;
      }
      else
      {
        assert(oracleResult.isFeasible());

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : oracleResult.points)
        {
          if (!isZero(point.objectiveValue - oracleQuery.minObjectiveValue))
          {
            resultPoints.push_back(point.vector);

            timeComponent = std::chrono::system_clock::now();
            affineComplement.add(*point.vector, 1, kernelDirectionColumn);
            result.timePointsRays += elapsedTime(timeComponent);
            ++result.lowerBound;
            std::cout << "  -> adding a point with objective value "
              << (double)oracleResult.points.front().objectiveValue << std::endl;
            success = true;
            break;
          }
        }
        if (success)
          continue;
      }
      
      auto vector = std::make_shared<sparse_vector<T>>(std::move(kernelDirectionVector));

      // We have to add an equation.
      resultEquations.push_back(Constraint<T>(-oracleQuery.minObjectiveValue, vector,
        -oracleQuery.minObjectiveValue));
      affineComplement.markEquation(kernelDirectionColumn);

      timeComponent = std::chrono::system_clock::now();
      redundancyCheck.add(resultEquations.back());
      result.timeEquations += elapsedTime(timeComponent);

      --result.upperBound;
      std::cout << "  -> adding an equation." << std::endl;
    }

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
  }

  AffineHullResult<double> affineHull(
    std::shared_ptr<Polyhedron<double, DoubleIsZero>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<double>>& knownEquations)
  {
    AffineHullResult<double> result;
    affineHullImplementation(polyhedron, result.points, result.rays, result.equations, query,
      knownEquations, DoubleIsZero(1.0e-9), result);
    return result;
  }

#if defined(IPO_WITH_GMP)

  AffineHullResult<rational> affineHull(
    std::shared_ptr<Polyhedron<rational, RationalIsZero>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<rational>>& knownEquations)
  {
    AffineHullResult<rational> result;
    affineHullImplementation(polyhedron, result.points, result.rays, result.equations, query,
      knownEquations, RationalIsZero(), result);
    return result;
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
