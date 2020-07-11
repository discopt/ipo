#include <ipo/ipo.hpp>

#include "lu.hpp"
#include "redundancy.hpp"

#include <iostream>

namespace ipo
{
  template <typename T, typename IsZero>
  class AffineComplement
  {
  public:
    AffineComplement(std::size_t numVariables, const IsZero isZero = IsZero())
      : _numVariables(numVariables), _isZero(isZero), _lu(isZero),
      _columns(numVariables + 1)
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
      for (std::size_t c = 0; c < _columns.size(); ++c)
      {
        if (_columns[c].basisIndex < std::numeric_limits<std::size_t>::max() || _columns[c].definesEquation)
          continue;

        if (bestColumn == std::numeric_limits<std::size_t>::max())
          bestColumn = c;
        else if (c == _columns.size() - 1)
          break;
        else if (_columns[c].entries.size() < _columns[bestColumn].entries.size())
          bestColumn = c;
      }
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
  };

  static
  double remainingTime(const std::chrono::time_point<std::chrono::system_clock>& started, double limit)
  {
    return limit - std::chrono::duration<double>(std::chrono::duration<double>(
      std::chrono::system_clock::now() - started)).count();
  }

  template <typename T, typename IsZero>
  int affineHullImplementation(std::shared_ptr<Polyhedron<T, IsZero>> polyhedron,
    std::vector<sparse_vector<T>>& innerPoints,
    std::vector<sparse_vector<T>>& innerRays,
    std::vector<Constraint<T>>& outerEquations,
    const std::vector<Constraint<T>>& knownEquations, double timeLimit, IsZero isZero)
  {
    auto timeStarted = std::chrono::system_clock::now();

    std::size_t n = polyhedron->space()->dimension();
    auto redundancyCheck = EquationRedundancyCheck<T, IsZero>(n, isZero);
    for (auto equation : knownEquations)
    {
      if (redundancyCheck.add(equation) == EQUATION_INCONSISTENT)
      {
        // TODO: By inspecting the multipliers one could in principle find a smaller inconsistent set.

        outerEquations.clear();
        for (std::size_t e = 0; e < redundancyCheck.rank(); ++e)
          outerEquations.push_back(redundancyCheck.getEquation(e));
        outerEquations.push_back(equation);
        return -1;
      }
    }

    std::cout << "Initial " << knownEquations.size() << " equations have rank "
      << redundancyCheck.rank() << "." << std::endl;

    auto affineComplement = AffineComplement<T, IsZero>(n, isZero);
    int lowerBound = -1;
    int upperBound = n - redundancyCheck.rank();
    T* objective = new T[n];
    while (lowerBound < upperBound)
    {
      std::cout << "Iteration of affine hull computation: " << lowerBound
        << " <= dim <= " << upperBound << "." << std::endl;

      if (remainingTime(timeStarted, timeLimit) <= 0)
        return AFFINEHULL_TIMEOUT;


      // If we have dim(P) rays but no points yet, we solve the feasibility problem.

      if (innerPoints.empty() && int(innerRays.size()) == upperBound)
      {
        typename OptimizationOracle<T>::Query query;
        for (std::size_t v = 0; v < n; ++v)
          objective[v] = 0;
        query.timeLimit = remainingTime(timeStarted, timeLimit);
        typename OptimizationOracle<T>::Result result = polyhedron->maximize(objective, query);

        if (result.isUnbounded())
          throw std::runtime_error("maximize(<zero vector>) returned unbounded result.");
        else if (result.isFeasible())
        {
          std::cout << "  -> adding last point" << std::endl;
          innerPoints.push_back(result.points.front().vector);
          ++lowerBound;
          return upperBound;
        }
        else
          throw std::runtime_error("maximize(<zero vector>) returned infeasible after some rays.");
      }

      sparse_vector<T> kernelDirectionVector;
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      while (kernelDirectionVector.empty())
      {
        kernelDirectionColumn = affineComplement.selectColumn();
        std::cout << "  Column " << kernelDirectionColumn << std::flush;

        affineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector);
        assert(!kernelDirectionVector.empty());

        std::cout << " has kernel vector " << kernelDirectionVector << std::flush;

        if (remainingTime(timeStarted, timeLimit) <= 0)
          return AFFINEHULL_TIMEOUT;

        if (redundancyCheck.test(kernelDirectionVector) == EQUATION_INDEPENDENT)
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

        if (remainingTime(timeStarted, timeLimit) <= 0)
          return AFFINEHULL_TIMEOUT;
      }

      // Prepare for maximize kernelDirectionVector.

      std::cout << "  Maximization" << std::flush;
      for (std::size_t v = 0; v < n; ++v)
        objective[v] = 0;
      for (const auto& iter : kernelDirectionVector)
        objective[iter.first] = iter.second;

      typename OptimizationOracle<T>::Query query;
      if (!innerPoints.empty())
      {
        query.minObjectiveValue = innerPoints.front() * kernelDirectionVector;
        std::cout << " with common objective value " << (double)query.minObjectiveValue;
      }
      query.timeLimit = remainingTime(timeStarted, timeLimit);
        if (query.timeLimit <= 0)
          return AFFINEHULL_TIMEOUT;
      typename OptimizationOracle<T>::Result result = polyhedron->maximize(objective, query);
      std::cout << " yields " << result.points.size() << " points and " << result.rays.size()
        << " rays." << std::endl;
      if (remainingTime(timeStarted, timeLimit) <= 0)
        return AFFINEHULL_TIMEOUT;

      if (result.isInfeasible())
      {
        outerEquations.clear();
        outerEquations.push_back(neverSatisfiedConstraint<T>());
        std::cout << "  -> infeasible." << std::endl;
      }
      else if (result.isUnbounded())
      {
        innerRays.push_back(result.rays.front().vector);
        affineComplement.add(result.rays.front().vector, 0, kernelDirectionColumn);
        std::cout << "  -> adding a ray." << std::endl;
        ++lowerBound;
      }
      else if (innerPoints.empty())
      {
        assert(result.isFeasible());

        // Add the maximizer.

        innerPoints.push_back(result.points.front().vector);
        affineComplement.add(innerPoints.back(), 1, n);
        ++lowerBound;

        // If another point has a different objective value, we also add that and continue.
        for (std::size_t p = 2; p < result.points.size(); ++p)
        {
          if (!isZero(result.points[p].objectiveValue - result.points.front().objectiveValue))
          {
            innerPoints.push_back(result.points[p].vector);
            affineComplement.add(innerPoints.back(), 1, kernelDirectionColumn);
            ++lowerBound;
            std::cout << "  -> adding two points with objective values "
              << (double)result.points.front().objectiveValue << " and " << (double)result.points[p].objectiveValue
              << "." << std::endl;
            break;
          }
        }
        if (innerPoints.size() == 2)
          continue;

        query.minObjectiveValue = result.points.front().objectiveValue;
        std::cout << "  -> adding a point with objective value "
          << (double)result.points.front().objectiveValue << std::endl;
      }
      else
      {
        assert(result.isFeasible());

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : result.points)
        {
          if (point.objectiveValue != query.minObjectiveValue)
          {
            innerPoints.push_back(point.vector);
            affineComplement.add(point.vector, 1, kernelDirectionColumn);
            ++lowerBound;
            std::cout << "  -> adding a point with objective value "
              << (double)result.points.front().objectiveValue << std::endl;
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
      query.minObjectiveValue *= -1;

      query.timeLimit = remainingTime(timeStarted, timeLimit);
      if (query.timeLimit <= 0)
        return AFFINEHULL_TIMEOUT;

      result = polyhedron->maximize(objective, query);
      std::cout << "  Minimization yields " << result.points.size() << " points and "
        << result.rays.size() << " rays." << std::endl;

      if (remainingTime(timeStarted, timeLimit) <= 0)
        return AFFINEHULL_TIMEOUT;

      if (result.isUnbounded())
      {
        innerRays.push_back(result.rays.front().vector);
        affineComplement.add(result.rays.front().vector, 0, kernelDirectionColumn);
        ++lowerBound;
        std::cout << "  -> adding a ray." << std::endl;
      }
      else
      {
        assert(result.isFeasible());

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : result.points)
        {
          if (!isZero(point.objectiveValue - query.minObjectiveValue))
          {
            innerPoints.push_back(point.vector);
            affineComplement.add(point.vector, 1, kernelDirectionColumn);
            ++lowerBound;
            std::cout << "  -> adding a point with objective value "
              << (double)result.points.front().objectiveValue << std::endl;
            success = true;
            break;
          }
        }
        if (success)
          continue;
      }

      // We have to add an equation.
      outerEquations.push_back(Constraint<T>(-query.minObjectiveValue, kernelDirectionVector,
        -query.minObjectiveValue));
      affineComplement.markEquation(kernelDirectionColumn);
      redundancyCheck.add(outerEquations.back());
      --upperBound;
      std::cout << "  -> adding an equation." << std::endl;
    }

    delete[] objective;

    return lowerBound;
  }

  int affineHull(std::shared_ptr<Polyhedron<double, DoubleIsZero>> polyhedron,
    std::vector<sparse_vector<double>>& innerPoints,
    std::vector<sparse_vector<double>>& innerRays,
    std::vector<Constraint<double>>& outerEquations,
    const std::vector<Constraint<double>>& knownEquations, double timeLimit)
  {
    return affineHullImplementation(polyhedron, innerPoints, innerRays, outerEquations,
      knownEquations, timeLimit, DoubleIsZero(1.0e-9));
  }

#if defined(IPO_WITH_GMP)

  int affineHull(std::shared_ptr<Polyhedron<rational, RationalIsZero>> polyhedron,
    std::vector<sparse_vector<rational>>& innerPoints,
    std::vector<sparse_vector<rational>>& innerRays,
    std::vector<Constraint<rational>>& outerEquations,
    const std::vector<Constraint<rational>>& knownEquations, double timeLimit)
  {
    return affineHullImplementation(polyhedron, innerPoints, innerRays, outerEquations,
      knownEquations, timeLimit, RationalIsZero());
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
