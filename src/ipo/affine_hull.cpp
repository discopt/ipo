#include <ipo/affine_hull.hpp>

#define IPO_DEBUG_AFFINE_HULL_CHECK // Uncomment to enable debug checks.
#define IPO_DEBUG_AFFINE_HULL_PRINT // Uncomment to print activity.

#include <ipo/arithmetic.hpp>

#include "lu.hpp"
#include "redundancy.hpp"

#include <iostream>

namespace ipo
{
  AffineHullQuery::AffineHullQuery(bool exact)
    : epsilonConstraints(exact ? 0.0 : 1.0e-9), epsilonSafety(exact ? 0.0 : 1.0e-6),
    epsilonFactorization(exact ? 0.0 : 1.0e-12), epsilonCoefficient(exact ? 0.0 : 1.0e-12),
    hybrid(true), timeLimit(std::numeric_limits<double>::infinity())
  {

  }
  
  template <typename T>
  class AffineComplement
  {
  public:
    AffineComplement(std::size_t numVariables)
      : _numVariables(numVariables), _lu(),
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

    
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
    void _debugAdd(const sparse_vector<T>& row, const T& last, std::size_t newBasicColumn,
      const std::vector<T>& kernelVector)
    {
      if (last == 0)
      {
        _debugDenseRays.push_back(std::vector<T>());
        for (std::size_t v = 0; v < numVariables(); ++v)
          _debugDenseRays.back().push_back(T(0));
        for (const auto& iter : row)
          _debugDenseRays.back()[iter.first] = iter.second;
      }
      else
      {
        _debugDensePoints.push_back(std::vector<T>());
        for (std::size_t v = 0; v < numVariables(); ++v)
          _debugDensePoints.back().push_back(T(0));
        for (const auto& iter : row)
          _debugDensePoints.back()[iter.first] = iter.second;
      }

      std::size_t countNonzeroProducts = 0;
      for (std::size_t p = 1; p < _debugDensePoints.size(); ++p)
      {
        T prod = 0;
        for (std::size_t v = 0; v < numVariables(); ++v)
          prod += (_debugDensePoints[p][v] - _debugDensePoints[0][v]) * kernelVector[v];
        if (fabs(convertNumber<double>(prod)) > 1.0e-12)
          ++countNonzeroProducts;
      }
      for (std::size_t r = 0; r < _debugDenseRays.size(); ++r)
      {
        T prod = 0;
        for (std::size_t v = 0; v < numVariables(); ++v)
          prod += _debugDenseRays[r][v] * kernelVector[v];
        if (fabs(convertNumber<double>(prod)) > 1.0e-12)
          ++countNonzeroProducts;
      }

      if (countNonzeroProducts > 1)
      {
        for (std::size_t p = 1; p < _debugDensePoints.size(); ++p)
        {
          T prod = 0;
          for (std::size_t v = 0; v < numVariables(); ++v)
            prod += (_debugDensePoints[p][v] - _debugDensePoints[0][v]) * kernelVector[v];
          if (fabs(convertNumber<double>(prod)) > 1.0e-12)
            std::cout << "Product for point " << p << "-0 is " << prod << std::endl;
        }
        for (std::size_t r = 0; r < _debugDenseRays.size(); ++r)
        {
          T prod = 0;
          for (std::size_t v = 0; v < numVariables(); ++v)
            prod += _debugDenseRays[r][v] * kernelVector[v];
          if (fabs(convertNumber<double>(prod)) > 1.0e-12)
            std::cout << "Product for ray " << r << " is " << prod << std::endl;
        }
      }
    }
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */

    void add(const sparse_vector<T>& row, const T& last, std::size_t newBasicColumn,
      double epsilonFactorization)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "\nAffineComplement::add(row=" << row << ", last=" << last << ", newBasicColumn="
        << newBasicColumn << ")." << std::endl;      
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

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
      _lu.extend(&newRow[0], &newColumn[0], newRow.back(), epsilonFactorization);

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

    void computeKernelVector(std::size_t column, sparse_vector<T>& vector,
      double epsilonCoefficient)
    {
      std::vector<T> rhs(rank(), T(0));
      const ColumnData& columnData = _columns[column];
      for (std::size_t i = 0; i < columnData.rows.size(); ++i)
        rhs[columnData.rows[i]] = columnData.entries[i];

      _lu.solveLeft(&rhs[0]);
      _lu.solveUpper(&rhs[0]);

      vector.clear();
      for (std::size_t c = 0; c < numVariables(); ++c)
      {
        if (c == column)
          vector.push_back(c, T(-1));
        else if (_columns[c].basisIndex < std::numeric_limits<std::size_t>::max()
          && fabs(convertNumber<double>(rhs[_columns[c].basisIndex])) > epsilonCoefficient)
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

#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
    std::vector<std::vector<T>> _debugDenseMatrix;
    std::vector<std::vector<T>> _debugDensePoints;
    std::vector<std::vector<T>> _debugDenseRays;
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */

    struct ColumnData
    {
      std::size_t basisIndex;
      std::vector<std::size_t> rows;
      std::vector<T> entries;
      bool definesEquation;
    };

    std::size_t _numVariables;
    IncrementalLUFactorization<T> _lu;
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

  template <typename T, typename U>
  static int filterGivenEquation(EquationRedundancyCheck<T>& redundancyCheck,
    const std::vector<Constraint<U>>& knownEquations,
    std::vector<Constraint<U>>& resultEquations, double epsilonConstraint,
    double epsilonFactorization, double epsilonCoefficient)
  {
    for (auto equation : knownEquations)
    {
      EquationRedundancy redundancy = redundancyCheck.add(convertConstraint<T>(equation),
        epsilonConstraint, epsilonFactorization, epsilonCoefficient);
      if (redundancy == EQUATION_INCONSISTENT)
      {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Found given equation inconsistency." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        resultEquations.clear();
        resultEquations.push_back(neverSatisfiedConstraint<U>());
        return -1;
      }
      else if (redundancy == EQUATION_INDEPENDENT)
        resultEquations.push_back(equation);
    }
    return redundancyCheck.numVariables() - redundancyCheck.rank();
  }

  template <typename T>
  void affineHullPure(std::shared_ptr<Polyhedron<T>> polyhedron,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultRays,
    std::vector<Constraint<T>>& resultEquations, const AffineHullQuery& query,
    const std::vector<Constraint<T>>& knownEquations, AffineHullResultCommon& result)
  {
    auto timeStarted = std::chrono::system_clock::now();
    auto timeComponent = std::chrono::system_clock::now();
    std::size_t n = polyhedron->space()->dimension();
    result.upperBound = n;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "affineHullPure() called for ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto redundancyCheck = EquationRedundancyCheck<T>(n);
    result.upperBound = filterGivenEquation(redundancyCheck, knownEquations, resultEquations,
      query.epsilonConstraints, query.epsilonFactorization, query.epsilonCoefficient);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "Initial " << knownEquations.size() << " equations have rank "
      << redundancyCheck.rank() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto affineComplement = AffineComplement<T>(n);
    std::vector<T> objective(n);
    while (result.lowerBound < result.upperBound)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "Iteration of affine hull computation: " << result.lowerBound
        << " <= dim <= " << result.upperBound << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        
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
        typename OptimizationOracle<T>::Response oracleResponse = polyhedron->maximize(
          &objective[0], oracleQuery);
        result.timeOracles += elapsedTime(timeComponent);

        if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
          throw std::runtime_error("maximize(<zero vector>) returned unbounded result.");
        else if (oracleResponse.outcome == OPTIMIZATION_FEASIBLE)
        {
          std::cout << "  -> adding last point" << std::endl;
          resultPoints.push_back(oracleResponse.points.front().vector);
          result.dimension = result.lowerBound = result.upperBound;
          return;
        }
        else
          throw std::runtime_error("maximize(<zero vector>) returned infeasible after some rays.");
      }

      sparse_vector<T> kernelDirectionVector;
      double kernelDirectionNorm;
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      while (kernelDirectionVector.empty())
      {
        timeComponent = std::chrono::system_clock::now();
        kernelDirectionColumn = affineComplement.selectColumn();
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "  Column " << kernelDirectionColumn << std::flush;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        affineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector,
          query.epsilonCoefficient);
        result.timeKernel += elapsedTime(timeComponent);
        assert(!kernelDirectionVector.empty());

        kernelDirectionNorm = euclideanNorm(kernelDirectionVector);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << " has kernel vector " << kernelDirectionVector << std::flush;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        
        if (elapsedTime(timeStarted) >= query.timeLimit)
        {
          result.dimension = AFFINEHULL_ERROR_TIMEOUT;
          return;
        }

        timeComponent = std::chrono::system_clock::now();
        EquationRedundancy redundancy = redundancyCheck.test(kernelDirectionVector,
          query.epsilonCoefficient);
        result.timeEquations += elapsedTime(timeComponent);
  
        if (redundancy == EQUATION_INDEPENDENT)
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << " which is independent of the equations." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          break;
        }
        else
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << " which is in the span of the equations." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
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
        std::cout << " with common objective value " << convertNumber<double>(oracleQuery.minObjectiveValue);
      }
      oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
      if (oracleQuery.timeLimit <= 0)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      timeComponent = std::chrono::system_clock::now();
      typename OptimizationOracle<T>::Response oracleResponse = polyhedron->maximize(
        &objective[0], oracleQuery);
      result.timeOracles += elapsedTime(timeComponent);

      std::cout << " yields " << oracleResponse << std::endl;
      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      if (oracleResponse.outcome == OPTIMIZATION_INFEASIBLE)
      {
        resultEquations.clear();
        resultEquations.push_back(neverSatisfiedConstraint<T>());
        result.dimension = result.upperBound = -1;
        std::cout << "  -> infeasible." << std::endl;
      }
      else if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
      {
        resultRays.push_back(oracleResponse.rays.front().vector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        affineComplement._debugAdd(*oracleResponse.rays.front().vector, 0, kernelDirectionColumn,
          objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        affineComplement.add(*oracleResponse.rays.front().vector, 0, kernelDirectionColumn,
          query.epsilonFactorization);
        result.timePointsRays += elapsedTime(timeComponent);
        std::cout << "  -> adding a ray." << std::endl;
        ++result.lowerBound;
      }
      else if (resultPoints.empty())
      {
        assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);
        if (oracleResponse.points.empty())
        {
          throw std::runtime_error(
            "Detected invalid oracle behavior: responded with feasible polyhedron but no points (no minObjectiveValue given).");
        }
        assert(!oracleResponse.points.empty());

        // Add the maximizer.

        resultPoints.push_back(oracleResponse.points.front().vector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        affineComplement._debugAdd(*resultPoints.back(), 1, n, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        affineComplement.add(*resultPoints.back(), 1, n, query.epsilonFactorization);
        result.timePointsRays += elapsedTime(timeComponent);
        ++result.lowerBound;

        // If another point has a different objective value, we also add that and continue.
        for (std::size_t p = 2; p < oracleResponse.points.size(); ++p)
        {
          double difference = fabs(convertNumber<double, mpq_class>(oracleResponse.points[p].objectiveValue
            - oracleResponse.points.front().objectiveValue));
          difference /= kernelDirectionNorm;
          if (difference > query.epsilonConstraints)
          {
            resultPoints.push_back(oracleResponse.points[p].vector);
            timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
            affineComplement._debugAdd(*resultPoints.back(), 1, kernelDirectionColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
            affineComplement.add(*resultPoints.back(), 1, kernelDirectionColumn,
              query.epsilonFactorization);
            result.timePointsRays += elapsedTime(timeComponent);
            ++result.lowerBound;
            std::cout << "  -> adding two points with objective values "
              << convertNumber<double>(oracleResponse.points.front().objectiveValue) << " and " 
              << convertNumber<double>(oracleResponse.points[p].objectiveValue) << "." << std::endl;
            break;
          }
        }
        if (resultPoints.size() == 2)
          continue;

        oracleQuery.minObjectiveValue = oracleResponse.points.front().objectiveValue;
        std::cout << "  -> adding a point with objective value "
          << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
      }
      else
      {
        assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : oracleResponse.points)
        {
          double difference = fabs(convertNumber<double, mpq_class>(point.objectiveValue - oracleQuery.minObjectiveValue));
          difference /= kernelDirectionNorm;
          if (difference > query.epsilonConstraints)
          {
            resultPoints.push_back(point.vector);
            timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
            affineComplement._debugAdd(*point.vector, 1, kernelDirectionColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
            affineComplement.add(*point.vector, 1, kernelDirectionColumn,
              query.epsilonFactorization);
            result.timePointsRays += elapsedTime(timeComponent);
            ++result.lowerBound;
            std::cout << "  -> adding a point with objective value "
              << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
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
      oracleResponse = polyhedron->maximize(&objective[0], oracleQuery);
      result.timeOracles += elapsedTime(timeComponent);

      std::cout << "  Minimization yields " << oracleResponse << "." << std::endl;

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
      {
        resultRays.push_back(oracleResponse.rays.front().vector);

        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        affineComplement._debugAdd(*oracleResponse.rays.front().vector, 0, kernelDirectionColumn,
          objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        affineComplement.add(*oracleResponse.rays.front().vector, 0, kernelDirectionColumn,
          query.epsilonFactorization);
        result.timePointsRays += elapsedTime(timeComponent);

        ++result.lowerBound;
        std::cout << "  -> adding a ray." << std::endl;
      }
      else
      {
        assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : oracleResponse.points)
        {
          double difference = fabs(convertNumber<double, mpq_class>(point.objectiveValue - oracleQuery.minObjectiveValue));
          difference /= kernelDirectionNorm;
          if (difference > query.epsilonConstraints)
          {
            resultPoints.push_back(point.vector);
            timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
            affineComplement._debugAdd(*point.vector, 1, kernelDirectionColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
            affineComplement.add(*point.vector, 1, kernelDirectionColumn,
              query.epsilonFactorization);
            result.timePointsRays += elapsedTime(timeComponent);
            ++result.lowerBound;
            std::cout << "  -> adding a point with objective value "
              << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
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
      redundancyCheck.add(resultEquations.back(), query.epsilonConstraints,
        query.epsilonFactorization, query.epsilonCoefficient);
      result.timeEquations += elapsedTime(timeComponent);

      --result.upperBound;
      std::cout << "  -> adding an equation." << std::endl;
    }

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
  }

  void affineHullHybrid(std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    std::vector<std::shared_ptr<sparse_vector<mpq_class>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<mpq_class>>>& resultRays,
    std::vector<Constraint<mpq_class>>& resultEquations, const AffineHullQuery& query,
    const std::vector<Constraint<mpq_class>>& knownEquations, AffineHullResultCommon& result)
  {
    auto timeStarted = std::chrono::system_clock::now();
    auto timeComponent = std::chrono::system_clock::now();
    std::size_t n = polyhedron->space()->dimension();
    result.upperBound = n;
    std::vector<std::shared_ptr<sparse_vector<mpq_class>>> candidatePoints;
    std::vector<std::shared_ptr<sparse_vector<mpq_class>>> candidateRays;
    std::vector<Constraint<mpq_class>> candidateEquations;
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "affineHullHybrid() called with ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto approximateRedundancyCheck = EquationRedundancyCheck<double>(n);
    result.upperBound = filterGivenEquation(approximateRedundancyCheck, knownEquations,
      candidateEquations, query.epsilonConstraints, query.epsilonFactorization,
      query.epsilonCoefficient);
    if (result.upperBound == -1)
    {
      AffineHullQuery exactQuery = query;
      exactQuery.timeLimit -= elapsedTime(timeStarted);
      exactQuery.epsilonCoefficient = 0.0;
      exactQuery.epsilonFactorization = 0.0;
      exactQuery.epsilonConstraints = 0.0;
      exactQuery.epsilonSafety = 0.0;
      affineHullPure(polyhedron, resultPoints, resultRays, resultEquations, query, knownEquations,
        result);
      return;
    }

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "Initial " << knownEquations.size() << " equations have approximate rank "
      << approximateRedundancyCheck.rank() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
    int candidateLowerBound = -1;
    int candidateUpperBound = n - approximateRedundancyCheck.rank();

    auto approximateAffineComplement = AffineComplement<double>(n);
    std::vector<double> approximateObjective(n);
    while (candidateLowerBound < candidateUpperBound)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "Iteration of approximate affine hull computation: " << result.lowerBound
        << " <= dim <= " << result.upperBound << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      // If we have dim(P) rays but no points yet, we solve the feasibility problem.

      if (candidateEquations.empty() && int(candidateRays.size()) == candidateUpperBound)
      {
        typename OptimizationOracle<mpq_class>::Query oracleQuery;
        for (std::size_t v = 0; v < n; ++v)
          approximateObjective[v] = 0;
        oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
        timeComponent = std::chrono::system_clock::now();
        typename OptimizationOracle<mpq_class>::Response oracleResponse
          = polyhedron->maximizeDouble(&approximateObjective[0], oracleQuery);
        result.timeOracles += elapsedTime(timeComponent);

        if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
          throw std::runtime_error("maximize(<zero vector>) returned unbounded result.");
        else if (oracleResponse.outcome == OPTIMIZATION_FEASIBLE)
        {
          std::cout << "  -> adding last point" << std::endl;
          candidatePoints.push_back(oracleResponse.points.front().vector);
          candidateLowerBound = candidateUpperBound;
          break;
        }
        else
          throw std::runtime_error("maximize(<zero vector>) returned infeasible after some rays.");
      }

      sparse_vector<double> kernelDirectionVector;
      double kernelDirectionNorm;
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      while (kernelDirectionVector.empty())
      {
        timeComponent = std::chrono::system_clock::now();
        kernelDirectionColumn = approximateAffineComplement.selectColumn();
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "  Column " << kernelDirectionColumn << std::flush;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        approximateAffineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector,
          query.epsilonCoefficient);
        result.timeKernel += elapsedTime(timeComponent);
        assert(!kernelDirectionVector.empty());

        kernelDirectionNorm = euclideanNorm(kernelDirectionVector);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << " has kernel vector " << kernelDirectionVector << std::flush;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        
        if (elapsedTime(timeStarted) >= query.timeLimit)
        {
          result.dimension = AFFINEHULL_ERROR_TIMEOUT;
          return;
        }

        timeComponent = std::chrono::system_clock::now();
        EquationRedundancy redundancy = approximateRedundancyCheck.test(kernelDirectionVector,
          query.epsilonCoefficient);
        result.timeEquations += elapsedTime(timeComponent);
  
        if (redundancy == EQUATION_INDEPENDENT)
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << " which is independent of the equations." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          break;
        }
        else
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << " which is in the span of the equations." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          approximateAffineComplement.markEquation(kernelDirectionColumn);
          kernelDirectionVector.clear();
        }

        if (elapsedTime(timeStarted) >= query.timeLimit)
        {
          result.dimension = AFFINEHULL_ERROR_TIMEOUT;
          return;
        }
      }

      // Prepare for maximizing kernelDirectionVector.

      std::cout << "  Maximization" << std::flush;
      for (std::size_t v = 0; v < n; ++v)
        approximateObjective[v] = 0;
      for (const auto& iter : kernelDirectionVector)
        approximateObjective[iter.first] = iter.second;

      typename OptimizationOracle<mpq_class>::Query oracleQuery;
      if (!candidatePoints.empty())
      {
        oracleQuery.minObjectiveValue = *candidatePoints.front() * kernelDirectionVector;
        std::cout << " with common objective value " << convertNumber<double>(
          oracleQuery.minObjectiveValue);
      }
      oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
      if (oracleQuery.timeLimit <= 0)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      timeComponent = std::chrono::system_clock::now();
      typename OptimizationOracle<mpq_class>::Response oracleResponse = polyhedron->maximizeDouble(
        &approximateObjective[0], oracleQuery);
      result.timeOracles += elapsedTime(timeComponent);

      std::cout << " yields " << oracleResponse << std::endl;
      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      if (oracleResponse.outcome == OPTIMIZATION_INFEASIBLE)
      {
        assert(resultEquations.empty());
        resultEquations.push_back(neverSatisfiedConstraint<mpq_class>());
        result.dimension = result.upperBound = -1;
        std::cout << "  -> infeasible." << std::endl;
      }
      else if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
      {
        candidateRays.push_back(oracleResponse.rays.front().vector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        approximateAffineComplement._debugAdd(convertTo<double>(*oracleResponse.rays.front().vector),
          0, kernelDirectionColumn, approximateObjective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        approximateAffineComplement.add(convertTo<double>(*oracleResponse.rays.front().vector), 0,
          kernelDirectionColumn, query.epsilonFactorization);
        result.timePointsRays += elapsedTime(timeComponent);
        std::cout << "  -> adding a ray." << std::endl;
        ++candidateLowerBound;
      }
      else if (resultPoints.empty())
      {
        assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);
        if (oracleResponse.points.empty())
        {
          throw std::runtime_error(
            "Detected invalid oracle behavior: responded with feasible polyhedron but no points (no minObjectiveValue given).");
        }
        assert(!oracleResponse.points.empty());

        // Add the maximizer.

        resultPoints.push_back(oracleResponse.points.front().vector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        approximateAffineComplement._debugAdd(convertTo<double>(*resultPoints.back()), 1, n,
          approximateObjective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        approximateAffineComplement.add(convertTo<double>(*resultPoints.back()), 1, n,
          query.epsilonFactorization);
        result.timePointsRays += elapsedTime(timeComponent);
        ++candidateLowerBound;

        // If another point has a different objective value, we also add that and continue.
        for (std::size_t p = 2; p < oracleResponse.points.size(); ++p)
        {
          double difference = fabs(convertNumber<double, mpq_class>(oracleResponse.points[p].objectiveValue
            - oracleResponse.points.front().objectiveValue));
          difference /= kernelDirectionNorm;
          if (difference > query.epsilonConstraints)
          {
            candidatePoints.push_back(oracleResponse.points[p].vector);
            timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
            approximateAffineComplement._debugAdd(convertTo<double>(*resultPoints.back()), 1,
              kernelDirectionColumn, approximateObjective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
            approximateAffineComplement.add(convertTo<double>(*resultPoints.back()), 1,
              kernelDirectionColumn, query.epsilonFactorization);
            result.timePointsRays += elapsedTime(timeComponent);
            ++candidateLowerBound;
            std::cout << "  -> adding two points with objective values "
              << convertNumber<double>(oracleResponse.points.front().objectiveValue) << " and " 
              << convertNumber<double>(oracleResponse.points[p].objectiveValue) << "." << std::endl;
            break;
          }
        }
        if (candidatePoints.size() == 2)
          continue;

        oracleQuery.minObjectiveValue = oracleResponse.points.front().objectiveValue;
        std::cout << "  -> adding a point with objective value "
          << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
      }
      else
      {
        assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : oracleResponse.points)
        {
          double difference = fabs(convertNumber<double, mpq_class>(point.objectiveValue - oracleQuery.minObjectiveValue));
          difference /= kernelDirectionNorm;
          if (difference > query.epsilonConstraints)
          {
            candidatePoints.push_back(point.vector);
            timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
            approximateAffineComplement._debugAdd(convertTo<double>(*point.vector), 1,
              kernelDirectionColumn, approximateObjective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
            approximateAffineComplement.add(convertTo<double>(*point.vector), 1, kernelDirectionColumn,
              query.epsilonFactorization);
            result.timePointsRays += elapsedTime(timeComponent);
            ++candidateLowerBound;
            std::cout << "  -> adding a point with objective value "
              << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
            success = true;
            break;
          }
        }
        if (success)
          continue;
      }

      // Prepare minimization.

      for (std::size_t v = 0; v < n; ++v)
        approximateObjective[v] *= -1;
      oracleQuery.minObjectiveValue *= -1;

      oracleQuery.timeLimit = query.timeLimit - elapsedTime(timeStarted);
      if (oracleQuery.timeLimit <= 0)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      timeComponent = std::chrono::system_clock::now();
      oracleResponse = polyhedron->maximizeDouble(&approximateObjective[0], oracleQuery);
      result.timeOracles += elapsedTime(timeComponent);

      std::cout << "  Minimization yields " << oracleResponse << "." << std::endl;

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
      {
        candidateRays.push_back(oracleResponse.rays.front().vector);

        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        approximateAffineComplement._debugAdd(convertTo<double>(*oracleResponse.rays.front().vector),
          0, kernelDirectionColumn, approximateObjective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        approximateAffineComplement.add(convertTo<double>(*oracleResponse.rays.front().vector), 0,
          kernelDirectionColumn, query.epsilonFactorization);
        result.timePointsRays += elapsedTime(timeComponent);

        ++candidateLowerBound;
        std::cout << "  -> adding a ray." << std::endl;
      }
      else
      {
        assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

        // Check if some vector has a different objective value.
        bool success = false;
        for (const auto& point : oracleResponse.points)
        {
          double difference = fabs(convertNumber<double, mpq_class>(point.objectiveValue - oracleQuery.minObjectiveValue));
          difference /= kernelDirectionNorm;
          if (difference > query.epsilonConstraints)
          {
            candidatePoints.push_back(point.vector);
            timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
            approximateAffineComplement._debugAdd(convertTo<double>(*point.vector), 1,
              kernelDirectionColumn, approximateObjective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
            approximateAffineComplement.add(convertTo<double>(*point.vector), 1, kernelDirectionColumn,
              query.epsilonFactorization);
            result.timePointsRays += elapsedTime(timeComponent);
            ++candidateLowerBound;
            std::cout << "  -> adding a point with objective value "
              << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
            success = true;
            break;
          }
        }
        if (success)
          continue;
      }

      auto vector = std::make_shared<sparse_vector<double>>(std::move(kernelDirectionVector));

      // We create the approximate equation.
      auto approximateEquation = Constraint<double>(-convertNumber<double>(oracleQuery.minObjectiveValue),
        vector, -convertNumber<double>(oracleQuery.minObjectiveValue));
      approximateAffineComplement.markEquation(kernelDirectionColumn);

      timeComponent = std::chrono::system_clock::now();
      approximateRedundancyCheck.add(approximateEquation, query.epsilonConstraints,
        query.epsilonFactorization, query.epsilonCoefficient);
      result.timeEquations += elapsedTime(timeComponent);

      --candidateUpperBound;
      std::cout << "  -> adding an equation." << std::endl;
    }

    assert(false);

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
  }

  AffineHullResult<double> affineHull(
    std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<double>>& knownEquations)
  {
    AffineHullResult<double> result;
    affineHullPure(polyhedron, result.points, result.rays, result.equations, query,
      knownEquations, result);
    return result;
  }

#if defined(IPO_WITH_GMP)

  AffineHullResult<mpq_class> affineHull(
    std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<mpq_class>>& knownEquations)
  {
    if (query.epsilonConstraints == 0.0 || query.epsilonCoefficient == 0.0
      || query.epsilonFactorization == 0.0 || query.epsilonSafety == 0.0)
    {
      AffineHullResult<mpq_class> result;
      affineHullPure(polyhedron, result.points, result.rays, result.equations, query,
        knownEquations, result);
      return result;
    }
    else
    {
      AffineHullResult<mpq_class> result;
      affineHullHybrid(polyhedron, result.points, result.rays, result.equations, query,
        knownEquations, result);
      return result;
    }
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
