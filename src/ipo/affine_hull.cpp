#include <ipo/affine_hull.hpp>

// #define IPO_DEBUG_AFFINE_HULL_CHECK // Uncomment to enable debug checks.
// #define IPO_DEBUG_AFFINE_HULL_PRINT // Uncomment to print activity.

#include <ipo/arithmetic.hpp>

#include "lu.hpp"
#include "redundancy.hpp"

#include <iostream>

namespace ipo
{
  AffineHullQuery::AffineHullQuery()
    : epsilonConstraints(1.0e-9), epsilonSafety(1.0e-6), epsilonFactorization(1.0e-12),
    epsilonCoefficient(1.0e-12), algorithm(AFFINEHULL_ALGORITHM_DIRECT),
    timeLimit(std::numeric_limits<double>::infinity())
  {

  }

  AffineHullQuery& AffineHullQuery::operator=(const AffineHullQuery& other)
  {
    epsilonConstraints = other.epsilonConstraints;
    epsilonSafety = other.epsilonSafety;
    epsilonCoefficient = other.epsilonCoefficient;
    epsilonFactorization = other.epsilonFactorization;
    algorithm = other.algorithm;
    timeLimit = other.timeLimit;
    return *this;
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
      const T* kernelVector)
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
  static int filterGivenEquation(EquationRedundancyCheck<U>& redundancyCheck,
    const std::vector<Constraint<T>>& knownEquations,
    std::vector<Constraint<U>>* pResultEquations, double epsilonConstraint,
    double epsilonFactorization, double epsilonCoefficient,
    EquationRedundancyCheck<double>* pApproximateRedundancyCheck = NULL)
  {
    for (auto equation : knownEquations)
    {
      auto convertedEquation = convertConstraint<U>(equation);
      EquationRedundancy redundancy = redundancyCheck.add(convertedEquation,
        epsilonConstraint, epsilonFactorization, epsilonCoefficient);
      if (pApproximateRedundancyCheck)
      {
        pApproximateRedundancyCheck->add(convertConstraint<double>(equation),
          epsilonConstraint, epsilonFactorization, epsilonCoefficient);
      }
      if (redundancy == EQUATION_INCONSISTENT)
      {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Found given equation inconsistency." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        if (pResultEquations)
        {
          pResultEquations->clear();
          pResultEquations->push_back(neverSatisfiedConstraint<U>());
        }
        return -1;
      }
      else if (redundancy == EQUATION_INDEPENDENT)
      {
        if (pResultEquations)
          pResultEquations->push_back(convertedEquation);
      }
    }
    return redundancyCheck.numVariables() - redundancyCheck.rank();
  }

#if defined(IPO_WITH_GMP)

  static OptimizationOracle<mpq_class>::Response maximize(
    std::shared_ptr<Polyhedron<mpq_class>> polyhedron, const mpq_class* objective,
    const typename OptimizationOracle<mpq_class>::Query& query)
  {
    return polyhedron->maximize(objective, query);
  }
#endif /* IPO_WITH_GMP */

  static OptimizationOracle<double>::Response maximize(
    std::shared_ptr<Polyhedron<double>> polyhedron, const double* objective,
    const typename OptimizationOracle<double>::Query& query)
  {
    return polyhedron->maximize(objective, query);
  }

#if defined(IPO_WITH_GMP)

  static OptimizationOracle<mpq_class>::Response maximize(
    std::shared_ptr<Polyhedron<mpq_class>> polyhedron, const double* objective,
    const typename OptimizationOracle<mpq_class>::Query& query)
  {
    return polyhedron->maximizeDouble(objective, query);
  }

#endif /* IPO_WITH_GMP */

  template <typename T, typename U>
  static void findLastPoint(std::shared_ptr<Polyhedron<T>> polyhedron, U* objective,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints, double timeLimit,
    double& timeOracles)
  {
    typename OptimizationOracle<T>::Query oracleQuery;
    for (std::size_t v = 0; v < polyhedron->space()->dimension(); ++v)
      objective[v] = 0;
    oracleQuery.timeLimit = timeLimit;
    auto timeComponent = std::chrono::system_clock::now();
    auto oracleResponse = maximize(polyhedron, objective, oracleQuery);
    timeOracles += elapsedTime(timeComponent);

    if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
      throw std::runtime_error("maximize(<zero vector>) returned unbounded result.");
    else if (oracleResponse.outcome == OPTIMIZATION_FEASIBLE)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "  -> adding last point" << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
      resultPoints.push_back(oracleResponse.points.front().vector);
      return;
    }
    else
      throw std::runtime_error("maximize(<zero vector>) returned infeasible after some rays.");
  }

  static const int INTERNAL_TIMEOUT = AFFINEHULL_ERROR_TIMEOUT;
  static const int INTERNAL_INFEASIBLE = -1;
  static const int INTERNAL_NONE = 0;
  static const int INTERNAL_INITIAL_ONE = 1;
  static const int INTERNAL_INITIAL_TWO = 2;
  static const int INTERNAL_FURTHER_ONE = 3;

  template <typename T, typename U>
  static int tryMaximization(std::shared_ptr<Polyhedron<T>> polyhedron,
    AffineComplement<U>& affineComplement, U* objective,
    const sparse_vector<U>& kernelDirectionVector, std::size_t kernelDirectionColumn,
    double kernelDirectionNorm, U& kernelDirectionValue,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultRays, double epsilonFactorization,
    double epsilonConstraints, double timeLimit, double& timeOracles, double& timePointsRays,
    AffineComplement<double>* pAffineComplement = NULL)
  {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "  Maximization" << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    for (std::size_t v = 0; v < polyhedron->space()->dimension(); ++v)
      objective[v] = 0;
    for (const auto& iter : kernelDirectionVector)
      objective[iter.first] = iter.second;

    typename OptimizationOracle<T>::Query oracleQuery;
    if (!resultPoints.empty())
    {
      oracleQuery.hasMinObjectiveValue = true;
      oracleQuery.minObjectiveValue = kernelDirectionValue;
    }
    if (timeLimit <= 0)
      return INTERNAL_TIMEOUT;
    oracleQuery.timeLimit = timeLimit;

    auto timeComponent = std::chrono::system_clock::now();
    auto oracleResponse = maximize(polyhedron, objective, oracleQuery);
    double timeOracle = elapsedTime(timeComponent);
    timeOracles += timeOracle;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "Result: " << oracleResponse << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    if (timeOracle >= timeLimit)
      return INTERNAL_TIMEOUT;

    if (oracleResponse.outcome == OPTIMIZATION_INFEASIBLE)
      return INTERNAL_INFEASIBLE;
    else if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
    {
      resultRays.push_back(oracleResponse.rays.front().vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertTo<U>(*oracleResponse.rays.front().vector), 0,
        kernelDirectionColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      affineComplement.add(convertTo<U>(*oracleResponse.rays.front().vector), 0,
        kernelDirectionColumn, epsilonFactorization);
      if (pAffineComplement)
      {
        pAffineComplement->add(convertTo<double>(*oracleResponse.rays.front().vector), 0,
          kernelDirectionColumn, epsilonFactorization);
      }
      timePointsRays += elapsedTime(timeComponent);
      if (resultPoints.empty() && !oracleResponse.points.empty())
      {
        resultPoints.push_back(oracleResponse.points.front().vector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        affineComplement._debugAdd(convertTo<U>(*oracleResponse.points.front().vector), 1,
          polyhedron->space()->dimension(), objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        affineComplement.add(convertTo<U>(*oracleResponse.points.front().vector), 1,
          polyhedron->space()->dimension(), epsilonFactorization);
        if (pAffineComplement)
        {
          pAffineComplement->add(convertTo<double>(*oracleResponse.points.front().vector), 1,
            polyhedron->space()->dimension(), epsilonFactorization);
        }
        timePointsRays += elapsedTime(timeComponent);
        return INTERNAL_INITIAL_TWO;
      }
      return resultPoints.empty() ? INTERNAL_INITIAL_ONE : INTERNAL_FURTHER_ONE;
    }
    else if (resultPoints.empty())
    {
      assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

      if (oracleResponse.points.empty())
      {
        throw std::runtime_error(
          "Detected invalid oracle behavior: responded with feasible polyhedron but no points.");
      }
      assert(!oracleResponse.points.empty());

      // Add the maximizer.

      resultPoints.push_back(oracleResponse.points.front().vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertTo<U>(*resultPoints.back()), 1,
        polyhedron->space()->dimension(), objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      affineComplement.add(convertTo<U>(*resultPoints.back()), 1, polyhedron->space()->dimension(),
        epsilonFactorization);
      if (pAffineComplement)
      {
        pAffineComplement->add(convertTo<double>(*resultPoints.back()), 1,
          polyhedron->space()->dimension(), epsilonFactorization);
      }
      timePointsRays += elapsedTime(timeComponent);

      // If another point has a different objective value, we also add that and continue.
      for (std::size_t p = 2; p < oracleResponse.points.size(); ++p)
      {
        double difference = fabs(convertNumber<double, T>(oracleResponse.points[p].objectiveValue
          - oracleResponse.points.front().objectiveValue));
        difference /= kernelDirectionNorm;
        if (difference > epsilonConstraints)
        {
          resultPoints.push_back(oracleResponse.points[p].vector);
          timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
          affineComplement._debugAdd(convertTo<U>(*resultPoints.back()), 1, kernelDirectionColumn,
            objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
          affineComplement.add(convertTo<U>(*resultPoints.back()), 1, kernelDirectionColumn,
            epsilonFactorization);
          if (pAffineComplement)
          {
            pAffineComplement->add(convertTo<double>(*resultPoints.back()), 1,
              kernelDirectionColumn, epsilonFactorization);
          }
          timePointsRays += elapsedTime(timeComponent);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "  -> adding two points with objective values "
            << convertNumber<double>(oracleResponse.points.front().objectiveValue) << " and "
            << convertNumber<double>(oracleResponse.points[p].objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          return INTERNAL_INITIAL_TWO;
        }
      }
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "  -> adding a point with objective value "
        << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
      kernelDirectionValue = convertNumber<U>(oracleResponse.points.front().objectiveValue);
      return INTERNAL_INITIAL_ONE;
    }
    else
    {
      assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

      // Check if some vector has a different objective value.
      for (const auto& point : oracleResponse.points)
      {
        double difference = fabs(convertNumber<double, T>(
          point.objectiveValue - oracleQuery.minObjectiveValue));
        difference /= kernelDirectionNorm;
        if (difference > epsilonConstraints)
        {
          resultPoints.push_back(point.vector);
          timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
          affineComplement._debugAdd(convertTo<U>(*point.vector), 1, kernelDirectionColumn,
            objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
          affineComplement.add(convertTo<U>(*point.vector), 1, kernelDirectionColumn,
            epsilonFactorization);
          if (pAffineComplement)
          {
            pAffineComplement->add(convertTo<double>(*point.vector), 1,
              kernelDirectionColumn, epsilonFactorization);
          }
          timePointsRays += elapsedTime(timeComponent);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "  -> adding a point with objective value "
            << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          return INTERNAL_FURTHER_ONE;
        }
      }
      return INTERNAL_NONE;
    }
    return std::numeric_limits<int>::max();
  }

  
  template <typename T, typename U>
  static int tryMinimization(std::shared_ptr<Polyhedron<T>> polyhedron,
    AffineComplement<U>& affineComplement, U* objective,
    const sparse_vector<U>& kernelDirectionVector, std::size_t kernelDirectionColumn,
    double kernelDirectionNorm, const U& kernelDirectionValue,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultRays, double epsilonFactorization,
    double epsilonConstraints, double timeLimit, double& timeOracles, double& timePointsRays,
    AffineComplement<double>* pAffineComplement = NULL)
  {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "  Minimization" << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    for (std::size_t v = 0; v < polyhedron->space()->dimension(); ++v)
      objective[v] = -objective[v];

#if !defined(NDEBUG)
    for (const auto& iter : kernelDirectionVector)
      assert(objective[iter.first] == -iter.second);
#endif /* !NDEBUG */

    typename OptimizationOracle<T>::Query oracleQuery;
    if (!resultPoints.empty())
    {
      oracleQuery.hasMinObjectiveValue = true;
      oracleQuery.minObjectiveValue = -kernelDirectionValue;
    }
    if (timeLimit <= 0)
      return INTERNAL_TIMEOUT;
    oracleQuery.timeLimit = timeLimit;

    auto timeComponent = std::chrono::system_clock::now();
    auto oracleResponse = maximize(polyhedron, objective, oracleQuery);
    double timeOracle = elapsedTime(timeComponent);
    timeOracles += timeOracle;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "Result: " << oracleResponse << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    if (timeOracle >= timeLimit)
      return INTERNAL_TIMEOUT;

    if (oracleResponse.outcome == OPTIMIZATION_INFEASIBLE)
      throw std::runtime_error("Oracle for minimization claims infeasible.");
    else if (oracleResponse.outcome == OPTIMIZATION_UNBOUNDED)
    {
      resultRays.push_back(oracleResponse.rays.front().vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertTo<U>(*oracleResponse.rays.front().vector), 0,
        kernelDirectionColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      affineComplement.add(convertTo<U>(*oracleResponse.rays.front().vector), 0,
        kernelDirectionColumn, epsilonFactorization);
      if (pAffineComplement)
      {
        pAffineComplement->add(convertTo<double>(*oracleResponse.rays.front().vector), 0,
          kernelDirectionColumn, epsilonFactorization);
      }
      timePointsRays += elapsedTime(timeComponent);
      return INTERNAL_FURTHER_ONE;
    }
    else
    {
      assert(!resultPoints.empty());
      assert(oracleResponse.outcome == OPTIMIZATION_FEASIBLE);

      // Check if some vector has a different objective value.
      if (oracleResponse.points.empty())
        return INTERNAL_NONE;

      double difference = fabs(convertNumber<double, T>(
        oracleResponse.points.front().objectiveValue - oracleQuery.minObjectiveValue));
      difference /= kernelDirectionNorm;
      if (difference <= epsilonConstraints)
        return INTERNAL_NONE;
  
      resultPoints.push_back(oracleResponse.points.front().vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertTo<U>(*oracleResponse.points.front().vector), 1,
        kernelDirectionColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      affineComplement.add(convertTo<U>(*oracleResponse.points.front().vector), 1,
        kernelDirectionColumn, epsilonFactorization);
      if (pAffineComplement)
      {
        pAffineComplement->add(convertTo<double>(*oracleResponse.points.front().vector), 1,
          kernelDirectionColumn, epsilonFactorization);
      }
      timePointsRays += elapsedTime(timeComponent);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "  -> adding a point with objective value "
          << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
      return INTERNAL_FURTHER_ONE;
    }
    return std::numeric_limits<int>::max();
  }

  template <typename T>
  static void addEquation(EquationRedundancyCheck<T>& redundancyCheck,
    AffineComplement<T>& affineComplement,
    std::vector<Constraint<T>>& resultEquations,
    const sparse_vector<T>& kernelDirectionVector, std::size_t kernelDirectionColumn,
    const T& kernelDirectionValue, double epsilonFactorization, double epsilonConstraints,
    double epsilonCoefficient, double& timeEquations,
    EquationRedundancyCheck<double>* pRedundancyCheck = NULL,
    AffineComplement<double>* pAffineComplement = NULL)
  {
    auto vector = std::make_shared<sparse_vector<T>>(std::move(kernelDirectionVector));
    resultEquations.push_back(Constraint<T>(kernelDirectionValue, vector, kernelDirectionValue));
    affineComplement.markEquation(kernelDirectionColumn);
    if (pAffineComplement)
      pAffineComplement->markEquation(kernelDirectionColumn);

    auto timeComponent = std::chrono::system_clock::now();
    redundancyCheck.add(resultEquations.back(), epsilonConstraints, epsilonFactorization,
      epsilonCoefficient);
    if (pRedundancyCheck)
    {
      pRedundancyCheck->add(convertConstraint<double>(resultEquations.back()), epsilonConstraints,
        epsilonCoefficient, epsilonCoefficient);
    }
    timeEquations += elapsedTime(timeComponent);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "  -> adding an equation." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
  }

  template <typename T>
  static void affineHullDirect(std::shared_ptr<Polyhedron<T>> polyhedron,
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
    std::cout << "affineHullDirect() called for ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto redundancyCheck = EquationRedundancyCheck<T>(n);
    result.upperBound = filterGivenEquation(redundancyCheck, knownEquations, &resultEquations,
      query.epsilonConstraints, query.epsilonFactorization, query.epsilonCoefficient);
    if (result.upperBound == -1)
    {
      result.dimension = -1;
      return;
    }

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
        findLastPoint(polyhedron, &objective[0], resultPoints,
          query.timeLimit - elapsedTime(timeStarted), result.timeOracles);
        result.dimension = result.lowerBound = result.upperBound;
        return;
      }

      sparse_vector<T> kernelDirectionVector;
      double kernelDirectionNorm = std::numeric_limits<double>::signaling_NaN();
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      T kernelDirectionValue = 0;
      while (kernelDirectionVector.empty())
      {
        timeComponent = std::chrono::system_clock::now();
        kernelDirectionColumn = affineComplement.selectColumn();
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "  Column " << kernelDirectionColumn << std::flush;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        affineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector,
          query.epsilonCoefficient);
        result.numKernel++;
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

      // Maximize kernelDirectionVector.

      int outcome = tryMaximization(polyhedron, affineComplement, &objective[0],
        kernelDirectionVector, kernelDirectionColumn, kernelDirectionNorm, kernelDirectionValue,
        resultPoints, resultRays, query.epsilonFactorization, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays);
      if (outcome == INTERNAL_TIMEOUT)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_INFEASIBLE)
      {
        result.dimension = -1;
        return;
      }
      else if (outcome == INTERNAL_INITIAL_ONE)
      {
        result.lowerBound++;
      }
      else if (outcome == INTERNAL_INITIAL_TWO)
      {
        result.lowerBound += 2;
        continue;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        result.lowerBound++;
        continue;
      }

      // Minimize kernelDirectionVector.

      outcome = tryMinimization(polyhedron, affineComplement, &objective[0],
        kernelDirectionVector, kernelDirectionColumn, kernelDirectionNorm, kernelDirectionValue,
        resultPoints, resultRays, query.epsilonFactorization, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays);
      if (outcome == INTERNAL_TIMEOUT)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        result.lowerBound++;
        continue;
      }

      // Add equation defined by kernelDirectionVector.

      addEquation(redundancyCheck, affineComplement, resultEquations, kernelDirectionVector,
        kernelDirectionColumn, kernelDirectionValue, query.epsilonFactorization,
        query.epsilonConstraints, query.epsilonCoefficient, result.timeEquations);
      --result.upperBound;
    }

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
  }

#if defined(IPO_WITH_GMP)

  static void affineHullMixed(std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    std::vector<std::shared_ptr<sparse_vector<mpq_class>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<mpq_class>>>& resultRays,
    std::vector<Constraint<mpq_class>>& resultEquations, const AffineHullQuery& query,
    const std::vector<Constraint<mpq_class>>& knownEquations, AffineHullResultCommon& result)
  {
    auto timeStarted = std::chrono::system_clock::now();
    auto timeComponent = std::chrono::system_clock::now();
    std::size_t n = polyhedron->space()->dimension();
    result.upperBound = n;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "affineHullMixed() called for ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto redundancyCheck = EquationRedundancyCheck<mpq_class>(n);
    auto approximateRedundancyCheck = EquationRedundancyCheck<double>(n);
    result.upperBound = filterGivenEquation(redundancyCheck, knownEquations, &resultEquations, 0,
      0, 0, &approximateRedundancyCheck);
    if (result.upperBound == -1)
    {
      result.dimension = -1;
      return;
    }

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "Initial " << knownEquations.size() << " equations have rank "
      << redundancyCheck.rank() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto affineComplement = AffineComplement<mpq_class>(n);
    auto approximateAffineComplement = AffineComplement<double>(n);
    std::vector<mpq_class> objective(n);
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
        findLastPoint(polyhedron, &objective[0], resultPoints,
          query.timeLimit - elapsedTime(timeStarted), result.timeOracles);
        result.dimension = result.lowerBound = result.upperBound;
        return;
      }

      sparse_vector<mpq_class> kernelDirectionVector;
      sparse_vector<double> approximateKernelDirectionVector;
      double kernelDirectionNorm = std::numeric_limits<double>::signaling_NaN();
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      mpq_class kernelDirectionValue = 0;
      while (kernelDirectionVector.empty())
      {
        assert(affineComplement.rank() == approximateAffineComplement.rank());

        timeComponent = std::chrono::system_clock::now();
        kernelDirectionColumn = approximateAffineComplement.selectColumn();

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Computing approximate kernel vector." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        approximateAffineComplement.computeKernelVector(kernelDirectionColumn,
          approximateKernelDirectionVector, query.epsilonCoefficient);
        result.timeKernel += elapsedTime(timeComponent);
        assert(!approximateKernelDirectionVector.empty());

        kernelDirectionNorm = euclideanNorm(approximateKernelDirectionVector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Checking approximate kernel vector." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        EquationRedundancy redundancy = approximateRedundancyCheck.test(
          approximateKernelDirectionVector, query.epsilonCoefficient);
        result.timeEquations += elapsedTime(timeComponent);

        if (redundancy != EQUATION_INDEPENDENT)
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Approximate kernel vector lies in the span." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          approximateAffineComplement.markEquation(kernelDirectionColumn);
          approximateKernelDirectionVector.clear();
          continue;
        }
        else
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Approximate kernel vector is independent." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        }

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Computing exact kernel vector." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        affineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector, 0);
        result.numKernel++;
        result.timeKernel += elapsedTime(timeComponent);
        assert(!kernelDirectionVector.empty());

        kernelDirectionNorm = euclideanNorm(kernelDirectionVector);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Checking exact kernel vector." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        if (elapsedTime(timeStarted) >= query.timeLimit)
        {
          result.dimension = AFFINEHULL_ERROR_TIMEOUT;
          return;
        }

        timeComponent = std::chrono::system_clock::now();
        redundancy = redundancyCheck.test(kernelDirectionVector, 0);
        result.timeEquations += elapsedTime(timeComponent);
  
        if (redundancy == EQUATION_INDEPENDENT)
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "Exact kernel vector is independent." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          break;
        }
        else
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "Exact kernel vector lies in the span." << std::endl;
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

      // Maximize kernelDirectionVector.

      int outcome = tryMaximization(polyhedron, affineComplement, &objective[0],
        kernelDirectionVector, kernelDirectionColumn, kernelDirectionNorm, kernelDirectionValue,
        resultPoints, resultRays, query.epsilonFactorization, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays,
        &approximateAffineComplement);
      if (outcome == INTERNAL_TIMEOUT)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_INFEASIBLE)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_INITIAL_ONE)
      {
        result.lowerBound++;
      }
      else if (outcome == INTERNAL_INITIAL_TWO)
      {
        result.lowerBound += 2;
        continue;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        result.lowerBound++;
        continue;
      }

      // Minimize kernelDirectionVector.

      outcome = tryMinimization(polyhedron, affineComplement, &objective[0],
        kernelDirectionVector, kernelDirectionColumn, kernelDirectionNorm, kernelDirectionValue,
        resultPoints, resultRays, query.epsilonFactorization, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays,
        &approximateAffineComplement);
      if (outcome == INTERNAL_TIMEOUT)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        result.lowerBound++;
        continue;
      }

      // Add equation defined by kernelDirectionVector.

      addEquation(redundancyCheck, affineComplement, resultEquations, kernelDirectionVector,
        kernelDirectionColumn, kernelDirectionValue, query.epsilonFactorization,
        query.epsilonConstraints, query.epsilonCoefficient, result.timeEquations,
        &approximateRedundancyCheck, &approximateAffineComplement);
      --result.upperBound;
    }

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
  }

  static void affineHullVerify(std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
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
    std::vector<Constraint<double>> approximateEquations;
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "affineHullVerify() called with ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto approximateRedundancyCheck = EquationRedundancyCheck<double>(n);
    auto exactRedundancyCheck = EquationRedundancyCheck<mpq_class>(n);
    result.upperBound = filterGivenEquation(approximateRedundancyCheck, knownEquations,
      &approximateEquations, query.epsilonConstraints, query.epsilonFactorization,
      query.epsilonCoefficient);
    if (result.upperBound == -1)
    {
      result.upperBound = filterGivenEquation(exactRedundancyCheck, knownEquations,
        &resultEquations, 0.0, 0.0, 0.0);
      if (result.upperBound == -1)
      {
        result.dimension = -1;
        return;
      }
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
      std::cout << "Iteration of approximate affine hull computation: " << candidateLowerBound
        << " <= dim <= " << candidateUpperBound << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }

      // If we have dim(P) rays but no points yet, we solve the feasibility problem.

      if (candidatePoints.empty() && int(candidateRays.size()) == candidateUpperBound)
      {
        findLastPoint(polyhedron, &approximateObjective[0], candidatePoints,
          query.timeLimit - elapsedTime(timeStarted), result.timeOracles);
        candidateLowerBound = candidateUpperBound;
        break;
      }

      sparse_vector<double> kernelDirectionVector;
      double kernelDirectionNorm = std::numeric_limits<double>::signaling_NaN();
      std::size_t kernelDirectionColumn = std::numeric_limits<std::size_t>::max();
      double kernelDirectionValue = 0.0;
      while (kernelDirectionVector.empty())
      {
        timeComponent = std::chrono::system_clock::now();
        kernelDirectionColumn = approximateAffineComplement.selectColumn();
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "  Column " << kernelDirectionColumn << std::flush;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        approximateAffineComplement.computeKernelVector(kernelDirectionColumn, kernelDirectionVector,
          query.epsilonCoefficient);
        result.numKernel++;
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

      // Maximize kernelDirectionVector.

      int outcome = tryMaximization(polyhedron, approximateAffineComplement,
        &approximateObjective[0], kernelDirectionVector, kernelDirectionColumn, kernelDirectionNorm,
        kernelDirectionValue, candidatePoints, candidateRays, query.epsilonFactorization,
        query.epsilonConstraints, query.timeLimit - elapsedTime(timeStarted), result.timeOracles,
        result.timePointsRays);
      if (outcome == INTERNAL_TIMEOUT)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_INFEASIBLE)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_INITIAL_ONE)
      {
        candidateLowerBound++;
      }
      else if (outcome == INTERNAL_INITIAL_TWO)
      {
        candidateLowerBound += 2;
        continue;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        candidateLowerBound++;
        continue;
      }

      // Minimize kernelDirectionVector.

      outcome = tryMinimization(polyhedron, approximateAffineComplement, &approximateObjective[0],
        kernelDirectionVector, kernelDirectionColumn, kernelDirectionNorm, kernelDirectionValue,
        candidatePoints, candidateRays, query.epsilonFactorization, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays);
      if (outcome == INTERNAL_TIMEOUT)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        candidateLowerBound++;
        continue;
      }

      // Add equation defined by kernelDirectionVector.

      addEquation(approximateRedundancyCheck, approximateAffineComplement, approximateEquations,
        kernelDirectionVector, kernelDirectionColumn, kernelDirectionValue,
        query.epsilonFactorization, query.epsilonConstraints, query.epsilonCoefficient,
        result.timeEquations);
      --candidateUpperBound;
    }

    assert(false);

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
  }

#endif /* IPO_WITH_GMP */

  AffineHullResult<double> affineHull(
    std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<double>>& knownEquations)
  {
    AffineHullResult<double> result;
    affineHullDirect(polyhedron, result.points, result.rays, result.equations, query,
      knownEquations, result);
    return result;
  }

#if defined(IPO_WITH_GMP)

  AffineHullResult<mpq_class> affineHull(
    std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<mpq_class>>& knownEquations)
  {
    AffineHullResult<mpq_class> result;
    if (query.algorithm == AFFINEHULL_ALGORITHM_DIRECT)
    {
      affineHullDirect(polyhedron, result.points, result.rays, result.equations, query,
        knownEquations, result);
    }
    else if (query.algorithm == AFFINEHULL_ALGORITHM_VERIFY)
    {
      affineHullVerify(polyhedron, result.points, result.rays, result.equations, query,
        knownEquations, result);
    }
    else if (query.algorithm == AFFINEHULL_ALGORITHM_MIXED)
    {
      affineHullMixed(polyhedron, result.points, result.rays, result.equations, query,
        knownEquations, result);
    }
    return result;
  }

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
