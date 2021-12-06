#include <ipo/affine_hull.hpp>

// #define IPO_DEBUG_AFFINE_HULL_CHECK // Uncomment to enable debug checks.
// // #define IPO_DEBUG_AFFINE_HULL_PRINT // Uncomment to print activity.

#include <ipo/arithmetic.hpp>

#include "lu.hpp"
#include "redundancy.hpp"

#include <iostream>
#include <algorithm>
#include <chrono>

namespace ipo
{
  static const int INTERNAL_NUMERICS = AFFINEHULL_ERROR_NUMERICS;
  static const int INTERNAL_TIMEOUT = AFFINEHULL_ERROR_TIMEOUT;
  static const int INTERNAL_INFEASIBLE = -1;
  static const int INTERNAL_OKAY = 0;
  static const int INTERNAL_NONE = 0;
  static const int INTERNAL_INITIAL_ONE = 1;
  static const int INTERNAL_INITIAL_TWO = 2;
  static const int INTERNAL_FURTHER_ONE = 3;

  AffineHullQuery::AffineHullQuery()
    : epsilonLinearDependence(1.0e-6), epsilonConstraints(1.0e-6), epsilonInvalidate(1.e-20),
    epsilonSafety(1.0e-6), epsilonFactorization(1.0e-20), epsilonCoefficient(1.0e-12),
    timeLimit(std::numeric_limits<double>::infinity())
  {

  }

  AffineHullQuery& AffineHullQuery::operator=(const AffineHullQuery& other)
  {
    epsilonConstraints = other.epsilonConstraints;
    epsilonInvalidate = other.epsilonInvalidate;
    epsilonSafety = other.epsilonSafety;
    epsilonCoefficient = other.epsilonCoefficient;
    epsilonFactorization = other.epsilonFactorization;
    timeLimit = other.timeLimit;
    return *this;
  }

  template <typename Number>
  AffineHullResult<Number>::AffineHullResult()
  {
    dimension = AFFINEHULL_ERROR_RUNNING;
    lowerBound = -1;
    upperBound = std::numeric_limits<int>::max();
    timeTotal = 0.0;
    timeOracles = 0.0;
    numKernel = 0;
    timeKernel = 0.0;
    timePointsRays = 0.0;
    timeEquations = 0.0;
  }

  template <typename Number>
  AffineHullResult<Number>::AffineHullResult(AffineHullResult<Number>&& other)
    : points(std::move(other.points)), rays(std::move(other.rays)),
    equations(std::move(other.equations))
  {
    dimension = other.dimension;
    lowerBound = other.lowerBound;
    upperBound = other.upperBound;
    timeTotal = other.timeTotal;
    timeOracles = other.timeOracles;
    numKernel = other.numKernel;
    timeKernel = other.timeKernel;
    timePointsRays = other.timePointsRays;
    timeEquations = other.timeEquations;
  }

  template <typename Number>
  AffineHullResult<Number>& AffineHullResult<Number>::operator=(const AffineHullResult<Number>& other)
  {
    points = other.points;
    rays = other.rays;
    equations = other.equations;
    dimension = other.dimension;
    lowerBound = other.lowerBound;
    upperBound = other.upperBound;
    timeTotal = other.timeTotal;
    timeOracles = other.timeOracles;
    numKernel = other.numKernel;
    timeKernel = other.timeKernel;
    timePointsRays = other.timePointsRays;
    timeEquations = other.timeEquations;
    return *this;
  }

  template <typename Number>
  AffineHullResult<Number>& AffineHullResult<Number>::operator=(AffineHullResult<Number>&& other)
  {
    points = std::move(other.points);
    rays = std::move(other.rays);
    equations = std::move(other.equations);
    dimension = other.dimension;
    lowerBound = other.lowerBound;
    upperBound = other.upperBound;
    timeTotal = other.timeTotal;
    timeOracles = other.timeOracles;
    numKernel = other.numKernel;
    timeKernel = other.timeKernel;
    timePointsRays = other.timePointsRays;
    timeEquations = other.timeEquations;
    return *this;
  }

  template class AffineHullResult<double>;

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

      if (countNonzeroProducts != 1)
      {
        T firstProd = 0;
        for (std::size_t p = 0; p < _debugDensePoints.size(); ++p)
        {
          T prod = 0;
          for (std::size_t v = 0; v < numVariables(); ++v)
            prod += _debugDensePoints[p][v] * kernelVector[v];
          if (p == 0)
            firstProd = prod;
          std::cout << "!!! Product for point " << p << " is " << firstProd << " + " 
            << (prod - firstProd) << "." << std::endl;
        }
        for (std::size_t r = 0; r < _debugDenseRays.size(); ++r)
        {
          T prod = 0;
          for (std::size_t v = 0; v < numVariables(); ++v)
            prod += _debugDenseRays[r][v] * kernelVector[v];
          std::cout << "!!! Product for ray " << r << " is " << prod << std::endl;
        }
      }
    }
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */

    int add(const sparse_vector<T>& row, const T& last, std::size_t newBasicColumn,
      double epsilonFactorization)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    add(row=" << row << ", last=" << last << ", newBasicColumn="
        << newBasicColumn << ")." << std::endl;      
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      // Add column to basis.
      assert(newBasicColumn < _columns.size());
      assert(!_columns[newBasicColumn].isBasic());
      _columns[newBasicColumn].basisIndex = _basisIndexToColumn.size();
      _basisIndexToColumn.push_back(newBasicColumn);

      // Copy part of new row to basis matrix.
      std::vector<T> newRow(rank(), T(0));
      for (const auto& iter : row)
      {
        if (_columns[iter.first].isBasic())
          newRow[_columns[iter.first].basisIndex] = iter.second;
      }
      if (last != 0 && _columns.back().isBasic())
        newRow[_columns.back().basisIndex] = last;

      // Copy new basic column.
      std::vector<T> newColumn(rank() - 1, T(0));
      for (std::size_t i = 0; i < _columns[newBasicColumn].rows.size(); ++i)
        newColumn[_columns[newBasicColumn].rows[i]] = _columns[newBasicColumn].entries[i];

      // Update LU decomposition.
      if (!_lu.extend(&newRow[0], &newColumn[0], newRow.back(), epsilonFactorization))
        return INTERNAL_NUMERICS;

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

      return INTERNAL_OKAY;
    }

    void computeKernelVector(std::size_t column, sparse_vector<T>& result, T& rhs,
      double epsilonCoefficient)
    {
      std::vector<T> vector(rank(), T(0));
      const ColumnData& columnData = _columns[column];
      for (std::size_t i = 0; i < columnData.rows.size(); ++i)
        vector[columnData.rows[i]] = columnData.entries[i];

      _lu.solveLeft(&vector[0]);
      _lu.solveUpper(&vector[0]);

      result.clear();
      for (std::size_t c = 0; c < numVariables(); ++c)
      {
        if (c == column)
          result.push_back(c, T(-1));
        else if (_columns[c].isBasic()
          && fabs(convertNumber<double>(vector[_columns[c].basisIndex])) > epsilonCoefficient)
        {
          assert(vector[_columns[c].basisIndex] != 0.0);
          result.push_back(c, vector[_columns[c].basisIndex]);
        }
      }
      if (column == numVariables())
        rhs = -1;
      else if (_columns.back().isBasic())
        rhs = -vector[_columns.back().basisIndex];
      else
        rhs = 0;
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

      bool isBasic() const
      {
        return basisIndex < std::numeric_limits<std::size_t>::max();
      }
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
    std::vector<Constraint<U>>* pResultEquations, double epsilonCoefficient,
    EquationRedundancyCheck<double>* pApproximateRedundancyCheck = NULL)
  {
    for (auto equation : knownEquations)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
      auto convertedEquation = convertConstraint<U>(equation);
      auto redundancy = redundancyCheck.test(convertedEquation, euclideanNorm(convertedEquation.vector()));
      if (redundancy.maxViolation < epsilonCoefficient && redundancy.rhs >= epsilonCoefficient)
      {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "  Given equation implies infeasibility with violation "
          << redundancy.maxViolation << " and rhs " << redundancy.rhs << ": " << equation
          << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        if (pResultEquations)
        {
          pResultEquations->clear();
          pResultEquations->push_back(neverSatisfiedConstraint<U>());
        }
        return INTERNAL_INFEASIBLE;
      }
      else if (redundancy.maxViolation > epsilonCoefficient)
      {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "  Given equation is independent with violation "
          << redundancy.maxViolation << " for coordinate " << redundancy.maxCoordinate << ": " 
          << equation << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        if (!redundancyCheck.add(convertedEquation, redundancy.maxCoordinate, 0.0))
          return INTERNAL_NUMERICS;
        if (pApproximateRedundancyCheck)
        {
          pApproximateRedundancyCheck->add(convertConstraint<double>(equation),
            redundancy.maxCoordinate, 0.0);
        }
        if (pResultEquations)
          pResultEquations->push_back(convertedEquation);
      }
    }
    return redundancyCheck.numVariables() - redundancyCheck.rank();
  }

  struct KernelVectorDouble
  {
    std::size_t column;
    std::shared_ptr<sparse_vector<double>> vector;
    double rhs;
    double norm;
    double maxViolation;
    std::size_t sparsity;
    std::size_t redundancyCoordinate;
    std::size_t redundancyRank;
    bool isLast;

    KernelVectorDouble() = default;

    KernelVectorDouble(std::size_t col, bool valid, bool last)
      : column(col), rhs(0.0), norm(1.0), maxViolation(std::numeric_limits<double>::max()),
      sparsity(1), redundancyCoordinate(std::numeric_limits<std::size_t>::max()),
      redundancyRank(std::numeric_limits<std::size_t>::max()), isLast(last)
    {
      if (valid && !last)
      {
        vector = std::make_shared<sparse_vector<double>>();
        vector->push_back(col, -1);
        redundancyCoordinate = col;
      }
    }

    void ensureValidity(const sparse_vector<double>& solution, double pointFactor, double epsilon)
    {
      if (vector)
      {
        double violation = fabs(*vector * solution - rhs * pointFactor);
        if (violation > epsilon * norm)
          vector.reset();
      }
    }
  };

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
  std::ostream& operator<<(std::ostream& stream, const KernelVectorDouble& kv)
  {
    stream << "Kernel vector for column " << kv.column;
    if (kv.vector)
    {
      double min = std::numeric_limits<double>::infinity();
      double max = 0.0;
      for (const auto& iter : *kv.vector)
      {
        double x = fabs(iter.second);
        min = std::min(min, x);
        max = std::max(max, x);
      }
      double quotient = (min > 0 && max < std::numeric_limits<double>::infinity())
        ? max/min : std::numeric_limits<double>::infinity();
      return stream << " with sparsity " << kv.sparsity
        << " and max violation " << kv.maxViolation << " attained at " << kv.redundancyCoordinate
        << " has condition number " << quotient << " and is " << *kv.vector << ".";
    }
    else
    {
      return stream << " with expected sparsity " << kv.sparsity
        << " and expected max violation " << kv.maxViolation << " attained at "
        << kv.redundancyCoordinate << " is invalid.";
    }
  }
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

  struct KernelVectorDoubleLess
  {
    double epsilon;

    bool operator()(const KernelVectorDouble& a, const KernelVectorDouble& b) const
    {
      if (a.isLast && !b.isLast)
        return true;
      if (!a.isLast && b.isLast)
        return false;

      if (a.maxViolation <= epsilon && b.maxViolation > epsilon)
        return true;
      if (a.maxViolation > epsilon && b.maxViolation <= epsilon)
        return false;

      if (a.maxViolation >= epsilon)
      {
        // If both have a safe violation, we prefer valid vectors, and among those the sparse ones.

        if (a.sparsity > b.sparsity)
          return true;
        if (a.sparsity < b.sparsity)
          return false;

        if (!a.vector && b.vector)
          return true;
        if (a.vector && !b.vector)
          return false;

        return false;
      }
      else
      {
        // If both don't have a safe violation we prefer invalid ones to validate them.
        // Among the valid ones, those with a largest violation are better.

        if (a.vector && !b.vector)
          return true;
        if (!a.vector && b.vector)
          return false;

        return a.maxViolation < b.maxViolation;
      }
    }
  };

#if defined(IPO_WITH_GMP)

  struct KernelVectorRational
  {
    std::size_t column;
    std::shared_ptr<sparse_vector<mpq_class>> vector;
    mpq_class rhs;
    std::shared_ptr<sparse_vector<double>> approximateVector;
    double approximateRhs;
    double norm;
    double maxViolation;
    std::size_t sparsity;
    std::size_t redundancyCoordinate;
    std::size_t redundancyRank;
    bool isLast;

    KernelVectorRational() = default;

    KernelVectorRational(std::size_t col, bool valid, bool last)
      : column(col), rhs(0), approximateRhs(1.0), norm(1.0),
      maxViolation(std::numeric_limits<double>::max()), sparsity(1),
      redundancyCoordinate(std::numeric_limits<std::size_t>::max()),
      redundancyRank(std::numeric_limits<std::size_t>::max()), isLast(last)
    {
      if (valid && !last)
      {
        vector = std::make_shared<sparse_vector<mpq_class>>();
        vector->push_back(col, mpq_class(-1));
        approximateVector = std::make_shared<sparse_vector<double>>();
        approximateVector->push_back(col, -1.0);        
        redundancyCoordinate = col;
      }
    }

    void ensureValidity(const sparse_vector<mpq_class>& solution, const mpq_class& pointFactor,
      double epsilon)
    {
      if (vector)
      {
        double violation = fabs(convertNumber<double, mpq_class>(
          *vector * solution - rhs * pointFactor));
        if (violation > epsilon * norm)
          vector.reset();
      }
    }
  };

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
  std::ostream& operator<<(std::ostream& stream, const KernelVectorRational& kv)
  {
    stream << "Kernel vector for column " << kv.column << ": ";
    if (kv.vector)
    {
      return stream << *kv.vector << " with sparsity " << kv.sparsity
        << " and max violation " << kv.maxViolation << " attained at " << kv.redundancyCoordinate;
    }
    else
    {
      return stream << "<invalid> with expected sparsity " << kv.sparsity
        << " and expected max violation " << kv.maxViolation;
    }
  }
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

  struct KernelVectorRationalLess
  {
    double epsilon;

    bool operator()(const KernelVectorRational& a, const KernelVectorRational& b) const
    {
      if (a.isLast && !b.isLast)
        return true;
      if (!a.isLast && b.isLast)
        return false;

      if (a.maxViolation <= epsilon && b.maxViolation > epsilon)
        return true;
      if (a.maxViolation > epsilon && b.maxViolation <= epsilon)
        return false;

      if (a.maxViolation >= epsilon)
      {
        // If both have a safe violation, we prefer valid vectors, and among those the sparse ones.

        if (a.sparsity > b.sparsity)
          return true;
        if (a.sparsity < b.sparsity)
          return false;

        if (!a.approximateVector && b.approximateVector)
          return true;
        if (a.approximateVector && !b.approximateVector)
          return false;

        if (!a.vector && b.vector)
          return true;
        if (a.vector && !b.vector)
          return false;

        return false;
      }
      else
      {
        // If both don't have a safe violation we prefer invalid ones to validate them.
        // Among the valid ones, those with a largest violation are better.

        if (a.approximateVector && !b.approximateVector)
          return true;
        if (!a.approximateVector && b.approximateVector)
          return false;

        return a.maxViolation < b.maxViolation;
      }
    }
  };

#endif /* IPO_WITH_GMP */

  template <typename P, typename T, typename U>
  static void findLastPoint(std::shared_ptr<P> polyhedron, U* objective,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints, double timeLimit,
    double& timeOracles)
  {
    typename P::OptOracle::Query oracleQuery;
    for (std::size_t v = 0; v < polyhedron->space()->dimension(); ++v)
      objective[v] = 0;
    oracleQuery.timeLimit = timeLimit;
    auto timeComponent = std::chrono::system_clock::now();
    auto oracleResponse = polyhedron->maximize(objective, oracleQuery);
    timeOracles += elapsedTime(timeComponent);

    if (oracleResponse.outcome == OptimizationOutcome::UNBOUNDED)
      throw std::runtime_error("maximize(<zero vector>) returned unbounded result.");
    else if (oracleResponse.outcome == OptimizationOutcome::FEASIBLE)
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

  static void stabilizeDirection(std::size_t n,
    const std::vector<std::shared_ptr<sparse_vector<double>>>& points,
    const std::vector<std::shared_ptr<sparse_vector<double>>>& rays,
    double epsilon, std::shared_ptr<sparse_vector<double>>& input, double& rhs)
  {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    Stabilizing vector of size " << input->size() << ": " << *input << "."
        << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    std::vector<double> sortedEntries;
    sortedEntries.reserve(input->size());
    for (const auto& iter : *input)
      sortedEntries.push_back(fabs(iter.second));
    std::sort(sortedEntries.begin(), sortedEntries.end());

    double threshold = sortedEntries.front();
    double lastThreshold = threshold;
    assert(threshold > 0);
    std::size_t lastSupport = std::numeric_limits<std::size_t>::max();
    std::vector<double> dense(n, 0.0);
    double lastAverage = rhs;
    while (lastSupport > 0)
    {
      threshold *= 4.0;
      std::size_t support = 0;
      for (const auto& iter : *input)
      {
        if (fabs(iter.second) > threshold)
        {
          dense[iter.first] = iter.second;
          ++support;
        }
        else
          dense[iter.first] = 0.0;
      }
      if (support == 0)
        break;
      if (support == lastSupport)
        continue;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "      Testing sparsified vector with threshold " << threshold << " and support "
        << support << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      double norm = euclideanNorm(&dense[0], n);
      if (norm == 0.0)
        continue;

      // Check points.

      bool failed = false;
      double min = std::numeric_limits<double>::infinity();
      double max = -std::numeric_limits<double>::infinity();
      for (const auto& point : points)
      {
        double prod = *point * &dense[0];
        min = std::min(min, prod);
        max = std::max(max, prod);
        if (max - min > 2 * norm * epsilon)
        {
          failed = true;
          break;
        }
      }

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "      Normalized discrepancy is " << ((max - min) / norm) << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      // Check rays.

      if (!failed)
      {
        for (const auto& ray : rays)
        {
          double prod = *ray * &dense[0];
          if (fabs(prod) / euclideanNorm(*ray) > epsilon)
          {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
            std::cout << "      Not orthogonal to some ray." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
            failed = true;
            break;
          }
        }
      }

      if (failed)
        break;

      lastThreshold = threshold;
      lastAverage = 0.5 * (min + max);
      lastSupport = support;
    }

    auto output = std::make_shared<sparse_vector<double>>();
    for (const auto& iter : *input)
    {
      if (fabs(iter.second) >= lastThreshold)
        output->push_back(iter.first, iter.second);
    }
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    if (input->size() != output->size())
    {
      std::cout << "    Stabilized vector has size " << output->size() << ": " << *output << "."
        << std::endl;
    }
    else
      std::cout << "    Did not stabilize." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    input = output;
    rhs = lastAverage;
  }

  template <typename T, typename KV>
  static void ensureValidity(std::vector<KV>& kernelVectors, const sparse_vector<T>& solution,
    const T& pointFactor, double epsilon, AffineComplement<T>& affineComplement,
    const std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints)
  {
    for (auto& kv : kernelVectors)
    {
      sparse_vector<T> oldVector;
      sparse_vector<T> newVector;
      T debugKernelRhs;
      if (kv.vector)
      {
        oldVector = *kv.vector;
        affineComplement.computeKernelVector(kv.column, newVector, debugKernelRhs, epsilon);
      }
      kv.ensureValidity(solution, pointFactor, epsilon);

#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      if (kv.vector && euclideanDistance(oldVector, newVector) > 1.e-7)
      {
        std::cout << "\n Kernel vector that was stored as valid was updated.\n";
        for (auto& vector : resultPoints)
        {
          std::cout << "Product of point with old vector is " << (*vector * oldVector)
            << " and with new vector: " << (*vector * newVector) << std::endl;
        }
        std::cout << "kv.column = " << kv.column << ", old vector is " << oldVector
          << ", new vector is " << newVector << std::endl;
        std::cout << "kv.rhs = " << kv.rhs << std::endl;
        std::cout << "Euclidean distance between the vectors is "
          << euclideanDistance(oldVector, newVector) << std::endl;
        assert(false);
      }
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
    }
  }

  template <typename P, typename T, typename U, typename KV>
  static int tryMaximization(std::shared_ptr<P> polyhedron,
    AffineComplement<U>& affineComplement, U* objective, std::vector<KV>& kernelVectors,
    const sparse_vector<U>& kernelVector, std::size_t kernelVectorColumn, double kernelVectorNorm,
    U& kernelVectorValue, std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultRays, const AffineHullQuery& query,
    double epsilonConstraints, double timeLimit, double& timeOracles, double& timePointsRays,
    AffineComplement<double>* pAffineComplement = NULL)
  {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "    Maximization" << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    for (std::size_t v = 0; v < polyhedron->space()->dimension(); ++v)
      objective[v] = 0;
    kernelVectorValue = 0;
    for (const auto& iter : kernelVector)
      objective[iter.first] = iter.second;

    typename P::OptOracle::Query oracleQuery;
    if (!resultPoints.empty())
    {
      // We set a minimum required objective value to the common value plus an epsilon for
      // (normalized) constraints.

      kernelVectorValue = kernelVector * *resultPoints.front();
      oracleQuery.setMinPrimalBound(kernelVectorValue + epsilonConstraints * kernelVectorNorm);
      // TODO: We actually desire a much bigger value.
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    Common objective of previous points is " << kernelVectorValue << "."
        << " Requiring " << oracleQuery.minPrimalBound() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
    }
    if (timeLimit <= 0)
      return INTERNAL_TIMEOUT;
    oracleQuery.timeLimit = timeLimit;

    auto timeComponent = std::chrono::system_clock::now();
    auto oracleResponse = polyhedron->maximize(objective, oracleQuery);
    double timeOracle = elapsedTime(timeComponent);
    timeOracles += timeOracle;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "    Result: " << oracleResponse << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    if (timeOracle >= timeLimit || oracleResponse.outcome == OptimizationOutcome::TIMEOUT)
      return INTERNAL_TIMEOUT;

    if (oracleResponse.outcome == OptimizationOutcome::INFEASIBLE)
      return INTERNAL_INFEASIBLE;
    else if (oracleResponse.outcome == OptimizationOutcome::UNBOUNDED)
    {
      const auto& firstRay = oracleResponse.rays.front();
      resultRays.push_back(firstRay.vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertTo<U>(*firstRay.vector), 0,
        kernelVectorColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      int error = affineComplement.add(convertSparseVector<U>(*firstRay.vector), 0,
        kernelVectorColumn, query.epsilonFactorization);
      if (error != INTERNAL_OKAY)
        return error;
      
      if (pAffineComplement)
      {
        error = pAffineComplement->add(convertSparseVector<double>(*firstRay.vector), 0,
          kernelVectorColumn, query.epsilonFactorization);
        if (error != INTERNAL_OKAY)
          return error;
      }
      timePointsRays += elapsedTime(timeComponent);

      // Invalidate kernel vectors.
      ensureValidity(kernelVectors, *firstRay.vector, T(0), firstRay.norm() * query.epsilonInvalidate, affineComplement, resultPoints);

      // Process points if present.

      if (resultPoints.empty() && !oracleResponse.points.empty())
      {
        const auto& bestPoint = oracleResponse.points.front();
        resultPoints.push_back(bestPoint.vector);
        timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
        affineComplement._debugAdd(convertSparseVector<U>(*bestPoint.vector), 1,
          polyhedron->space()->dimension(), objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
        int error = affineComplement.add(convertSparseVector<U>(*bestPoint.vector), 1,
          polyhedron->space()->dimension(), query.epsilonFactorization);
        if (error != INTERNAL_OKAY)
          return error;
        if (pAffineComplement)
        {
          error = pAffineComplement->add(convertSparseVector<double>(*bestPoint.vector), 1,
            polyhedron->space()->dimension(), query.epsilonFactorization);
          if (error != INTERNAL_OKAY)
            return error;
        }
        timePointsRays += elapsedTime(timeComponent);
        
        // Invalidate kernel vectors.
        ensureValidity(kernelVectors, *bestPoint.vector, T(1), query.epsilonInvalidate,
          affineComplement, resultPoints);

        return INTERNAL_INITIAL_TWO;
      }
      return resultPoints.empty() ? INTERNAL_INITIAL_ONE : INTERNAL_FURTHER_ONE;
    }
    else if (resultPoints.empty())
    {
      assert(oracleResponse.outcome == OptimizationOutcome::FEASIBLE);

      if (oracleResponse.points.empty())
      {
        throw std::runtime_error(
          "Detected invalid oracle behavior: responded with feasible polyhedron but no points.");
      }
      assert(!oracleResponse.points.empty());

      // Add the maximizer.

      const auto& bestPoint = oracleResponse.points.front();
      resultPoints.push_back(bestPoint.vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertSparseVector<U>(*bestPoint.vector), 1,
        polyhedron->space()->dimension(), objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      int error = affineComplement.add(convertSparseVector<U>(*bestPoint.vector), 1,
        polyhedron->space()->dimension(), query.epsilonFactorization);
      if (error != INTERNAL_OKAY)
        return error;
      if (pAffineComplement)
      {
        error = pAffineComplement->add(convertSparseVector<double>(*bestPoint.vector), 1,
          polyhedron->space()->dimension(), query.epsilonFactorization);
        if (error != INTERNAL_OKAY)
          return error;
      }
      timePointsRays += elapsedTime(timeComponent);

      // Invalidate kernel vectors.
      ensureValidity(kernelVectors, *bestPoint.vector, T(1), query.epsilonInvalidate, affineComplement, resultPoints);

      // If another point has a different objective value, we also add that and continue.
      for (std::size_t p = oracleResponse.points.size() - 1; p > 0; --p)
      {
        const auto& otherPoint = oracleResponse.points[p];
        double difference = fabs(convertNumber<double, T>(otherPoint.objectiveValue
          - bestPoint.objectiveValue));
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "    Two initial points found. Absolute objective value difference: "
          << difference << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        difference /= kernelVectorNorm;
        // TODO: In fact we desire a much larger difference.
        if (difference > 2 * epsilonConstraints)
        {
          resultPoints.push_back(otherPoint.vector);
          timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
          affineComplement._debugAdd(convertSparseVector<U>(*otherPoint.vector), 1, kernelVectorColumn,
            objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
          int error = affineComplement.add(convertSparseVector<U>(*otherPoint.vector), 1, kernelVectorColumn,
            query.epsilonFactorization);
          if (error != INTERNAL_OKAY)
            return error;
          if (pAffineComplement)
          {
            error = pAffineComplement->add(convertSparseVector<double>(*otherPoint.vector), 1,
              kernelVectorColumn, query.epsilonFactorization);
            if (error != INTERNAL_OKAY)
              return error;
          }
          timePointsRays += elapsedTime(timeComponent);

          // Invalidate kernel vectors.
          ensureValidity(kernelVectors, *otherPoint.vector, T(1), query.epsilonInvalidate, affineComplement, resultPoints);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "    -> adding two points with objective values "
            << convertNumber<double>(bestPoint.objectiveValue) << " and "
            << convertNumber<double>(otherPoint.objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          return INTERNAL_INITIAL_TWO;
        }
      }
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    -> adding a point with objective value "
        << convertNumber<double>(bestPoint.objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
      kernelVectorValue = convertNumber<U>(bestPoint.objectiveValue);
      return INTERNAL_INITIAL_ONE;
    }
    else
    {
      assert(oracleResponse.outcome == OptimizationOutcome::FEASIBLE);

      // Check if some vector has a different objective value.
      for (const auto& point : oracleResponse.points)
      {
        assert(oracleQuery.hasMinPrimalBound());
        double difference = fabs(convertNumber<double, T>(
          point.objectiveValue - kernelVectorValue));
        difference /= kernelVectorNorm;
        // TODO: In fact we desire a larger difference.
        if (difference > 2 * epsilonConstraints)
        {
          resultPoints.push_back(point.vector);
          timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
          affineComplement._debugAdd(convertSparseVector<U>(*point.vector), 1, kernelVectorColumn,
            objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
          int error = affineComplement.add(convertSparseVector<U>(*point.vector), 1, kernelVectorColumn,
            query.epsilonFactorization);
          if (error != INTERNAL_OKAY)
            return error;
          if (pAffineComplement)
          {
            error = pAffineComplement->add(convertSparseVector<double>(*point.vector), 1,
              kernelVectorColumn, query.epsilonFactorization);
            if (error != INTERNAL_OKAY)
              return error;
          }
          timePointsRays += elapsedTime(timeComponent);

          // Invalidate kernel vectors.
          ensureValidity(kernelVectors, *point.vector, T(1), query.epsilonInvalidate, affineComplement, resultPoints);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "    -> adding a point with objective value "
            << convertNumber<double>(oracleResponse.points.front().objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
          return INTERNAL_FURTHER_ONE;
        }
      }
      return INTERNAL_NONE;
    }
    return std::numeric_limits<int>::max();
  }

  template <typename P, typename T, typename U, typename KV>
  static int tryMinimization(std::shared_ptr<P> polyhedron,
    AffineComplement<U>& affineComplement, U* objective, std::vector<KV>& kernelVectors,
    const sparse_vector<U>& kernelVector, std::size_t kernelVectorColumn,
    double kernelVectorNorm, const U& kernelVectorValue,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultPoints,
    std::vector<std::shared_ptr<sparse_vector<T>>>& resultRays, const AffineHullQuery& query,
    double epsilonConstraints, double timeLimit, double& timeOracles, double& timePointsRays,
    AffineComplement<double>* pAffineComplement = NULL)
  {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "    Minimization" << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    for (std::size_t v = 0; v < polyhedron->space()->dimension(); ++v)
    {
      objective[v] = -objective[v];
#if !defined(NDEBUG)
      assert(objective[v] == -kernelVector.find(v, 0));
#endif /* !NDEBUG */
    }

    typename P::OptOracle::Query oracleQuery;
    if (!resultPoints.empty())
    {
      oracleQuery.setMinPrimalBound(-kernelVectorValue + epsilonConstraints * kernelVectorNorm);
      // TODO: We actually desire a much larger value.
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    Common objective of previous points is " << -kernelVectorValue << "."
        << " Requiring " << oracleQuery.minPrimalBound() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    }
    if (timeLimit <= 0)
      return INTERNAL_TIMEOUT;
    oracleQuery.timeLimit = timeLimit;

    auto timeComponent = std::chrono::system_clock::now();
    auto oracleResponse = polyhedron->maximize(objective, oracleQuery);
    double timeOracle = elapsedTime(timeComponent);
    timeOracles += timeOracle;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "    Result: " << oracleResponse << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    if (timeOracle >= timeLimit)
      return INTERNAL_TIMEOUT;

    if (oracleResponse.outcome == OptimizationOutcome::INFEASIBLE)
      throw std::runtime_error("Oracle for minimization claims infeasible.");
    else if (oracleResponse.outcome == OptimizationOutcome::UNBOUNDED)
    {
      const auto& firstRay = oracleResponse.rays.front();
      resultRays.push_back(firstRay.vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertSparseVector<U>(*firstRay.vector), 0,
        kernelVectorColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      int error = affineComplement.add(convertSparseVector<U>(*firstRay.vector), 0,
        kernelVectorColumn, query.epsilonFactorization);
      if (error != INTERNAL_OKAY)
        return error;
      if (pAffineComplement)
      {
        error = pAffineComplement->add(convertSparseVector<double>(*firstRay.vector), 0,
          kernelVectorColumn, query.epsilonFactorization);
        if (error != INTERNAL_OKAY)
          return error;
      }
      timePointsRays += elapsedTime(timeComponent);

      // Invalidate kernel vectors.
      ensureValidity(kernelVectors, *firstRay.vector, T(0), firstRay.norm() * query.epsilonInvalidate, affineComplement, resultPoints);

      return INTERNAL_FURTHER_ONE;
    }
    else
    {
      assert(!resultPoints.empty());
      assert(oracleResponse.outcome == OptimizationOutcome::FEASIBLE);

      // Check if some vector has a different objective value.
      if (oracleResponse.points.empty())
        return INTERNAL_NONE;

      const auto& bestPoint = oracleResponse.points.front();
      double difference = fabs(convertNumber<double, T>(
        bestPoint.objectiveValue - -kernelVectorValue));
      difference /= kernelVectorNorm;
      if (difference <= epsilonConstraints)
        return INTERNAL_NONE;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    Normalized objective value difference is " << difference << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      resultPoints.push_back(bestPoint.vector);
      timeComponent = std::chrono::system_clock::now();
#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      affineComplement._debugAdd(convertSparseVector<U>(*bestPoint.vector), 1,
        kernelVectorColumn, objective);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */
      int error = affineComplement.add(convertSparseVector<U>(*bestPoint.vector), 1,
        kernelVectorColumn, query.epsilonFactorization);
      if (error != INTERNAL_OKAY)
        return error;
      if (pAffineComplement)
      {
        error = pAffineComplement->add(convertSparseVector<double>(*bestPoint.vector), 1,
          kernelVectorColumn, query.epsilonFactorization);
        if (error != INTERNAL_OKAY)
          return error;
      }
      timePointsRays += elapsedTime(timeComponent);

      // Invalidate kernel vectors.
      ensureValidity(kernelVectors, *bestPoint.vector, T(1), query.epsilonInvalidate, affineComplement, resultPoints);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    -> adding a point with objective value "
          << convertNumber<double>(bestPoint.objectiveValue) << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
      return INTERNAL_FURTHER_ONE;
    }
    return std::numeric_limits<int>::max();
  }

  template <typename T>
  static int addEquation(EquationRedundancyCheck<T>& redundancyCheck,
    std::size_t redundancyCoordinate, AffineComplement<T>& affineComplement,
    std::vector<Constraint<T>>& resultEquations,
    std::shared_ptr<sparse_vector<T>> kernelDirectionVector, std::size_t kernelDirectionColumn,
    const T& kernelDirectionValue, double epsilonFactorization,
    double epsilonCoefficient, double& timeEquations,
    EquationRedundancyCheck<double>* pRedundancyCheck = NULL,
    AffineComplement<double>* pAffineComplement = NULL)
  {
    resultEquations.push_back(Constraint<T>(kernelDirectionValue, kernelDirectionVector, kernelDirectionValue));
    affineComplement.markEquation(kernelDirectionColumn);
    if (pAffineComplement)
      pAffineComplement->markEquation(kernelDirectionColumn);
    auto timeComponent = std::chrono::system_clock::now();
    if (!redundancyCheck.add(resultEquations.back(), redundancyCoordinate, epsilonFactorization))
      return INTERNAL_NUMERICS;
    if (pRedundancyCheck)
    {
      pRedundancyCheck->add(convertConstraint<double>(resultEquations.back()),
        redundancyCoordinate, epsilonFactorization);
    }
    timeEquations += elapsedTime(timeComponent);
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "    -> adding equation " << (resultEquations.size() - 1) << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    return INTERNAL_OKAY;
  }

  template <>
  AffineHullResult<double> affineHull(std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<double>>& knownEquations)
  {
    AffineHullResult<double> result;
    auto timeStarted = std::chrono::system_clock::now();
    auto timeComponent = std::chrono::system_clock::now();
    std::size_t n = polyhedron->space()->dimension();
    result.upperBound = n;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "affineHull<double>() called for ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto redundancyCheck = EquationRedundancyCheck<double>(n);
    int filterResult =  filterGivenEquation(redundancyCheck, knownEquations, &result.equations,
      query.epsilonCoefficient);
    if (filterResult == INTERNAL_INFEASIBLE)
    {
      result.upperBound = -1;
      result.dimension = -1;
      return result;
    }
    else if (filterResult == INTERNAL_NUMERICS)
    {
      result.dimension = AFFINEHULL_ERROR_NUMERICS;
      return result;
    }
    else
    {
      result.upperBound = filterResult;
    }

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "  Initial " << knownEquations.size() << " equations have rank "
      << redundancyCheck.rank() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto affineComplement = AffineComplement<double>(n);
    std::vector<double> objective(n);
    KernelVectorDoubleLess kernelVectorLess;
    kernelVectorLess.epsilon = query.epsilonLinearDependence;
    std::vector<KernelVectorDouble> kernelVectors;
    for (std::size_t v = 0; v <= n; ++v)
      kernelVectors.push_back(KernelVectorDouble(v, redundancyCheck.rank() == 0, v == n));

    // Main iteration loop.

    while (result.lowerBound < result.upperBound)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "  Iteration: " << result.lowerBound << " <= dim <= " << result.upperBound
        << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      verifyAffineHullResult(polyhedron, result);
#endif /* IPO_DEBUG_AFFINE_HULL_CHECK */

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return result;
      }

      // If we have dim(P) rays but no points yet, we solve the feasibility problem.

      if (result.points.empty() && int(result.rays.size()) == result.upperBound)
      {
        findLastPoint(polyhedron, &objective[0], result.points,
          query.timeLimit - elapsedTime(timeStarted), result.timeOracles);
        result.dimension = result.lowerBound = result.upperBound;
        return result;
      }

      std::make_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);
      while (true)
      {
        std::pop_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "    Top priority kernel vector: " << kernelVectors.back() << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        bool updateVector = !kernelVectors.back().vector;
        if (updateVector)
        {
          timeComponent = std::chrono::system_clock::now();
          kernelVectors.back().vector = std::make_shared<sparse_vector<double>>();
          affineComplement.computeKernelVector(kernelVectors.back().column,
            *kernelVectors.back().vector, kernelVectors.back().rhs, query.epsilonCoefficient);
          result.numKernel++;
          result.timeKernel += elapsedTime(timeComponent);
          assert(!kernelVectors.back().vector->empty());

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "    Computed kernel vector is " << *kernelVectors.back().vector
            << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

          stabilizeDirection(n, result.points, result.rays, query.epsilonConstraints,
            kernelVectors.back().vector, kernelVectors.back().rhs);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "    Stabilized kernel vector is " << *kernelVectors.back().vector
            << " with rhs " << kernelVectors.back().rhs << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

            
          if (elapsedTime(timeStarted) >= query.timeLimit)
          {
            result.dimension = AFFINEHULL_ERROR_TIMEOUT;
            return result;
          }

          kernelVectors.back().sparsity = kernelVectors.back().vector->size();
          kernelVectors.back().norm = euclideanNorm(*kernelVectors.back().vector);
        }

        bool updateRedundancy = updateVector
          || kernelVectors.back().redundancyRank != redundancyCheck.rank();
        if (updateRedundancy)
        {
          timeComponent = std::chrono::system_clock::now();
          auto redundancy = redundancyCheck.test(*kernelVectors.back().vector,
            kernelVectors.back().norm);
          result.timeEquations += elapsedTime(timeComponent);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "    Redundancy checked: violation is " << redundancy.maxViolation 
            << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

          if (elapsedTime(timeStarted) >= query.timeLimit)
          {
            result.dimension = AFFINEHULL_ERROR_TIMEOUT;
            return result;
          }

          kernelVectors.back().maxViolation = redundancy.maxViolation;
          kernelVectors.back().redundancyCoordinate = redundancy.maxCoordinate;
          kernelVectors.back().redundancyRank = redundancyCheck.rank();
        }
        else
        {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "    Redundancy must NOT be checked." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
        }

        if (!updateRedundancy)
          break;

        std::push_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);
      }
      KernelVectorDouble kernelVector = kernelVectors.back();
      kernelVectors.pop_back();

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    Selected kernel vector: " << kernelVector << std::endl;
      sparse_vector<double> debugKernelVector;
      double debugKernelRhs;
      affineComplement.computeKernelVector(kernelVector.column, debugKernelVector, debugKernelRhs, query.epsilonFactorization);
      std::cout << "    Verification reads: " << debugKernelVector << std::endl;
      std::cout << "    Product with all points are:";
      for (const auto& point : result.points)
      {
        std::cout << " " << (*point * debugKernelVector);
      }
      std::cout << "\n    Expected right-hand side is " << kernelVector.rhs << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      double kernelVectorValue = 0.0;

      // Maximize kernelVectorValue.

      int outcome = tryMaximization(polyhedron, affineComplement, &objective[0],
        kernelVectors, *kernelVector.vector, kernelVector.column, kernelVector.norm,
        kernelVectorValue, result.points, result.rays, query, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays);
      if (outcome < INTERNAL_OKAY)
      {
        result.dimension = outcome;
        return result;
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
        kernelVectors, *kernelVector.vector, kernelVector.column, kernelVector.norm,
        kernelVectorValue, result.points, result.rays, query, query.epsilonConstraints,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays);
      if (outcome < INTERNAL_OKAY)
      {
        result.dimension = outcome;
        return result;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        result.lowerBound++;
        continue;
      }

      // Add equation defined by kernelDirectionVector.

#if defined(IPO_DEBUG_AFFINE_HULL_CHECK)
      std::vector<double> direction(n, 0.0);
      for (const auto& iter : *kernelVector.vector)
        direction[iter.first] = iter.second;
      RealOptimizationOracle::Query verifyQuery;
      auto verifyResult = polyhedron->maximize(&direction[0], verifyQuery);
      std::cout << "    Maximization of equation normal yields " << verifyResult << std::endl;
      for (std::size_t v = 0; v < n; ++v)
        direction[v] = -direction[v];
      verifyResult = polyhedron->maximize(&direction[0], verifyQuery);
      std::cout << "    Maximization of negated equation normal yields " << verifyResult << std::endl;
      auto redTest = redundancyCheck.test(*kernelVector.vector, 1);
      std::cout << "    Repeated redundancy check yields unnormalized violation " << redTest.maxViolation << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      assert(kernelVector.redundancyCoordinate < n);
      int error = addEquation(redundancyCheck, kernelVector.redundancyCoordinate, affineComplement,
        result.equations, kernelVector.vector, kernelVector.column, kernelVectorValue,
        query.epsilonFactorization, query.epsilonCoefficient,
        result.timeEquations);
      if (error != INTERNAL_OKAY)
      {
        result.dimension = error;
        return result;
      }
      --result.upperBound;
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "    ";
      polyhedron->space()->printConstraint(std::cout, result.equations.back());
      std::cout << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */
    }

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
    return result;
  }

  bool verifyAffineHullResult(std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullResult<double>& result)
  {
    std::size_t n = polyhedron->space()->dimension();

    // Check that points are affinely independent.

    std::vector<std::vector<double>> matrix(result.points.size() + result.rays.size());
    for (std::size_t p = 0; p < result.points.size(); ++p)
    {
      matrix[p].resize(n+1, 0.0);
      for (const auto iter : *result.points[p])
        matrix[p][iter.first] = iter.second;
      matrix[p].back() = 1.0;
    }
    for (std::size_t r = 0; r < result.rays.size(); ++r)
    {
      matrix[result.points.size() + r].resize(n+1, 0.0);
      for (const auto iter : *result.rays[r])
        matrix[result.points.size() + r][iter.first] = iter.second;
    }

    int innerDimension = rowEchelon(n+1, matrix) - 1;
    if (innerDimension != result.lowerBound)
    {
      std::cerr << "The " << matrix.size() << " points and rays have affine dimension "
        << innerDimension << std::endl;
      return false;
    }

    // Check that equations are independent.

    matrix.clear();
    matrix.resize(result.equations.size());
    for (std::size_t e = 0; e < result.equations.size(); ++e)
    {
      matrix[e].resize(n, 0.0);
      for (const auto& iter : result.equations[e].vector())
        matrix[e][iter.first] = iter.second;
    }

    int outerDimension = n - rowEchelon(n, matrix);
    if (outerDimension != result.upperBound)
    {
      std::cerr << "The " << matrix.size() << " equations form a system of rank "
        << (n - outerDimension) << " and define a subspace of dimension " << outerDimension
        << "." << std::endl;
      return false;
    }

    // Check satisfaction of equations.

    std::size_t worstEquation = 0;
    std::size_t worstPoint = std::numeric_limits<std::size_t>::max();
    std::size_t worstRay = std::numeric_limits<std::size_t>::max();
    double worstViolation = 0.0;
    for (std::size_t e = 0; e < result.equations.size(); ++e)
    {
      const Constraint<double>& equation = result.equations[e];
      if (equation.type() != ConstraintType::EQUATION)
      {
        std::cerr << "Constraint " << e << " is not an equation!" << std::endl;
        return false;
      }
      for (std::size_t p = 0; p < result.points.size(); ++p)
      {
        double violation = fabs(equation.vector() * *result.points[p] - equation.rhs());
        violation /= euclideanNorm(equation.vector());
        if (violation > worstViolation)
        {
          worstViolation = violation;
          worstEquation = e;
          worstPoint = p;
          worstRay = std::numeric_limits<std::size_t>::max();
        }
      }
      for (std::size_t r = 0; r < result.rays.size(); ++r)
      {
        double violation = fabs(equation.vector() * *result.rays[r]);
        violation /= euclideanNorm(equation.vector()) * euclideanNorm(*result.rays[r]);
        if (violation > worstViolation)
        {
          worstViolation = violation;
          worstEquation = e;
          worstPoint = std::numeric_limits<std::size_t>::max();
          worstRay = r;
        }
      }
    }
    if (worstViolation > 1.0e-6)
    {
      std::cerr << "Worst violation attained by ";
      if (worstPoint < std::numeric_limits<std::size_t>::max())
        std::cerr << "point " << worstPoint << " (" << *result.points[worstPoint] << ")";
      else
        std::cerr << "ray " << worstRay << " (" << *result.rays[worstRay] << ")";
      std::cerr << " and equation " << worstEquation << " (" << result.equations[worstEquation]
        << ") with violation " << worstViolation << "." << std::endl;
    }
    return worstViolation <= 1.0e-6;
  }

  template AffineHullResult<double> affineHull(
    std::shared_ptr<Polyhedron<double>> polyhedron,
    const AffineHullQuery& query,
    const std::vector<Constraint<double>>& knownEquations);

#if defined(IPO_WITH_GMP)

  template class AffineHullResult<mpq_class>;

  template <>
  IPO_EXPORT
  AffineHullResult<mpq_class> affineHull(std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    const AffineHullQuery& query, const std::vector<Constraint<mpq_class>>& knownEquations)
  {
    AffineHullResult<mpq_class> result;
    auto timeStarted = std::chrono::system_clock::now();
    auto timeComponent = std::chrono::system_clock::now();
    std::size_t n = polyhedron->space()->dimension();
    result.upperBound = n;

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "affineHull<Rational>() called for ambient dimension " << n << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto redundancyCheck = EquationRedundancyCheck<mpq_class>(n);
    auto approximateRedundancyCheck = EquationRedundancyCheck<double>(n);
    int filterResult = filterGivenEquation(redundancyCheck, knownEquations, &result.equations, 0,
      &approximateRedundancyCheck);
    if (filterResult == INTERNAL_INFEASIBLE)
    {
      result.upperBound = -1;
      result.dimension = -1;
      return result;
    }
    else if (filterResult == INTERNAL_NUMERICS)
    {
      result.dimension = AFFINEHULL_ERROR_NUMERICS;
      return result;
    }
    else
    {
      result.upperBound = filterResult;
    }

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
    std::cout << "Initial " << knownEquations.size() << " equations have rank "
      << redundancyCheck.rank() << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

    auto affineComplement = AffineComplement<mpq_class>(n);
    auto approximateAffineComplement = AffineComplement<double>(n);
    std::vector<mpq_class> objective(n);
    KernelVectorRationalLess kernelVectorLess;
    kernelVectorLess.epsilon = query.epsilonLinearDependence;
    std::vector<KernelVectorRational> kernelVectors;
    for (std::size_t v = 0; v <= n; ++v)
      kernelVectors.push_back(KernelVectorRational(v, redundancyCheck.rank() == 0, v == n));
    while (result.lowerBound < result.upperBound)
    {
#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "Iteration of affine hull computation: " << result.lowerBound
        << " <= dim <= " << result.upperBound << "." << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      if (elapsedTime(timeStarted) >= query.timeLimit)
      {
        result.dimension = AFFINEHULL_ERROR_TIMEOUT;
        return result;
      }

      // If we have dim(P) rays but no points yet, we solve the feasibility problem.

      if (result.points.empty() && int(result.rays.size()) == result.upperBound)
      {
        findLastPoint(polyhedron, &objective[0], result.points,
          query.timeLimit - elapsedTime(timeStarted), result.timeOracles);
        result.dimension = result.lowerBound = result.upperBound;
        return result;
      }

      std::make_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);
      while (true)
      {
        std::pop_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
        std::cout << "Top priority kernel vector: " << kernelVectors.back() << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

        bool updateVector = !kernelVectors.back().vector;
        if (updateVector)
        {
//           success = false;

          bool updateApproximateVector = !kernelVectors.back().approximateVector;
          if (updateApproximateVector)
          {
            timeComponent = std::chrono::system_clock::now();
            kernelVectors.back().approximateVector = std::make_shared<sparse_vector<double>>();
            approximateAffineComplement.computeKernelVector(kernelVectors.back().column,
              *kernelVectors.back().approximateVector, kernelVectors.back().approximateRhs,
              query.epsilonCoefficient);
            result.numKernel++;
            result.timeKernel += elapsedTime(timeComponent);
            assert(!kernelVectors.back().approximateVector->empty());

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
            std::cout << "Computed approximate kernel vector is "
              << *kernelVectors.back().approximateVector << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

            if (elapsedTime(timeStarted) >= query.timeLimit)
            {
              result.dimension = AFFINEHULL_ERROR_TIMEOUT;
              return result;
            }

            kernelVectors.back().sparsity = kernelVectors.back().approximateVector->size();
            kernelVectors.back().norm = euclideanNorm(*kernelVectors.back().approximateVector);

            timeComponent = std::chrono::system_clock::now();
            auto redundancy = approximateRedundancyCheck.test(*kernelVectors.back().approximateVector,
              kernelVectors.back().norm);
            result.timeEquations += elapsedTime(timeComponent);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
            std::cout << "Redundancy check yields approximate violation " << redundancy.maxViolation
              << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

            if (elapsedTime(timeStarted) >= query.timeLimit)
            {
              result.dimension = AFFINEHULL_ERROR_TIMEOUT;
              return result;
            }

            kernelVectors.back().maxViolation = redundancy.maxViolation;
            kernelVectors.back().redundancyCoordinate = redundancy.maxCoordinate;

            std::push_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);
            continue;
          }
          else
          {
            timeComponent = std::chrono::system_clock::now();
            kernelVectors.back().vector = std::make_shared<sparse_vector<mpq_class>>();
            affineComplement.computeKernelVector(kernelVectors.back().column,
              *kernelVectors.back().vector, kernelVectors.back().rhs, 0);
            *kernelVectors.back().approximateVector = convertSparseVector<double>(
              *kernelVectors.back().vector);
            kernelVectors.back().approximateRhs = convertNumber<double>(
              kernelVectors.back().rhs);
            result.numKernel++;
            result.timeKernel += elapsedTime(timeComponent);
            assert(!kernelVectors.back().vector->empty());

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
            std::cout << "Computed exact kernel vector is "
              << *kernelVectors.back().vector << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

            if (elapsedTime(timeStarted) >= query.timeLimit)
            {
              result.dimension = AFFINEHULL_ERROR_TIMEOUT;
              return result;
            }

            kernelVectors.back().sparsity = kernelVectors.back().vector->size();
            kernelVectors.back().norm = euclideanNorm(*kernelVectors.back().vector);
          }
        }

        bool updateRedundancy = updateVector
          || kernelVectors.back().redundancyRank != redundancyCheck.rank();
        if (updateRedundancy)
        {
          timeComponent = std::chrono::system_clock::now();
          auto redundancy = redundancyCheck.test(*kernelVectors.back().vector,
            kernelVectors.back().norm);
          result.timeEquations += elapsedTime(timeComponent);

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
          std::cout << "Redundancy check yields exact violation " << redundancy.maxViolation << std::endl;
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

          if (elapsedTime(timeStarted) >= query.timeLimit)
          {
            result.dimension = AFFINEHULL_ERROR_TIMEOUT;
            return result;
          }

          kernelVectors.back().maxViolation = redundancy.maxViolation;
          kernelVectors.back().redundancyCoordinate = redundancy.maxCoordinate;
          kernelVectors.back().redundancyRank = redundancyCheck.rank();
        }

        if (!updateRedundancy)
          break;

        std::push_heap(kernelVectors.begin(), kernelVectors.end(), kernelVectorLess);
      }
      KernelVectorRational kernelVector = kernelVectors.back();
      kernelVectors.pop_back();

#if defined(IPO_DEBUG_AFFINE_HULL_PRINT)
      std::cout << "Selected kernel vector has maximum violation " << kernelVector.maxViolation << std::endl;
      assert(kernelVector.maxViolation > 0);
#endif /* IPO_DEBUG_AFFINE_HULL_PRINT */

      mpq_class kernelVectorValue = 0;

      int outcome = tryMaximization(polyhedron, affineComplement, &objective[0],
        kernelVectors, *kernelVector.vector, kernelVector.column, kernelVector.norm,
        kernelVectorValue, result.points, result.rays, query, 0.0,
        query.timeLimit - elapsedTime(timeStarted), result.timeOracles, result.timePointsRays,
        &approximateAffineComplement);
      if (outcome < INTERNAL_OKAY)
      {
        result.dimension = outcome;
        return result;
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

      outcome = tryMinimization(polyhedron, affineComplement, &objective[0], kernelVectors,
        *kernelVector.vector, kernelVector.column, kernelVector.norm, kernelVectorValue,
        result.points, result.rays, query, 0.0, query.timeLimit - elapsedTime(timeStarted),
        result.timeOracles, result.timePointsRays, &approximateAffineComplement);
      if (outcome < INTERNAL_OKAY)
      {
        result.dimension = outcome;
        return result;
      }
      else if (outcome == INTERNAL_FURTHER_ONE)
      {
        result.lowerBound++;
        continue;
      }

      // Add equation defined by kernelDirectionVector.

      assert(kernelVector.redundancyCoordinate < n);
      int error = addEquation(redundancyCheck, kernelVector.redundancyCoordinate, affineComplement,
        result.equations, kernelVector.vector, kernelVector.column, kernelVectorValue,
        query.epsilonFactorization, query.epsilonCoefficient, result.timeEquations,
        &approximateRedundancyCheck, &approximateAffineComplement);
      if (error == INTERNAL_NUMERICS)
      {
        result.dimension = AFFINEHULL_ERROR_NUMERICS;
        return result;
      }
      else
      {
        assert(!error);
      }
      --result.upperBound;
    }

    result.dimension = result.lowerBound;
    result.timeTotal = elapsedTime(timeStarted);
    return result;
  }

  template AffineHullResult<mpq_class> affineHull(
    std::shared_ptr<Polyhedron<mpq_class>> polyhedron,
    const AffineHullQuery& query,
    const std::vector<Constraint<mpq_class>>& knownEquations);

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
