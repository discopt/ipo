#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/data.hpp>

#include <memory>

namespace ipo
{
  /**
  * \brief Base class for all IPO oracles.
  */

  class Oracle
  {
  public:
    /**
     * \brief Constructs an oracle with given \p name.
     */

    IPO_EXPORT
    Oracle(const std::string& name);

    /**
     * \brief Returns the oracle's name.
     */

    IPO_EXPORT
    inline const std::string& name() const
    {
      return _name;
    }

    /**
     * \brief Returns the oracle's ambient space.
     */

    IPO_EXPORT
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

  protected:
    std::string _name;
    std::shared_ptr<Space> _space;
  };

  class OptimizationOracle: public Oracle
  {
  public:
    /**
     * \brief Structure for storing the query.
     **/

    struct Query
    {
#if defined(IPO_WITH_GMP)
      /// Are rational points/rays requested?
      bool rational;
#endif /* IPO_WITH_GMP */
      /// The caller will only use points having at least this value.
      double minObjectiveValue;
      /// Minimum number of solutions for early termination.
      std::size_t minNumSolutions;
      /// Maximum number of solutions to return.
      std::size_t maxNumSolutions;
      /// Time limit
      double timeLimit;

      /**
       * \brief Constructs the query structure.
       */

      IPO_EXPORT
      Query();

      /**
       * \brief Clears the query data.
       */

      IPO_EXPORT
      void reset();
    };

    /**
     * \brief Structure for storing the query result.
     **/

    struct Result
    {
      struct Point
      {
        Value objectiveValue;
        Vector vector;

        Point(const Vector& vector);

        Point(const Vector& vector, const Value& value);

        inline bool operator<(const Point& other) const
        {
          return objectiveValue > other.objectiveValue;
        }
      };

      /// Whether the time limit was reached.
      bool hitTimeLimit;
      /// Array with objective values and vectors of all points.
      std::vector<Point> points;
      /// Array with all rays.
      std::vector<Vector> rays;
      /// Upper bound on the optimal solution value.
      Value dualBound;

      /**
       * \brief Constructs the result structure.
       */

      IPO_EXPORT
      Result();

      /**
       * \brief Clears the result data.
       */

      IPO_EXPORT
      void reset();

      /**
       * \brief Returns \c true if the problem was infeasible.
       */

      IPO_EXPORT
      inline bool isInfeasible() const
      {
        return points.empty() && rays.empty();
      }

      /**
       * \brief Returns \c true if a point was found.
       */

      IPO_EXPORT
      inline bool isFeasible() const
      {
        return !points.empty();
      }

      /**
       * \brief Returns \c true if the problem was unbounded.
       */

      IPO_EXPORT
      inline bool isUnbounded() const
      {
        return !rays.empty();
      }

      /**
       * \brief Sorts the points by decreasing objective value.
       * 
       * Sorts the points by decreasing objective value.
       */

      IPO_EXPORT
      void sortPoints();

      /**
       * \brief Checks that points are sorted by decreasing objective value.
       * 
       * Checks that points are sorted by decreasing objective value.
       */

      IPO_EXPORT
      bool checkPointsSorted() const;
    };

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    OptimizationOracle(const std::string& name);

    /**
     * \brief Returns true iff the oracle is exact.
     *
     * Returns true iff the oracle is exact in the sense of being able to return solutions as
     * exact rational vectors.
     */

    IPO_EXPORT
    virtual bool isExact() const = 0;

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     * \param result Structure for returning the result.
     **/

    virtual void maximize(const double* objectiveVector, const Query& query, Result& result) = 0;

#if defined(IPO_WITH_GMP)

    /**
     * \brief Maximize a rational objective vector.
     *
     * Maximize a rational objective vector. The default implementation just converts the
     * objective vector to a floating-point vector and calls \ref maximize.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     * \param result Structure for returning the result.
     **/

    IPO_EXPORT
    virtual void maximize(const mpq_class* objectiveVector, const Query& query, Result& result);

#endif /* IPO_WITH_GMP */

  };

  /**
   * \brief Base class for a separation oracle.
   *
   * Base class for a separation oracle. When queried with a point, the oracle returns any number
   * (including none) of (less-than-or-equal) inequalities.
   */

  class SeparationOracle: public Oracle
  {
  public:
    /**
     * \brief Structure for storing the query.
     **/

    struct Query
    {
#if defined(IPO_WITH_GMP)
      /// Are rational inequalities in addition to the floating-point ones requested?
      bool rational;
#endif /* IPO_WITH_GMP */

      /// Maximum number of solutions to return.
      std::size_t maxNumInequalities;
      /// Time limit
      double timeLimit;

      /**
       * \brief Constructs the query structure.
       */

      IPO_EXPORT
      Query();

      /**
       * \brief Clears the query data.
       */

      IPO_EXPORT
      void reset();
    };

    /**
     * \brief Structure for storing the query result.
     **/

    struct Result
    {
      /// Whether the time limit was reached.
      bool hitTimeLimit;
      /// Array with constraints.
      std::vector<Constraint> constraints;

      /**
       * \brief Constructs the result structure.
       */

      IPO_EXPORT
      Result();

      /**
       * \brief Clears the result data.
       */

      IPO_EXPORT
      void reset();
    };

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    SeparationOracle(const std::string& name);

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \param result Structure for returning the result.
     **/

    IPO_EXPORT
    virtual void getInitial(const Query& query, Result& result);

    /**
     * \brief Separates a point/ray with floating-point coordinates.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param query Structure for query.
     * \param isPoint Whether a point shall be separated.
     * \param result Structure for returning the result.
     *
     * \returns \c true if and only if the point/ray was separated.
     **/

    virtual bool separate(const double* vector, bool isPoint, const Query& query,
      Result& result) = 0;

#if defined(IPO_WITH_GMP)

    /**
     * \brief Separates a point/ray with rational coordinates.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param query Structure for query.
     * \param isPoint Whether a point shall be separated.
     * \param result Structure for returning the result.
     *
     * \returns \c true if and only if the point/ray was separated.
     **/

    IPO_EXPORT
    virtual bool separate(const mpq_class* vector, bool isPoint, const Query& query,
      Result& result);

#endif /* IPO_WITH_GMP */

  };

} /* namespace ipo */
