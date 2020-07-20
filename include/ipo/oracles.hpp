#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/sparse_vector.hpp>
#include <ipo/rational.hpp>
#include <ipo/constraint.hpp>

#include <memory>

namespace ipo
{
  /**
  * \brief Base class for all IPO oracles.
  */

  template <typename T>
  class Oracle
  {
  public:
    /**
     * \brief Constructs an oracle with given \p name.
     */

    IPO_EXPORT
    Oracle(const std::string& name)
      : _name(name), _space(nullptr)
    {

    }

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

  enum OptimizationOutcome
  {
    OPTIMIZATION_TIMEOUT = -1,
    OPTIMIZATION_FOUNDNOTHING= 0,
    OPTIMIZATION_UNBOUNDED = 1,
    OPTIMIZATION_INFEASIBLE = 2,
    OPTIMIZATION_FEASIBLE = 3
  };

  template <typename T>
  class OptimizationOracle: public Oracle<T>
  {
  public:
    /**
     * \brief Structure for storing the query.
     **/

    struct Query
    {
      bool hasMinObjectiveValue;
      /**
       * \brief Strict lower bound on objective value of returned points.
       * 
       * Strict lower bound on objective value of returned points. The oracle must not return points
       * with a value less than or equal to this number.
       */
      T minObjectiveValue;
      /// Maximum number of solutions to return.
      std::size_t maxNumSolutions;
      /// Time limit
      double timeLimit;

      /**
       * \brief Constructs the query structure.
       */

      IPO_EXPORT
      Query()
        : hasMinObjectiveValue(false), minObjectiveValue(0),
        maxNumSolutions(std::numeric_limits<std::size_t>::max()),
        timeLimit(std::numeric_limits<double>::infinity())
      {

      }

    };

    /**
     * \brief Structure for storing the query result.
     **/

    struct Response
    {
      struct Point
      {
        T objectiveValue;
        std::shared_ptr<sparse_vector<T>> vector;

        IPO_EXPORT
        Point(std::shared_ptr<sparse_vector<T>> vec)
          : objectiveValue(-std::numeric_limits<double>::signaling_NaN()), vector(vec)
        {

        }

        IPO_EXPORT
        Point(std::shared_ptr<sparse_vector<T>> vec, const T& value)
          : objectiveValue(value), vector(vec)
        {

        }

        IPO_EXPORT
        inline bool operator<(const Point& other) const
        {
          return objectiveValue > other.objectiveValue;
        }
      };

      struct Ray
      {
        std::shared_ptr<sparse_vector<T>> vector;

        IPO_EXPORT
        Ray(std::shared_ptr<sparse_vector<T>> vec)
          : vector(vec)
        {

        }

      };

      /// Outcome of query.
      OptimizationOutcome outcome;
      /// Lower bound on the optimum.
      T primalBound;
      /// Whether we know a finite upper bound.
      bool hasDualBound;
      /// Upper bound on the optimum.
      T dualBound;
      /// Array with objective values and vectors of all points.
      std::vector<Point> points;
      /// Array with all rays.
      std::vector<Ray> rays;
      /// Whether the time limit was reached.
      bool hitTimeLimit;

      /**
       * \brief Constructs the result structure.
       */

      IPO_EXPORT
      Response()
        : outcome(OPTIMIZATION_FOUNDNOTHING), primalBound(0), hasDualBound(false), dualBound(0),
        hitTimeLimit(false)
      {

      }

      /**
       * \brief Move-construcs the result structure.
       */

      IPO_EXPORT
      Response(Response&& other)
        : outcome(other.outcome), primalBound(std::move(other.primalBound)),
        hasDualBound(other.hasDualBound), dualBound(std::move(other.dualBound)),
        points(std::move(other.points)), rays(std::move(other.rays)), 
        hitTimeLimit(other.hitTimeLimit)
      {

      }

      /**
       * \brief Move-assignment operator.
       *
       * Move-assignment operator.
       */

      IPO_EXPORT
      inline Response& operator=(Response&& other)
      {
        outcome = other.outcome;
        primalBound = std::move(other.primalBound);
        hasDualBound = other.hasDualBound;
        dualBound = std::move(other.dualBound);
        points = std::move(other.points);
        rays = std::move(other.rays);
        hitTimeLimit = other.hitTimeLimit;
        return *this;
      }

      IPO_EXPORT
      inline bool wasSuccessful() const
      {
        return outcome == OPTIMIZATION_INFEASIBLE || outcome == OPTIMIZATION_UNBOUNDED
          || outcome == OPTIMIZATION_FEASIBLE;
      }

      /**
       * \brief Sorts the points by decreasing objective value.
       *
       * Sorts the points by decreasing objective value.
       */

      IPO_EXPORT
      void sortPoints()
      {
        std::sort(points.begin(), points.end());
      }

      /**
       * \brief Checks that points are sorted by decreasing objective value.
       *
       * Checks that points are sorted by decreasing objective value.
       */

      IPO_EXPORT
      bool checkPointsSorted() const
      {
        return std::is_sorted(points.begin(), points.end());
      }
    };

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    OptimizationOracle(const std::string& name)
      : Oracle<T>(name)
    {

    }

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual Response maximizeDouble(const double* objectiveVector, const Query& query)
    {
      std::vector<T> tempObjectiveVector(this->space()->dimension());
      for (std::size_t i = 0; i < this->space()->dimension(); ++i)
        tempObjectiveVector[i] = objectiveVector[i];
      return maximize(&tempObjectiveVector[0], query);
    }

    /**
     * \brief Maximize an objective vector of type \ref T.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual Response maximize(const T* objectiveVector, const Query& query) = 0;

  };

  /**
   * \brief Base class for a separation oracle.
   *
   * Base class for a separation oracle. When queried with a point, the oracle returns any number
   * (including none) of (less-than-or-equal) inequalities.
   */

  template <typename T>
  class SeparationOracle: public Oracle<T>
  {
  public:
    /**
     * \brief Structure for storing the query.
     **/

    struct Query
    {
      /// Maximum number of solutions to return.
      std::size_t maxNumInequalities;
      /// Time limit
      double timeLimit;

      /**
       * \brief Constructs the query structure.
       */

      IPO_EXPORT
      Query()
        : maxNumInequalities(std::numeric_limits<std::size_t>::max()),
        timeLimit(std::numeric_limits<double>::infinity())
      {

      }
    };

    /**
     * \brief Structure for storing the query result.
     **/

    struct Response
    {
      /// Array with constraints.
      std::vector<Constraint<T>> constraints;
      /// Whether the time limit was reached.
      bool hitTimeLimit;

      /**
       * \brief Constructs the result structure.
       */

      IPO_EXPORT
      Response()
        : hitTimeLimit(false)
      {

      }

      /**
       * \brief Move-constructs result structure.
       */

      IPO_EXPORT
      Response(Response&& other)
        : constraints(std::move(other.constraints)), hitTimeLimit(other.hitTimeLimit)
      {

      }

      /**
       * \brief Move-assignment operator.
       *
       * Move-assignment operator.
       */

      IPO_EXPORT
      inline Response& operator=(Response&& other)
      {
        constraints = std::move(other.constraints);
        hitTimeLimit = other.hitTimeLimit;
        return *this;
      }

      /**
       * \brief Assignment operator.
       *
       * Assignment operator.
       */

      IPO_EXPORT
      inline Response& operator=(const Response& other)
      {
        constraints = other.constraints;
        hitTimeLimit = other.hitTimeLimit;
        return *this;
      }

      /**
       * \brief Returns the number of generated constraints.
       */

      IPO_EXPORT
      inline std::size_t numConstraints() const
      {
        return constraints.size();
      }

    };

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    SeparationOracle(const std::string& name)
      : Oracle<T>(name)
    {

    }

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
    virtual Response getInitial(const Query& query)
    {
      return Response();
    }

    /**
     * \brief Separates a point/ray with floating-point coordinates.
     *
     * Separates a point/ray with floating-point coordinates. This default implementation converts
     * the double objective to one of type \ref T.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param query Structure for query.
     * \param isPoint Whether a point shall be separated.
     * \param result Structure for returning the result.
     *
     * \returns \c true if and only if the point/ray was separated.
     **/

    IPO_EXPORT
    virtual Response separateDouble(const double* vector, bool isPoint, const Query& query)
    {
      std::vector<T> tempVector(this->space()->dimension());
      for (std::size_t i = 0; i < this->space()->dimension(); ++i)
        tempVector[i] = vector[i];
      return separate(&tempVector[0], isPoint, query);
    }

    /**
     * \brief Separates a point/ray of the corresponding type.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param query Structure for query.
     * \param isPoint Whether a point shall be separated.
     * \param result Structure for returning the result.
     *
     * \returns \c true if and only if the point/ray was separated.
     **/

    IPO_EXPORT
    virtual Response separate(const T* vector, bool isPoint, const Query& query) = 0;

  };

  IPO_EXPORT
  inline std::ostream& operator<<(std::ostream& stream,
    OptimizationOracle<double>::Response& response)
  {
    switch(response.outcome)
    {
    case OPTIMIZATION_TIMEOUT:
      return stream << "{timeout error}";
    case OPTIMIZATION_FOUNDNOTHING:
      return stream << "{found nothing}";
    case OPTIMIZATION_UNBOUNDED:
      return stream << "{unbounded, " << response.rays.size() << " rays, " << response.points.size()
        << " points}";
    case OPTIMIZATION_INFEASIBLE:
      return stream << "{infeasible}";
    case OPTIMIZATION_FEASIBLE:
      assert(response.rays.empty());
      stream << "{feasible, " << response.rays.size() << " rays, " << response.points.size()
        << " points, " << double(response.primalBound) << " <= OPT";
      if (response.hasDualBound)
        stream << " <= " << double(response.dualBound);
      return stream << "}";
    default:
      return stream << "{unknown error}";
    }
  }

  IPO_EXPORT
  inline std::ostream& operator<<(std::ostream& stream,
    OptimizationOracle<ipo::rational>::Response& response)
  {
    switch(response.outcome)
    {
    case OPTIMIZATION_TIMEOUT:
      return stream << "{timeout error}";
    case OPTIMIZATION_FOUNDNOTHING:
      return stream << "{found nothing}";
    case OPTIMIZATION_UNBOUNDED:
      return stream << "{unbounded, " << response.rays.size() << " rays, " << response.points.size()
        << " points}";
    case OPTIMIZATION_INFEASIBLE:
      return stream << "{infeasible}";
    case OPTIMIZATION_FEASIBLE:
      assert(response.rays.empty());
      stream << "{feasible, " << response.rays.size() << " rays, " << response.points.size()
        << " points, " << double(response.primalBound) << " <= OPT";
      if (response.hasDualBound)
        stream << " <= " << double(response.dualBound);
      return stream << "}";
    default:
      return stream << "{unknown error}";
    }
  }

} /* namespace ipo */
