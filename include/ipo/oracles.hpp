#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/sparse_vector.hpp>
#include <ipo/constraint.hpp>

#include <memory>

namespace ipo
{

  /// Outcome of a query to an OptimizationOracle.
  enum class OptimizationOutcome
  {
    /// A timeout occured.
    TIMEOUT = -1,
    /// The oracle did not find a solution with objective value at least CommonOptimizationQuery::minPrimalBound.
    FOUND_NOTHING = 0,
    /// The oracle returned an unbounded ray.
    UNBOUNDED = 1,
    /// The oracle asserted that the polyhedron is empty.
    INFEASIBLE = 2,
    /// The oracle found a soluton of with objective value at least CommonOptimizationQuery::minPrimalBound.
    FEASIBLE = 3
  };

  /**
   * \brief Structure for storing additional information of a query to an optimization oracle.
   **/

  template <typename NumberType>
  struct OptimizationQuery
  {
  public:
    typedef NumberType Number;

    /// Maximum number of solutions to return.
    std::size_t maxNumSolutions;
    /// Time limit for this oracle call.
    double timeLimit;

  protected:
    /// Whether CommonOptimizationQuery::minPrimalBound has a meaning.
    bool _hasMinPrimalBound;
    /**
     * \brief Threshold for primal bound allowing the oracle to terminate.
     *
     * The oracle shall terminate if the primal bound is greater than or equal to this value.
     **/
    Number _minPrimalBound;
    /// Whether CommonOptimizationQuery::maxDualBound has a meaning.
    bool _hasMaxDualBound;
    /**
     * \brief Threshold for dual bound allowing the oracle to terminate.
     * 
     * The oracle shall terminate if the dual bound is less than or equal to this value.
     **/
    Number _maxDualBound;

  public:
    /**
     * \brief Constructs the query structure.
     **/

    IPO_EXPORT
    OptimizationQuery()
      : maxNumSolutions(std::numeric_limits<std::size_t>::max()),
      timeLimit(std::numeric_limits<double>::infinity()), _hasMinPrimalBound(false),
      _minPrimalBound(0), _hasMaxDualBound(false), _maxDualBound(0)
    {

    }

    IPO_EXPORT
    bool hasMinPrimalBound() const
    {
      return _hasMinPrimalBound;
    }

    IPO_EXPORT
    const Number& minPrimalBound() const
    {
      assert(_hasMinPrimalBound);
      return _minPrimalBound;
    }

    IPO_EXPORT
    void setMinPrimalBound(const Number& minPrimalBound)
    {
      _hasMinPrimalBound = true;
      _minPrimalBound = minPrimalBound;
    }

    IPO_EXPORT
    void removeMinPrimalBound()
    {
      _hasMinPrimalBound = false;
    }

    IPO_EXPORT
    bool hasMaxDualBound() const
    {
      return _hasMaxDualBound;
    }

    IPO_EXPORT
    const Number& maxDualBound() const
    {
      assert(_hasMaxDualBound);
      return _maxDualBound;
    }

    IPO_EXPORT
    void setMaxDualBound(const Number& maxDualBound)
    {
      _hasMaxDualBound = true;
      _maxDualBound = maxDualBound;
    }

    IPO_EXPORT
    void removeMaxDualBound()
    {
      _hasMaxDualBound = false;
    }
  };

  /**
   * \brief Structure for storing the response of an optimization oracle.
   **/

  template <typename NumberType>
  struct OptimizationResponse
  {
    typedef NumberType Number;

    /**
     * \brief Structure for a returned point.
     **/
    
    struct Point
    {
      /// Objective value of the point.
      Number objectiveValue;
      /// Shared pointer to the actual point.
      std::shared_ptr<sparse_vector<Number>> vector;

      /**
       * \brief Constructs the point.
       * 
       * Constructs the point. Point::objectiveValue has to be set properly by the caller.
       * 
       * \param vec Shared pointer to the point.
       **/

      IPO_EXPORT
      Point(std::shared_ptr<sparse_vector<Number>> vec)
        : objectiveValue(-std::numeric_limits<double>::signaling_NaN()), vector(vec)
      {

      }

      /**
       * \brief Constructs the point.
       * 
       * \param vec Shared pointer to the point.
       * \param value Objective value of the point.
       **/

      IPO_EXPORT
      Point(std::shared_ptr<sparse_vector<Number>> vec, const Number& value)
        : objectiveValue(value), vector(vec)
      {

      }

      /**
       * \brief Comparison for points. Points with larger Point::objectiveValue are considered
       *        smaller.
       **/

      IPO_EXPORT
      inline bool operator<(const Point& other) const
      {
        return objectiveValue > other.objectiveValue;
      }
    };

    /**
     * \brief Structure for a returned unbounded ray.
     **/

    struct Ray
    {
    public:
      /// Shared pointer to the ray.
      std::shared_ptr<sparse_vector<Number>> vector;

      /**
       * \brief Constructs the ray.
       * 
       * Constructs the ray. Ray::_norm is computed.
       * 
       * \param vec Shared pointer to the ray.
       **/

      IPO_EXPORT
      Ray(std::shared_ptr<sparse_vector<Number>> vec)
        : vector(vec), _norm(euclideanNorm(*vec))
      {

      }

      /**
       * \brief Constructs the ray.
       * 
       * Constructs the ray.
       * 
       * \param vec Shared pointer to the ray.
       * \param norm Euclidean norm of the ray.
       **/

      IPO_EXPORT
      Ray(std::shared_ptr<sparse_vector<Number>> vec, double norm)
        : vector(vec), _norm(norm)
      {
        assert(norm == euclideanNorm(*vec));
      }

      /**
       * \brief Returns the Euclidean norm of the ray.
       **/

      IPO_EXPORT
      inline double norm() const
      {
        return _norm;
      }

    private:
      /// Euclidean norm of the ray.
      double _norm;
    };

    /// Outcome of query.
    OptimizationOutcome outcome;
    /// Whether we know a finite lower bound.
    bool hasPrimalBound;
    /// Best-known lower bound on the optimum.
    Number primalBound;
    /// Whether we know a finite upper bound.
    bool hasDualBound;
    /// Upper bound on the optimum.
    Number dualBound;
    /// Array with all points.
    std::vector<Point> points;
    /// Array with all rays.
    std::vector<Ray> rays;
    /// Whether the time limit was reached.
    bool hitTimeLimit;

    /**
     * \brief Constructs the result structure.
     **/

    IPO_EXPORT
    OptimizationResponse()
      : outcome(OptimizationOutcome::FOUND_NOTHING), hasPrimalBound(false), primalBound(0), hasDualBound(false),
      dualBound(0), hitTimeLimit(false)
    {

    }

    /**
     * \brief Move-construcs the response.
     **/

    IPO_EXPORT
    OptimizationResponse(OptimizationResponse&& other)
      : outcome(other.outcome), hasPrimalBound(other.hasPrimalBound), primalBound(std::move(other.primalBound)),
      hasDualBound(other.hasDualBound), dualBound(std::move(other.dualBound)), points(std::move(other.points)),
      rays(std::move(other.rays)), hitTimeLimit(other.hitTimeLimit)
    {

    }

    /**
     * \brief Move-assignment operator.
     **/

    IPO_EXPORT
    inline OptimizationResponse<Number>& operator=(OptimizationResponse<Number>&& other)
    {
      outcome = other.outcome;
      hasPrimalBound = other.hasPrimalBound;
      primalBound = std::move(other.primalBound);
      hasDualBound = other.hasDualBound;
      dualBound = std::move(other.dualBound);
      points = std::move(other.points);
      rays = std::move(other.rays);
      hitTimeLimit = other.hitTimeLimit;
      return *this;
    }

    /**
      * \brief Copy-assignment operator.
      */

    IPO_EXPORT
    inline OptimizationResponse<Number>& operator=(const OptimizationResponse<Number>& other)
    {
      outcome = other.outcome;
      hasPrimalBound = other.hasPrimalBound;
      primalBound = other.primalBound;
      hasDualBound = other.hasDualBound;
      dualBound = other.dualBound;
      points = other.points;
      rays = other.rays;
      hitTimeLimit = other.hitTimeLimit;
      return *this;
    }

    /**
     * \brief Whether the query was successful, i.e., no timeout or other error occured.
     */

    IPO_EXPORT
    inline bool wasSuccessful() const
    {
      return outcome == OptimizationOutcome::INFEASIBLE
        || outcome == OptimizationOutcome::UNBOUNDED
        || outcome == OptimizationOutcome::FEASIBLE;
    }
  };

  /**
   * \brief Base class for all IPO oracles.
   */

  template <typename NumberType>
  class Oracle
  {
  public:
    typedef NumberType Number;

    /**
     * \brief Constructs an oracle with given \p name.
     *
     * The parent constructor must set CommonOracle::_space properly.
     *
     * \param name Name of the oracle.
     */

    IPO_EXPORT
    Oracle(const std::string& name = "")
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
    /// Oracle's name.
    std::string _name;

    /// Oracle's ambient space.
    std::shared_ptr<Space> _space;
  };

  /**
   * \brief Base class for optimization oracles in double arithmetic.
   */

  template <typename NumberType>
  class OptimizationOracle : public Oracle<NumberType>
  {
  public:
    typedef NumberType Number;

    /// The query structure for the oracle.
    typedef OptimizationQuery<NumberType> Query;
    /// The response structure for the oracle.
    typedef OptimizationResponse<NumberType> Response;

    /**
     * \brief Constructs the oracle.
     * 
     * Constructs the oracle. The parent constructor must set CommonOracle::_space properly.
     * 
     * \param name Name of the oracle.
     */

    IPO_EXPORT
    OptimizationOracle(const std::string& name = "");

    /**
     * \brief Maximize an objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximize(const Number* objectiveVector, const Query& query) = 0;

    /**
     * \brief Maximize a double objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximizeDouble(const double* objectiveVector, const Query& query);
  };

  /**
   * \brief Structure for storing the query to a separation oracle.
   **/

  struct SeparationQuery
  {
    /// Maximum number of solutions to return.
    std::size_t maxNumInequalities;
    /// Time limit
    double timeLimit;

    /**
      * \brief Constructs the query structure.
      */

    IPO_EXPORT
    SeparationQuery()
      : maxNumInequalities(std::numeric_limits<std::size_t>::max()),
      timeLimit(std::numeric_limits<double>::infinity())
    {

    }
  };

  /**
   * \brief Structure for storing the response of a separation oracle.
   **/

  template <typename NumberType>
  struct SeparationResponse
  {
    typedef NumberType Number;

    /// Array with constraints.
    std::vector<Constraint<Number>> constraints;
    /// Whether the time limit was reached.
    bool hitTimeLimit;

    /**
     * \brief Constructs the result structure.
     **/

    IPO_EXPORT
    SeparationResponse()
      : hitTimeLimit(false)
    {

    }

    /**
     * \brief Move-constructs response.
     **/

    IPO_EXPORT
    SeparationResponse(SeparationResponse&& other)
      : constraints(std::move(other.constraints)), hitTimeLimit(other.hitTimeLimit)
    {

    }

    /**
     * \brief Move-assignment operator.
     */

    IPO_EXPORT
    inline SeparationResponse& operator=(SeparationResponse&& other)
    {
      constraints = std::move(other.constraints);
      hitTimeLimit = other.hitTimeLimit;
      return *this;
    }

    /**
     * \brief Assignment operator.
     */

    IPO_EXPORT
    inline SeparationResponse& operator=(const SeparationResponse& other)
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
   * \brief Prints the response of a real optimization oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const OptimizationResponse<double>& response);

  /**
   * \brief Base class for separation oracles.
   *
   * When queried with a point, the oracle returns any number (including none) of less-than-or-equal
   * inequalities.
   */

  template <typename NumberType>
  class SeparationOracle: public Oracle<NumberType>
  {
  public:
    typedef NumberType Number;

    /// Query structure.
    typedef SeparationQuery Query;
    /// Response structure.
    typedef SeparationResponse<Number> Response;

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    SeparationOracle(const std::string& name);

    /**
     * \brief Returns initially known inequalities.
     * 
     * The default implementation returns nothing.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
    virtual Response getInitial(const Query& query = Query());

    /**
     * \brief Separates a point/ray.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    IPO_EXPORT
    virtual Response separate(const Number* vector, bool isPoint, const Query& query = Query()) = 0;

    /**
     * \brief Separates a point/ray given in floating-point arithmetic.
     *
     * The default implementation converts the point/ray to a suitable arithmetic and calls 
     * \ref SeparationOracle::separate.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    IPO_EXPORT
    virtual Response separateDouble(const double* vector, bool isPoint, const Query& query);
  };

  /**
   * \brief Prints the response of a real separation oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const SeparationResponse<double>& response);

#if defined(IPO_WITH_GMP)

  /**
   * \brief Prints the response of a rational optimization oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const OptimizationResponse<mpq_class>& response);

  /**
   * \brief Prints the response of a rational separation oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream, const SeparationResponse<mpq_class>& response);

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
