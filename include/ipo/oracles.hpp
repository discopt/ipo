#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/sparse_vector.hpp>
#include <ipo/constraint.hpp>

#include <memory>

namespace ipo
{
  /**
  * \brief Base class for all IPO oracles.
  */

  template <typename R>
  class CommonOracle
  {
  public:
    /**
     * \brief Constructs an oracle with given \p name.
     * 
     * The parent constructor must set CommonOracle::_space properly.
     * 
     * \param name Name of the oracle.
     */

    IPO_EXPORT
    CommonOracle(const std::string& name)
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

  /// Outcome of a query to a RealOptimizationOracle or RationalOptimizationOracle.
  enum class OptimizationOutcome
  {
    /// A timeout occured.
    TIMEOUT = -1,
    FOUND_NOTHING = 0,
    UNBOUNDED = 1,
    INFEASIBLE = 2,
    FEASIBLE = 3
  };

  /**
   * \brief Structure for storing additional information of a query to an optimization oracle.
   **/

  template <typename R>
  struct CommonOptimizationQuery
  {
    /// Whether CommonOptimizationQuery::minObjectiveValue has a meaning.
    bool hasMinObjectiveValue;
    /**
      * \brief Strict lower bound on objective value of returned points.
      * 
      * Strict lower bound on objective value of returned points. The oracle must not return points
      * with a value less than or equal to this number.
      */
    R minObjectiveValue;
    /// Maximum number of solutions to return.
    std::size_t maxNumSolutions;
    /// Time limit for this oracle call.
    double timeLimit;

    /**
      * \brief Constructs the query structure.
      */

    IPO_EXPORT
    CommonOptimizationQuery()
      : hasMinObjectiveValue(false), minObjectiveValue(0),
      maxNumSolutions(std::numeric_limits<std::size_t>::max()),
      timeLimit(std::numeric_limits<double>::infinity())
    {

    }
  };

  /**
   * \brief Structure for storing the response of an optimization oracle.
   **/

  template <typename R>
  struct CommonOptimizationReponse
  {
    /**
     * \brief Structure for a returned point.
     */
    
    struct Point
    {
      /// Objective value of the point.
      R objectiveValue;
      /// Shared pointer to the actual point.
      std::shared_ptr<sparse_vector<R>> vector;

      /**
       * \brief Constructs the point.
       * 
       * Constructs the point. Point::objectiveValue has to be set properly by the caller.
       * 
       * \param vec Shared pointer to the point.
       */

      IPO_EXPORT
      Point(std::shared_ptr<sparse_vector<R>> vec)
        : objectiveValue(-std::numeric_limits<double>::signaling_NaN()), vector(vec)
      {

      }

      /**
       * \brief Constructs the point.
       * 
       * \param vec Shared pointer to the point.
       * \param value Objective value of the point.
       */

      IPO_EXPORT
      Point(std::shared_ptr<sparse_vector<R>> vec, const R& value)
        : objectiveValue(value), vector(vec)
      {

      }

      /**
       * \brief Comparison for points. Points with larger Point::objectiveValue are considered
       *        smaller.
       */

      IPO_EXPORT
      inline bool operator<(const Point& other) const
      {
        return objectiveValue > other.objectiveValue;
      }
    };

    /**
     * \brief Structure for a returned unbounded ray.
     */

    struct Ray
    {
    public:
      /// Shared pointer to the ray.
      std::shared_ptr<sparse_vector<R>> vector;

      /**
       * \brief Constructs the ray.
       * 
       * Constructs the ray. Ray::_norm is computed.
       * 
       * \param vec Shared pointer to the ray.
       */

      IPO_EXPORT
      Ray(std::shared_ptr<sparse_vector<R>> vec)
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
       */

      IPO_EXPORT
      Ray(std::shared_ptr<sparse_vector<R>> vec, double norm)
        : vector(vec), _norm(norm)
      {
        assert(norm == euclideanNorm(*vec));
      }

      /**
       * \brief Returns the Euclidean norm of the ray.
       */

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
    /// Best-known lower bound on the optimum.
    R primalBound;
    /// Whether we know a finite upper bound.
    bool hasDualBound;
    /// Upper bound on the optimum.
    R dualBound;
    /// Array with all points.
    std::vector<Point> points;
    /// Array with all rays.
    std::vector<Ray> rays;
    /// Whether the time limit was reached.
    bool hitTimeLimit;

    /**
      * \brief Constructs the result structure.
      */

    IPO_EXPORT
    CommonOptimizationReponse()
      : outcome(OptimizationOutcome::FOUND_NOTHING), primalBound(0), hasDualBound(false),
      dualBound(0), hitTimeLimit(false)
    {

    }

    /**
      * \brief Move-construcs the response.
      */

    IPO_EXPORT
    CommonOptimizationReponse(CommonOptimizationReponse&& other)
      : outcome(other.outcome), primalBound(std::move(other.primalBound)),
      hasDualBound(other.hasDualBound), dualBound(std::move(other.dualBound)),
      points(std::move(other.points)), rays(std::move(other.rays)), hitTimeLimit(other.hitTimeLimit)
    {

    }

    /**
      * \brief Move-assignment operator.
      */

    IPO_EXPORT
    inline CommonOptimizationReponse<R>& operator=(CommonOptimizationReponse<R>&& other)
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
   * \brief Base class for optimization oracles in double arithmetic.
   */

  class RealOptimizationOracle : public CommonOracle<double>
  {
  public:
    /// The query structure for the oracle.
    typedef CommonOptimizationQuery<double> Query;
    /// The response structure for the oracle.
    typedef CommonOptimizationReponse<double> Response;

    /**
     * \brief Constructs the oracle.
     * 
     * Constructs the oracle. The parent constructor must set CommonOracle::_space properly.
     * 
     * \param name Name of the oracle.
     */

    IPO_EXPORT
    RealOptimizationOracle(const std::string& name);

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximize(const double* objectiveVector, const Query& query) = 0;
  };

  /**
   * \brief Structure for storing the query to a separation oracle.
   **/

  struct CommonSeparationQuery
  {
    /// Maximum number of solutions to return.
    std::size_t maxNumInequalities;
    /// Time limit
    double timeLimit;

    /**
      * \brief Constructs the query structure.
      */

    IPO_EXPORT
    CommonSeparationQuery()
      : maxNumInequalities(std::numeric_limits<std::size_t>::max()),
      timeLimit(std::numeric_limits<double>::infinity())
    {

    }
  };

  /**
   *\brief Structure for storing the response of a separation oracle.
   **/

  template <typename R>
  struct CommonSeparationResponse
  {
    /// Array with constraints.
    std::vector<Constraint<R>> constraints;
    /// Whether the time limit was reached.
    bool hitTimeLimit;

    /**
      * \brief Constructs the result structure.
      */

    IPO_EXPORT
    CommonSeparationResponse()
      : hitTimeLimit(false)
    {

    }

    /**
      * \brief Move-constructs response.
      */

    IPO_EXPORT
    CommonSeparationResponse(CommonSeparationResponse&& other)
      : constraints(std::move(other.constraints)), hitTimeLimit(other.hitTimeLimit)
    {

    }

    /**
      * \brief Move-assignment operator.
      */

    IPO_EXPORT
    inline CommonSeparationResponse& operator=(CommonSeparationResponse&& other)
    {
      constraints = std::move(other.constraints);
      hitTimeLimit = other.hitTimeLimit;
      return *this;
    }

    /**
      * \brief Assignment operator.
      */

    IPO_EXPORT
    inline CommonSeparationResponse& operator=(const CommonSeparationResponse& other)
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
  std::ostream& operator<<(std::ostream& stream,
    const CommonOptimizationReponse<double>& response);

  /**
   * \brief Base class for separation oracles in floating-point arithmetic.
   *
   * When queried with a point, the oracle returns any number (including none) of less-than-or-equal
   * inequalities.
   */

  class RealSeparationOracle: public CommonOracle<double>
  {
  public:
    /// Query structure.
    typedef CommonSeparationQuery Query;
    /// Response structure.
    typedef CommonSeparationResponse<double> Response;

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    RealSeparationOracle(const std::string& name);

    /**
     * \brief Returns initially known inequalities.
     * 
     * The default implementation returns nothing.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
    virtual Response getInitial(const Query& query);

    /**
     * \brief Separates a point/ray given in floating-point arithmetic.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    IPO_EXPORT
    virtual Response separate(const double* vector, bool isPoint, const Query& query) = 0;

  };

  /**
   * \brief Prints the response of a real separation oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream,
    const CommonSeparationResponse<double>& response);

#if defined(IPO_WITH_GMP)
  /**
   * \brief Base class for optimization oracles in rational arithmetic.
   */

  class RationalOptimizationOracle : public CommonOracle<mpq_class>
  {
  public:
    /// The query structure for the oracle.
    typedef CommonOptimizationQuery<mpq_class> Query;
    /// The response structure for the oracle.
    typedef CommonOptimizationReponse<mpq_class> Response;

    /**
     * \brief Constructs the oracle.
     * 
     * Constructs the oracle. The parent constructor must set CommonOracle::_space properly.
     * 
     * \param name Name of the oracle.
     */

    IPO_EXPORT
    RationalOptimizationOracle(const std::string& name);

    /**
     * \brief Maximize a rational objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximize(const mpq_class* objectiveVector, const Query& query) = 0;

    /**
     * \brief Maximize a floating-point objective vector.
     * 
     * The default implementation creates a rational
     * objective from the floating-point one and calls RationalOptimizationOracle::maximize().
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximize(const double* objectiveVector, const Query& query);
  };

  /**
   * \brief Prints the response of a rational optimization oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream,
    const CommonOptimizationReponse<mpq_class>& response);

  /**
   * \brief Base class for separation oracles in rational arithmetic.
   *
   * When queried with a point, the oracle returns any number (including none) of less-than-or-equal
   * inequalities.
   */

  class RationalSeparationOracle: public CommonOracle<mpq_class>
  {
  public:
    /// Query structure.
    typedef CommonSeparationQuery Query;
    /// Response structure.
    typedef CommonSeparationResponse<mpq_class> Response;

    /**
     * \brief Constructs the oracle.
     */

    IPO_EXPORT
    RationalSeparationOracle(const std::string& name);

    /**
     * \brief Returns initially known inequalities.
     * 
     * The default implementation returns nothing.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
    virtual Response getInitial(const Query& query);

    /**
     * \brief Separates a point/ray given in rational arithmetic.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    IPO_EXPORT
    virtual Response separate(const mpq_class* vector, bool isPoint, const Query& query) = 0;

    /**
     * \brief Separates a point/ray given in floating-point arithmetic.
     *
     * The default implementation converts the point/ray to rational arithmetic and calls 
     * RationalSeparationOracle::separate.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    IPO_EXPORT
    virtual Response separate(const double* vector, bool isPoint, const Query& query);

  };

  /**
   * \brief Prints the response of a rational separation oracle.
   */

  IPO_EXPORT
  std::ostream& operator<<(std::ostream& stream,
    const CommonSeparationResponse<mpq_class>& response);

#endif /* IPO_WITH_GMP */

} /* namespace ipo */
