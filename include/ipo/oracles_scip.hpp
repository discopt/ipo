#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#ifdef IPO_WITH_SCIP

#include <ipo/oracles_mip.hpp>
#include <ipo/constraint.hpp>
#include <unordered_map>

// This is necessary due to a bug in SCIP. Whether some functionality is in a macro or not depends
// on how you include it.
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
#endif

namespace ipo
{
  // Forward declarations.

  class SCIPOptimizationOracleDouble;
  class SCIPSeparationOracleDouble;

#if defined(IPO_WITH_GMP)
  typedef RationalMIPExtendedOptimizationOracle SCIPOptimizationOracleRational;

  typedef RationalMIPExtendedSeparationOracle SCIPSeparationOracleRational;
#endif /* IPO_WITH_GMP */

  /**
   * \brief A SCIP solver instance.
   *
   * A SCIP solver instance. Via \ref getOptimizationOracle and \ref getSeparationOracle one we
   * create corresponding oracles.
   */

  class SCIPSolver: public std::enable_shared_from_this<SCIPSolver>
  {
  public:

    /**
     * \brief Constructs a solver instance from a \c SCIP struct.
     *
     * Constructs a solver instance from a \c SCIP struct. Before using it, a shared_ptr has
     * to point to it.
     */

    IPO_EXPORT
    SCIPSolver(SCIP* scip);

    /**
     * \brief Constructs a solver instance from any file that SCIP can read. Before using it,
     * a shared_ptr has to point to it.
     */

    IPO_EXPORT
    SCIPSolver(const std::string& fileName);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    ~SCIPSolver();

    /**
     * \brief Returns the name of the instance.
     */

    IPO_EXPORT
    inline const std::string& name() const
    {
      return _name;
    }

    /**
     * \brief Returns the associated ambient space.
     */

    IPO_EXPORT
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    /**
     * \brief Returns the original (maximization) objective of the instance as a dense vector.
     */

    IPO_EXPORT
    inline const double* instanceObjective() const
    {
      return _instanceObjective;
    }

    /**
     * \brief Returns an optimization oracle for the polyhedron.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPOptimizationOracleDouble> getOptimizationOracleDouble()
    {
      return getOptimizationOracleDouble(
        std::make_shared<Constraint<double>>(alwaysSatisfiedConstraint<double>()));
    }

    /**
     * \brief Returns an optimization oracle for the requested \p face.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPOptimizationOracleDouble> getOptimizationOracleDouble(
      std::shared_ptr<Constraint<double>> face)
    {
      return std::make_shared<SCIPOptimizationOracleDouble>(shared_from_this(), face);
    }

#if defined(IPO_WITH_GMP)

    /**
     * \brief Returns an optimization oracle for the polyhedron.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPOptimizationOracleRational> getOptimizationOracleRational()
    {
      return getOptimizationOracleRational(
        std::make_shared<Constraint<rational>>(alwaysSatisfiedConstraint<rational>()));
    }

    /**
     * \brief Returns an optimization oracle for the requested \p face.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPOptimizationOracleRational> getOptimizationOracleRational(
      std::shared_ptr<Constraint<rational>> face)
    {
      auto approximateFace = std::make_shared<Constraint<double>>(constraintToDouble(*face));
      auto approximateOracle = getOptimizationOracleDouble(approximateFace);
      return std::make_shared<SCIPOptimizationOracleRational>(_extender, approximateOracle, face);
    }

#endif /* IPO_WITH_GMP*/

    /**
     * \brief Returns a separation oracle for the polyhedron.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracleDouble> getSeparationOracleDouble()
    {
      return getSeparationOracleDouble(
        std::make_shared<Constraint<double>>(alwaysSatisfiedConstraint<double>()));
    }

    /**
     * \brief Returns a separation oracle for the \p face.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracleDouble> getSeparationOracleDouble(std::shared_ptr<Constraint<double>> face)
    {
      return std::make_shared<SCIPSeparationOracleDouble>(shared_from_this(), face);
    }

#if defined(IPO_WITH_GMP)

    /**
     * \brief Returns a separation oracle for the polyhedron.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracleRational> getSeparationOracleRational()
    {
      return getSeparationOracleRational(
        std::make_shared<Constraint<rational>>(alwaysSatisfiedConstraint<rational>()));
    }

    /**
     * \brief Returns a separation oracle for the \p face.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracleRational> getSeparationOracleRational(
      std::shared_ptr<Constraint<rational>> face)
    {
      auto approximateFace = std::make_shared<Constraint<double>>(constraintToDouble(*face));
      auto approximateOracle = getSeparationOracleDouble(approximateFace);
      return std::make_shared<SCIPSeparationOracleRational>(approximateOracle, face);
    }

#endif /* IPO_WITH_GMP */

  protected:

    friend SCIPOptimizationOracleDouble;
    friend SCIPSeparationOracleDouble;

    /**
     * \brief Initializes the solver data.
     */

    void initialize();

    /**
     * \brief Makes the given \p face the current one.
     */

    IPO_EXPORT
    void setFace(std::shared_ptr<Constraint<double>> face);

  protected:
    /// Actual SCIP instance.
    SCIP* _scip;
    /// Maps coordinates to SCIP variables.
    std::vector<SCIP_VAR*> _variables;
    /// Maps SCIP variables to coordinates.
    std::unordered_map<SCIP_VAR*, std::size_t> _variablesToCoordinates;
    double* _instanceObjective;
    std::string _name;
    std::shared_ptr<Space> _space;
    std::shared_ptr<Constraint<double>> _currentFace;
    std::unordered_map<std::shared_ptr<Constraint<double>>, SCIP_CONS*> _faceConstraints;

#if defined(IPO_WITH_GMP)
    std::shared_ptr<RationalMIPExtender> _extender;
#endif /* IPO_WITH_GMP */
  };

  /**
  * \brief OptimizationOracle based on the SCIP solver.
  */

  class SCIPOptimizationOracleDouble: public OptimizationOracle<double>
  {
  public:

    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param face The face we are optimizing over.
     */

    IPO_EXPORT
    SCIPOptimizationOracleDouble(std::shared_ptr<SCIPSolver> solver,
      std::shared_ptr<Constraint<double>> face);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SCIPOptimizationOracleDouble();


    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<double>::Result maximizeDouble(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query) override
    {
      return maximize(objectiveVector, query);
    }

    /**
     * \brief Maximize an objective vector of type double.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<double>::Result maximize(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query);

  protected:
    friend SCIPSolver;

    /// The solver instance
    std::shared_ptr<SCIPSolver> _solver;
    /// The index of the face we are optimizing over.
    std::shared_ptr<Constraint<double>> _face;
  };

  /**
  * \brief SeparationOracle for the LP relaxation based on the SCIP solver.
  */

  class SCIPSeparationOracleDouble: public SeparationOracle<double>
  {
  public:
    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param faceIndex Indexes the face we are separating for.
     */

    IPO_EXPORT
    SCIPSeparationOracleDouble(std::shared_ptr<SCIPSolver> solver,
      std::shared_ptr<Constraint<double>> face);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SCIPSeparationOracleDouble();

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
    SeparationOracle<double>::Result getInitial(const SeparationOracle<double>::Query& query)
      override;

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

    IPO_EXPORT
    virtual SeparationOracle<double>::Result separateDouble(const double* vector, bool isPoint,
      const SeparationOracle<double>::Query& query)
    {
      return separate(vector, isPoint, query);
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
    virtual SeparationOracle<double>::Result separate(const double* vector, bool isPoint,
      const SeparationOracle<double>::Query& query);

  protected:
    friend SCIPSolver;

    /// The solver instance
    std::shared_ptr<SCIPSolver> _solver;
    /// The index of the face we are separating for.
    std::shared_ptr<Constraint<double>> _face;
  };

} /* namespace ipo */

#endif /* IPO_WITH_SCIP */
