#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#ifdef IPO_WITH_SCIP

#include <ipo/oracles.hpp>
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

#ifdef IPO_WITH_SOPLEX
#include <ipo/make_rational_soplex.hpp>
#endif

namespace ipo
{
  // Forward declarations.

  class SCIPOptimizationOracle;
  class SCIPSeparationOracle;

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
    inline std::shared_ptr<SCIPOptimizationOracle> getOptimizationOracle()
    {
      return getOptimizationOracle(alwaysSatisfiedConstraint());
    }

    /**
     * \brief Returns an optimization oracle for the requested \p face.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPOptimizationOracle> getOptimizationOracle(Constraint constraint)
    {
      return std::make_shared<SCIPOptimizationOracle>(shared_from_this(), constraint);
    }

    /**
     * \brief Returns a separation oracle for the polyhedron.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracle> getSeparationOracle()
    {
      return getSeparationOracle(alwaysSatisfiedConstraint());
    }

    /**
     * \brief Returns a separation oracle for the requested \p face.
     */

    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracle> getSeparationOracle(Constraint constraint)
    {
      return std::make_shared<SCIPSeparationOracle>(shared_from_this(), constraint);
    }

  protected:

    friend SCIPOptimizationOracle;
    friend SCIPSeparationOracle;

    /**
     * \brief Initializes the solver data.
     */

    void initialize();

    /**
     * \brief Makes the given \p face the current one.
     */

    IPO_EXPORT
    void setFace(const Constraint& constraint);

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    void makeRational(OptimizationOracle::Result& result, const mpq_class* objectiveVector);

    void makeRational(OptimizationOracle::Result& result, const double* objectiveVector);
#endif
    
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
    Constraint _currentFace;
    std::unordered_map<Constraint, SCIP_CONS*, HashConstraint> _faceConstraints;
#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    MakeRationalSolver* _makeRationalSolver;
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
  };

  /**
  * \brief OptimizationOracle based on the SCIP solver.
  */

  class SCIPOptimizationOracle: public OptimizationOracle
  {
  public:

    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param faceIndex The face we are optimizing over.
     */

    IPO_EXPORT
    SCIPOptimizationOracle(std::shared_ptr<SCIPSolver> solver, Constraint face);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SCIPOptimizationOracle();

    /**
     * \brief Returns true iff the oracle is exact.
     *
     * Returns true iff the oracle is exact, i.e., upon request it can return solutions as exact
     * rational vectors.
     */

    IPO_EXPORT
    bool isExact() const override;

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     * \param result Structure for returning the result.
     **/

    IPO_EXPORT
    void maximize(const double* objectiveVector, const Query& query, Result& result) override;

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
    void maximize(const mpq_class* objectiveVector, const Query& query, Result& result) override;

#endif /* IPO_WITH_GMP */

  protected:
    friend SCIPSolver;

    /// The solver instance
    std::shared_ptr<SCIPSolver> _solver;
    /// The index of the face we are optimizing over.
    Constraint _face;
  };

  /**
  * \brief SeparationOracle for the LP relaxation based on the SCIP solver.
  */

  class SCIPSeparationOracle: public SeparationOracle
  {
  public:
    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param faceIndex Indexes the face we are separating for.
     */

    IPO_EXPORT
    SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver, Constraint face);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SCIPSeparationOracle();

    /**
     * \brief Separates a point/ray with floating-point coordinates.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param result Structure for returning the result.
     * \param timeLimit Time limit for this call (in seconds).
     *
     * \returns \c true if and only if the point/ray was separated.
     */

    IPO_EXPORT
    virtual bool separate(const double* vector, bool isPoint, const SeparationOracle::Query& query,
      SeparationOracle::Result& result);

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
    virtual bool separate(const mpq_class* vector, bool isPoint,
      const SeparationOracle::Query& query, SeparationOracle::Result& result);

#endif /* IPO_WITH_GMP */

  protected:
    friend SCIPSolver;

    /// The solver instance
    std::shared_ptr<SCIPSolver> _solver;
    /// The index of the face we are separating for.
    Constraint _face;
  };

} /* namespace ipo */

#endif /* IPO_WITH_SCIP */
