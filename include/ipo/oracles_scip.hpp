#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#ifdef IPO_WITH_SCIP

#include <ipo/oracles.hpp>
#include <ipo/mip.hpp>
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
  
  template <typename NumberType>
  class SCIPOptimizationOracle;

  template <typename NumberType>
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
    SCIPSolver(SCIP*&& scip);

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

    template <typename NumberType>
    IPO_EXPORT
    inline std::shared_ptr<SCIPOptimizationOracle<NumberType>> getOptimizationOracle()
    {
      return getOptimizationOracle<NumberType>(alwaysSatisfiedConstraint<NumberType>());
    }

    /**
     * \brief Returns an optimization oracle for the requested \p face.
     */

    template <typename NumberType>
    IPO_EXPORT
    std::shared_ptr<SCIPOptimizationOracle<NumberType>> getOptimizationOracle(
      const Constraint<NumberType>& face);

    /**
     * \brief Returns a separation oracle for the polyhedron.
     */

    template <typename NumberType>
    IPO_EXPORT
    inline std::shared_ptr<SCIPSeparationOracle<NumberType>> getSeparationOracle()
    {
      return getSeparationOracle<NumberType>(alwaysSatisfiedConstraint<NumberType>());
    }

    /**
     * \brief Returns a separation oracle for the \p face.
     */

    template <typename NumberType>
    IPO_EXPORT
    std::shared_ptr<SCIPSeparationOracle<NumberType>> getSeparationOracle(
      const Constraint<NumberType>& face);

  protected:

    friend SCIPOptimizationOracle<double>;
    friend SCIPSeparationOracle<double>;

#if defined(IPO_WITH_GMP)
    friend SCIPOptimizationOracle<mpq_class>;
    friend SCIPSeparationOracle<mpq_class>;
#endif /* IPO_WITH_GMP */

    /**
     * \brief Initializes the solver data.
     */

    void initialize();

    /**
     * \brief Adds the given \p face to the list of known faces.
     */

    IPO_EXPORT
    void addFace(Constraint<double>* face);

    /**
     * \brief Removes the given \p face from the list of known faces.
     */

    IPO_EXPORT
    void deleteFace(Constraint<double>* face);
    /**
     * \brief Makes the given \p face the current one.
     */

    IPO_EXPORT
    void selectFace(Constraint<double>* face);


  protected:
    struct BoundLimits
    {
      double minPrimalBound;
      double maxDualBound;
    };
    
    /// Actual SCIP instance.
    SCIP* _scip;
    /// Maps coordinates to SCIP variables.
    std::vector<SCIP_VAR*> _variables;
    /// Maps SCIP variables to coordinates.
    std::unordered_map<SCIP_VAR*, std::size_t> _variablesToCoordinates;
    double* _instanceObjective;
    std::string _name;
    std::shared_ptr<Space> _space;
    Constraint<double>* _currentFace;
    std::unordered_map<Constraint<double>*, SCIP_CONS*> _faceConstraints;
    BoundLimits _boundLimits;

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)
    RationalMIPExtender* _extender;
#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
  };

  /**
  * \brief OptimizationOracle based on the SCIP solver.
  */

  template <>
  class SCIPOptimizationOracle<double>: public OptimizationOracle<double>
  {
  public:

    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param face The face we are optimizing over.
     */

    IPO_EXPORT
    SCIPOptimizationOracle(std::shared_ptr<SCIPSolver> solver,
      const Constraint<double>& face);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SCIPOptimizationOracle();

    /**
     * \brief Maximize an objective vector of type double.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<double>::Response maximize(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query) override;

  protected:
    friend SCIPSolver;

    /// The solver instance
    std::shared_ptr<SCIPSolver> _solver;
    /// The index of the face we are optimizing over.
    Constraint<double> _face;
  };

  template <>
  class SCIPOptimizationOracle<mpq_class>: public RationalMIPExtendedOptimizationOracle
  {
  public:
    IPO_EXPORT
    SCIPOptimizationOracle(RationalMIPExtender* extender, std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const Constraint<mpq_class>& face)
      : RationalMIPExtendedOptimizationOracle(extender, approximateOracle, face)
    {
      
    }
  };

  /**
  * \brief SeparationOracle for the LP relaxation based on the SCIP solver.
  */

  template <typename NumberType>
  class SCIPSeparationOracle: public SeparationOracle<NumberType>
  {
  public:
    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param faceIndex Indexes the face we are separating for.
     */

    IPO_EXPORT
    SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver,
      const Constraint<NumberType>& face);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SCIPSeparationOracle();

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
    SeparationResponse<NumberType> getInitial(const SeparationQuery& query = SeparationQuery()) override;

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
    virtual SeparationResponse<NumberType> separate(const NumberType* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery());

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
    virtual SeparationResponse<NumberType> separateDouble(const double* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery());

  protected:
    friend SCIPSolver;

    /// The solver instance
    std::shared_ptr<SCIPSolver> _solver;
    /// The index of the face we are separating for.
    Constraint<NumberType> _face;
    Constraint<double> _approximateFace;
  };

} /* namespace ipo */

#endif /* IPO_WITH_SCIP */
