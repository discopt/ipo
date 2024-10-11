#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#if defined(IPO_DOUBLE_MIP_GUROBI) || defined(IPO_RATIONAL_MIP_GUROBI)

#include <ipo/oracles.hpp>
#include <ipo/mip.hpp>
#include <ipo/constraint.hpp>
#include <unordered_map>
#include <gurobi_c.h>

namespace ipo
{
  // Forward declarations.
  
  template <typename NumberType>
  class GurobiOptimizationOracle;

  template <typename NumberType>
  class GurobiSeparationOracle;

  /**
   * \brief A Gurobi solver instance.
   *
   * A Gurobi solver instance. Via \ref getOptimizationOracle and \ref getSeparationOracle one we
   * create corresponding oracles.
   */

  class GurobiSolver: public std::enable_shared_from_this<GurobiSolver>
  {
  public:

    /**
     * \brief Constructs a solver instance from a \c Gurobi pointer.
     *
     * Before using it, a shared_ptr has to point to it.
     */

    GurobiSolver(GRBmodel*&& model);

    /**
     * \brief Constructs a solver instance from any file that Gurobi can read.
     *
     * Before using it, a shared_ptr has to point to it.
     */

    GurobiSolver(const std::string& fileName);

    /**
     * \brief Destructor.
     */

    ~GurobiSolver();

    /**
     * \brief Returns the name of the instance.
     */

    inline const std::string& name() const
    {
      return _name;
    }

    /**
     * \brief Returns the associated ambient space.
     */

    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    /**
     * \brief Returns the original (maximization) objective of the instance as a dense vector.
     */

    inline const double* instanceObjective() const
    {
      return _instanceObjective;
    }

    /**
     * \brief Returns the original (maximization) objective offset of the instance.
     */

    inline double instanceObjectiveOffset() const
    {
      return _instanceObjective[_space->dimension()];
    }

    /**
     * \brief Returns an optimization oracle for the polyhedron.
     */

    template <typename NumberType>
    inline std::shared_ptr<GurobiOptimizationOracle<NumberType>> getOptimizationOracle()
    {
      return getOptimizationOracle<NumberType>(alwaysSatisfiedConstraint<NumberType>());
    }

    /**
     * \brief Returns an optimization oracle for the requested \p face.
     */

    template <typename NumberType>
    std::shared_ptr<GurobiOptimizationOracle<NumberType>> getOptimizationOracle(
      const Constraint<NumberType>& face);

    /**
     * \brief Returns a separation oracle for the polyhedron.
     */

    template <typename NumberType>
    inline std::shared_ptr<GurobiSeparationOracle<NumberType>> getSeparationOracle()
    {
      return getSeparationOracle<NumberType>(alwaysSatisfiedConstraint<NumberType>());
    }

    /**
     * \brief Returns a separation oracle for the \p face.
     */

    template <typename NumberType>
    std::shared_ptr<GurobiSeparationOracle<NumberType>> getSeparationOracle(
      const Constraint<NumberType>& face);

  protected:

#if defined(IPO_DOUBLE_MIP_GUROBI)

    friend GurobiOptimizationOracle<double>;
    friend GurobiSeparationOracle<double>;

#endif /* IPO_DOUBLE_MIP_GUROBI */
    
#if defined(IPO_RATIONAL_MIP_GUROBI)

    friend GurobiOptimizationOracle<rational>;
    friend GurobiSeparationOracle<rational>;

#endif /* IPO_RATIONAL_MIP_GUROBI */

    /**
     * \brief Initializes the solver data.
     */

    IPO_NO_EXPORT
    void initialize();

    /**
     * \brief Adds the given \p face to the list of known faces.
     */

    void addFace(Constraint<double>* face);

    /**
     * \brief Removes the given \p face from the list of known faces.
     */

    void deleteFace(Constraint<double>* face);
    /**
     * \brief Makes the given \p face the current one.
     */

    void selectFace(Constraint<double>* face);


  protected:
    struct BoundLimits
    {
      double minPrimalBound;
      double maxDualBound;
    };

    /// Gurobi environment (might be \c NULL).
    GRBenv* _env;
    /// Actual Gurobi instance.
    GRBmodel* _model;
    /// Maps Gurobi variables to coordinates.
    std::unordered_map<std::size_t, std::size_t> _variablesToCoordinates;
    double* _instanceObjective;
    std::string _name;
    std::shared_ptr<Space> _space;
    Constraint<double>* _currentFace;
    BoundLimits _boundLimits;

#if defined(IPO_RATIONAL_MIP_GUROBI)
    RationalMIPExtender* _extender;
#endif /* IPO_RATIONAL_MIP_GUROBI */
  };

#if defined(IPO_DOUBLE_MIP_GUROBI)

  /**
  * \brief OptimizationOracle based on the Gurobi solver.
  */

  template <>
  class GurobiOptimizationOracle<double>: public OptimizationOracle<double>
  {
  public:

    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param face The face we are optimizing over.
     */

    GurobiOptimizationOracle(std::shared_ptr<GurobiSolver> solver,
      const Constraint<double>& face);

    /**
     * \brief Destructor.
     */

    virtual ~GurobiOptimizationOracle();

    /**
     * \brief Maximize an objective vector of type double.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    virtual OptimizationOracle<double>::Response maximize(const double* objectiveVector,
      const OptimizationOracle<double>::Query& query) override;

  protected:
    friend GurobiSolver;

    /// The solver instance
    std::shared_ptr<GurobiSolver> _solver;
    /// The index of the face we are optimizing over.
    Constraint<double> _face;
  };

#endif /* IPO_DOUBLE_MIP_GUROBI */

#if defined(IPO_RATIONAL_MIP_GUROBI)

  template <>
  class GurobiOptimizationOracle<rational>: public RationalMIPExtendedOptimizationOracle
  {
  public:
    GurobiOptimizationOracle(RationalMIPExtender* extender, std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const Constraint<rational>& face)
      : RationalMIPExtendedOptimizationOracle(extender, approximateOracle, face)
    {
      
    }
  };

#endif /* IPO_RATIONAL_MIP_GUROBI */

  /**
  * \brief SeparationOracle for the LP relaxation based on the Gurobi solver.
  */

  template <typename NumberType>
  class GurobiSeparationOracle: public SeparationOracle<NumberType>
  {
  public:
    /**
     * \brief Constructs oracle using the \p solver.
     *
     * \param solver The solver instance that is used to answer the queries.
     * \param faceIndex Indexes the face we are separating for.
     */

    GurobiSeparationOracle(std::shared_ptr<GurobiSolver> solver,
      const Constraint<NumberType>& face);

    /**
     * \brief Destructor.
     */

    virtual ~GurobiSeparationOracle();

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

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

    virtual SeparationResponse<NumberType> separate(const NumberType* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery()) override;

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

    virtual SeparationResponse<NumberType> separateDouble(const double* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery()) override;

  protected:
    friend GurobiSolver;

    /// The solver instance
    std::shared_ptr<GurobiSolver> _solver;
    /// The index of the face we are separating for.
    Constraint<NumberType> _face;
    Constraint<double> _approximateFace;
  };

} /* namespace ipo */

#endif /* IPO_DOUBLE_MIP_GUROBI || IPO_RATIONAL_MIP_GUROBI */
