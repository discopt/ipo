#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>
#include <ipo/sparse_vector.hpp>

#include <memory>

namespace ipo
{
  template <typename NumberType>
  class Projection
  {
  public:
    typedef NumberType Number;

    /**
     * \brief Constructs a projection map with target dimension 0.
     * 
     * Use \ref addVariable or \ref addVariables to build a non-trivial projection map.
     */
    IPO_EXPORT
    Projection();

    /**
     * \brief Constructs an orthogonal projection map onto variables \p variableIndices of \p space.
     */
    IPO_EXPORT
    Projection(std::shared_ptr<Space> space, const std::vector<std::size_t>& variableIndices);

    /**
     * \brief Constructs an orthogonal projection map onto variables of \p space based on the regular expresion \p regex.
     */
    IPO_EXPORT
    Projection(std::shared_ptr<Space> space, const std::string& regex);

    /**
     * \brief Destructor.
     */
    IPO_EXPORT
    ~Projection();

    /**
     * \brief Adds variable \p variableIndex from \p space by means of an orthogonal projection.
     **/

    IPO_EXPORT
    void addVariable(const std::string& name, const sparse_vector<Number>& linear, const Number& constant);

    /**
     * \brief Adds variable \p variableIndex from \p space by means of an orthogonal projection.
     **/

    IPO_EXPORT
    void addVariable(std::shared_ptr<Space> space, std::size_t variableIndex);

    /**
     * \brief Adds variables \p variableIndices from \p space by means of orthogonal projection.
     **/

    IPO_EXPORT
    std::size_t addVariables(std::shared_ptr<Space> space, const std::vector<std::size_t>& variableIndices);
    /**
     * \brief Adds variables of \p space based on the regular expression \p regex by means of orthogonal projection.
     *
     * \returns Number of matches.
     **/

    IPO_EXPORT
    std::size_t addVariables(std::shared_ptr<Space> space, const std::string& regex);

    /**
     * \brief Returns the target space.
     */

    IPO_EXPORT
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    IPO_EXPORT
    std::shared_ptr<sparse_vector<Number>> projectPoint(std::shared_ptr<sparse_vector<Number>> point);

    IPO_EXPORT
    std::shared_ptr<sparse_vector<Number>> projectRay(std::shared_ptr<sparse_vector<Number>> ray);

    /**
     * \brief Lifts an \p objective vector to \p liftedObjective.
     * 
     * \returns Lifted objective value of points that are projected to the origin.
     **/

    IPO_EXPORT
    Number liftObjective(const Number* objective, std::vector<Number>& liftedObjective);

  protected:
    std::shared_ptr<Space> _space;                  /**< Defined space. */
    std::vector<sparse_vector<Number> > _mapLinear; /**< Vector containing the linear part of the projection map. */
    std::vector<Number> _mapConstant;               /**< Vector containing the absolute part of the projection map. */
  };

  std::vector<sparse_vector<double>> projectionEquations(const std::vector<sparse_vector<double>>& equations);

#if defined(IPO_WITH_GMP)
  std::vector<sparse_vector<mpq_class>> projectionEquations(const std::vector<sparse_vector<mpq_class>>& equations);
#endif /* IPO_WITH_GMP */


  template <typename NumberType>
  class ProjectionOptimizationOracle: public OptimizationOracle<NumberType>
  {
  public:
    typedef NumberType Number;

    /**
     * \brief Constructor for given projection and source oracle.
     *
     * Constructor for given \c projection and source \c oracle.
     *
     * \param oracle     Oracle in the source space.
     * \param projection Projection map.
     * \param name       Name of the new oracle.
     */

    IPO_EXPORT
    ProjectionOptimizationOracle(std::shared_ptr<OptimizationOracle<NumberType>> sourceOracle,
      std::shared_ptr<Projection<NumberType>> projection,
      const std::string& name = "");

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~ProjectionOptimizationOracle();

    IPO_EXPORT
    inline std::shared_ptr<Projection<NumberType>> projection()
    {
      return _projection;
    }

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual OptimizationResponse<Number> maximize(const Number* objectiveVector,
      const OptimizationQuery<NumberType>& query);

  protected:

    std::shared_ptr<Projection<Number>> _projection;
    std::shared_ptr<OptimizationOracle<Number>> _sourceOracle; /**< Source oracle. */
  };

} /* namespace ipo */

