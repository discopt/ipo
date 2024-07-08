#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>
#include <ipo/sparse_vector.hpp>

#include <memory>

namespace ipo
{
  /**
   * \brief Projection map 
   *
   * \see \ref projectionEquations
   */
  
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

    Projection();

    /**
     * \brief Constructs an orthogonal projection map onto variables \p variableIndices of \p space.
     */

    Projection(std::shared_ptr<Space> space, const std::vector<std::size_t>& variableIndices);

    /**
     * \brief Constructs an orthogonal projection map onto variables of \p space based on the regular expresion \p regex.
     */

    Projection(std::shared_ptr<Space> space, const std::string& regex);

    /**
     * \brief Destructor.
     */

    ~Projection();

    /**
     * \brief Linear part of the map.
     */

    inline const sparse_vector<Number>& mapLinear(std::size_t variable) const
    {
      return _mapLinear[variable];
    }

    /**
     * \brief Constant part of the map.
     */

    inline const Number& mapConstant(std::size_t variable) const
    {
      return _mapConstant[variable];
    }

    /**
     * \brief Adds variable \p variableIndex from \p space by means of an orthogonal projection.
     **/

    void addVariable(const std::string& name, const sparse_vector<Number>& linear, const Number& constant);

    /**
     * \brief Adds variable \p variableIndex from \p space by means of an orthogonal projection.
     **/

    void addVariable(std::shared_ptr<Space> space, std::size_t variableIndex);

    /**
     * \brief Adds variables \p variableIndices from \p space by means of orthogonal projection.
     **/

    std::size_t addVariables(std::shared_ptr<Space> space, const std::vector<std::size_t>& variableIndices);

    /**
     * \brief Adds variables of \p space based on the regular expression \p regex by means of orthogonal projection.
     *
     * \returns Number of matches.
     **/

    std::size_t addVariables(std::shared_ptr<Space> space, const std::string& regex);

    /**
     * \brief Returns the target space.
     */

    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    std::shared_ptr<sparse_vector<Number>> projectPoint(std::shared_ptr<sparse_vector<Number>> point);

    std::shared_ptr<sparse_vector<Number>> projectRay(std::shared_ptr<sparse_vector<Number>> ray);

    /**
     * \brief Lifts an \p objective vector to \p liftedObjective.
     * 
     * \returns Lifted objective value of points that are projected to the origin.
     **/

    Number liftObjective(const Number* objective, std::vector<Number>& liftedObjective);

  protected:
    std::shared_ptr<Space> _space;                  /**< Defined space. */
    std::vector<sparse_vector<Number> > _mapLinear; /**< Vector containing the linear part of the projection map. */
    std::vector<Number> _mapConstant;               /**< Vector containing the absolute part of the projection map. */
  };

  /**
   * \brief Projects a set of \p equations via \p projection.
   *
   * Computes a system of equations in the image space of \p projection that defines the projection of the affine space
   * defined by \p equations.
   *
   * Suppose \p equations are \f$ Ax = b \f$ and \p projection is \f$ x \mapsto y = Tx + t \f$.
   * Hence, we want to find out the orthogonal projection of the affine space defined by
   *
   * Iy - Tx = t
   *      Ax = b
   *
   * This is achieved by carrying out pivots on columns corresponding to x-variables until the corresponding matrix is
   * upper triangular. After that, the rows with only nonzeros on y-variables define the projected space.
   **/

  template <typename Number>
  std::vector<Constraint<Number>> projectionEquations(std::shared_ptr<Projection<Number>> projection,
    const std::vector<Constraint<Number>>& equations);

  /**
   * \brief Projects a set of \p inequalities via \p projection.
   *
   * Computes a system of inequalities in the image space of \p projection that defines a relaxation of the projection
   * of the polyhedron defined by \p constraints.
   *
   * Suppose a constraint from \p constraints reads \f$ a^\intercal x \leq \beta \f$ and \p projection is
   * \f$ x \mapsto y = Tx + t \f$ then we obtain \f$ c^\intercal y \leq \delta \f$ if
   * \f$ c^\intercal T = a^\intercal \f$ and \f$ \delta - c^\intercal t = \beta \f$ holds.
   **/

  template <typename Number>
  std::vector<Constraint<Number>> projectionCompact(std::shared_ptr<Projection<Number>> projection,
    const std::vector<Constraint<Number>>& constraints, bool onlyEquations = false);

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

    ProjectionOptimizationOracle(std::shared_ptr<OptimizationOracle<NumberType>> sourceOracle,
      std::shared_ptr<Projection<NumberType>> projection,
      const std::string& name = "");

    /**
     * \brief Destructor.
     */

    virtual ~ProjectionOptimizationOracle();

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

    virtual OptimizationResponse<Number> maximize(const Number* objectiveVector,
      const OptimizationQuery<NumberType>& query);

  protected:

    std::shared_ptr<Projection<Number>> _projection;
    std::shared_ptr<OptimizationOracle<Number>> _sourceOracle; /**< Source oracle. */
  };

} /* namespace ipo */

