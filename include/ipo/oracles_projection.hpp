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

  protected:
    std::shared_ptr<Space> _space;                  /**< Defined space. */
    std::vector<sparse_vector<Number> > _mapLinear; /**< Vector containing the linear part of the projection map. */
    std::vector<Number> _mapConstant;               /**< Vector containing the absolute part of the projection map. */
  };

  std::vector<sparse_vector<double>> projectionEquations(const std::vector<sparse_vector<double>>& equations);

#if defined(IPO_WITH_GMP)
  std::vector<sparse_vector<mpq_class>> projectionEquations(const std::vector<sparse_vector<mpq_class>>& equations);
#endif /* IPO_WITH_GMP */

  class ProjectionRealOptimizationOracle: public OptimizationOracle<double>
  {
  public:

    /**
     * \brief Constructor for given projection and source oracle.
     *
     * Constructor for given \c projection and source \c oracle.
     *
     * \param name       Name of the new oracle.
     * \param projection Projection map.
     * \param oracle     Oracle in the source space.
     */

    IPO_EXPORT
    ProjectionRealOptimizationOracle(std::shared_ptr<OptimizationOracle<double>> sourceOracle, const std::string& name = "");

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~ProjectionRealOptimizationOracle();

    /**
     * \brief Adds variable index \p sourceVariableIndex to the oracle.
     **/

    IPO_EXPORT
    void addVariable(std::size_t sourceVariableIndex);

    /**
     * \brief Adds all those variables from the source oracle that match \p regex.
     *
     * \returns Number of matches.
     **/

    IPO_EXPORT
    std::size_t addVariables(const std::string& regex);

  protected:

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximize(const double* objectiveVector, const Query& query);

  protected:

    std::shared_ptr<OptimizationOracle<double>> _sourceOracle; /**< Source oracle. */
    std::vector<sparse_vector<double> > _projectionLinear; /**< Vector containing the linear part of the projection. */
    std::vector<double> _projectionConstant; /**< Vector containing the absolute part of the projection. */
  };

#if defined(IPO_WITH_GMP)

  class ProjectionRationalOptimizationOracle: public OptimizationOracle<mpq_class>
  {
  public:

    /**
     * \brief Constructor for given projection and source oracle.
     *
     * Constructor for given \c projection and source \c oracle.
     *
     * \param name       Name of the new oracle.
     * \param projection Projection map.
     * \param oracle     Oracle in the source space.
     */

    IPO_EXPORT
    ProjectionRationalOptimizationOracle(std::shared_ptr<OptimizationOracle<mpq_class>> sourceOracle,
      const std::string& name = "");

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~ProjectionRationalOptimizationOracle();

    /**
     * \brief Adds variable index \p sourceVariableIndex to the oracle.
     **/

    IPO_EXPORT
    void addVariable(std::size_t sourceVariableIndex);

    /**
     * \brief Adds all those variables from the source oracle that match \p regex.
     *
     * \returns Number of matches.
     **/

    IPO_EXPORT
    std::size_t addVariables(const std::string& regex);

  protected:

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual Response maximize(const mpq_class* objectiveVector, const Query& query);

  protected:

    std::shared_ptr<OptimizationOracle<mpq_class>> _sourceOracle; /**< Source oracle. */
    std::vector<sparse_vector<mpq_class> > _projectionLinear; /**< Vector containing the linear part of the projection. */
    std::vector<mpq_class> _projectionConstant; /**< Vector containing the absolute part of the projection. */
  };

#endif /* IPO_WITH_GMP */

} /* namespace ipo */

