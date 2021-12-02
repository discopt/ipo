#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>
#include <ipo/sparse_vector.hpp>

#include <memory>

namespace ipo
{

  class ProjectionRealOptimizationOracle: public RealOptimizationOracle
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

    ProjectionRealOptimizationOracle(std::shared_ptr<RealOptimizationOracle> sourceOracle, const std::string& name = "");

    /**
     * \brief Destructor.
     */

    virtual ~ProjectionRealOptimizationOracle();

    /**
     * \brief Adds variable index \p sourceVariableIndex to the oracle.
     **/

    void addVariable(std::size_t sourceVariableIndex);

    /**
     * \brief Adds all those variables from the source oracle that match \p regex.
     **/

    void addVariables(const std::string& regex);

  protected:

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    virtual Response maximize(const double* objectiveVector, const Query& query);

  protected:

    std::shared_ptr<RealOptimizationOracle> _sourceOracle; /**< Source oracle. */
    std::vector<sparse_vector<double> > _projectionLinear; /**< Vector containing the linear part of the projection. */
    std::vector<double> _projectionConstant; /**< Vector containing the absolute part of the projection. */
  };

#if defined(IPO_WITH_GMP)

  class ProjectionRationalOptimizationOracle: public RationalOptimizationOracle
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

    ProjectionRationalOptimizationOracle(std::shared_ptr<RationalOptimizationOracle> sourceOracle,
      const std::string& name = "");

    /**
     * \brief Destructor.
     */

    virtual ~ProjectionRationalOptimizationOracle();

    /**
     * \brief Adds variable index \p sourceVariableIndex to the oracle.
     **/

    void addVariable(std::size_t sourceVariableIndex);

    /**
     * \brief Adds all those variables from the source oracle that match \p regex.
     **/

    void addVariables(const std::string& regex);

  protected:

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    virtual Response maximize(const mpq_class* objectiveVector, const Query& query);

  protected:

    std::shared_ptr<RationalOptimizationOracle> _sourceOracle; /**< Source oracle. */
    std::vector<sparse_vector<mpq_class> > _projectionLinear; /**< Vector containing the linear part of the projection. */
    std::vector<mpq_class> _projectionConstant; /**< Vector containing the absolute part of the projection. */
  };

#endif /* IPO_WITH_GMP */

} /* namespace ipo */

