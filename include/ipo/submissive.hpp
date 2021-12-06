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
  class SubmissiveOptimizationOracle: public OptimizationOracle<NumberType>
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
    SubmissiveOptimizationOracle(std::shared_ptr<OptimizationOracle<NumberType>> sourceOracle,
      const std::string& name = "");

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~SubmissiveOptimizationOracle();

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

    std::shared_ptr<OptimizationOracle<Number>> _sourceOracle; /**< Source oracle. */
  };

} /* namespace ipo */
