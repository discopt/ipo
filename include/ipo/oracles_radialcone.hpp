#pragma once

#include <ipo/oracles.hpp>

namespace ipo
{

  template <typename NumberType>
  class RadialConeOptimizationOracle: virtual public OptimizationOracle<NumberType>
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

    RadialConeOptimizationOracle(std::shared_ptr<OptimizationOracle<Number>> sourceOracle,
      std::shared_ptr<sparse_vector<Number>> apex, const std::string& name = "");

    /**
     * \brief Destructor.
     */

    virtual ~RadialConeOptimizationOracle();

    /**
     * \brief Returns whether one can exploit a trust region.
     */

    bool isTrustRegionCapable() const;

    /**
     * \brief Enables usage of a trust region.
     *
     * The trust region's size will be \p initialSize. Upon each failure, its size it increased by a factor of
     * \p growthRate. If its size reaches \p maximumSize then it is disabled. Every \p resetIteratons iterations (if
     * nonzero) the size is reset to \p initialSize.
     */

    void enableManhattanTrustRegion(const Number& initialSize, const Number& maximumSize, const Number& growthRate,
      size_t resetIterations);

    /**
     * \brief Disables usage of a trust region.
     */

    void disableManhattanTrustRegion();

    /**
     * \brief Maximize an objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    virtual OptimizationResponse<Number> maximize(const Number* objectiveVector,
      const OptimizationQuery<Number>& query);

  protected:
    /// Source oracle.
    std::shared_ptr<OptimizationOracle<Number>> _sourceOracle;
    /// Source oracle as a trust region oracle.
    TrustRegionOptimizationOracle<Number>* _trustRegionOracle;
    /// Apex of the cone.
    std::shared_ptr<sparse_vector<Number>> _apex;
    /// Whether trust regions are enabled.
    bool _trustRegionEnabled;
    /// Current trust region.
    ManhattanTrustRegion<Number> _trustRegion;
    std::size_t _iteration;
    /// Number of iterations for resetting.
    std::size_t _trustRegionResetIterations;
    /// Initial size of the trust region.
    Number _trustRegionInitialSize;
    /// Maximum size of the trust region.
    Number _trustRegionMaximumSize;
    /// Grwoth rate for trust region size.
    Number _trustRegionGrowthRate;
  };


} /* namespace ipo */
