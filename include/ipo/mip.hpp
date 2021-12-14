#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#if defined(IPO_RATIONAL_MIP) && defined(IPO_RATIONAL_LP)

#include <memory>

#include <ipo/space.hpp>
#include <ipo/oracles.hpp>

namespace ipo
{
  class RationalMIPExtender
  {
  public:

    RationalMIPExtender(const std::vector<bool>& integrality, const std::vector<std::pair<double, double>>& bounds);

    virtual ~RationalMIPExtender();

    void addConstraint(const Constraint<rational>& constraint);

    void addConstraint(const Constraint<double>& constraint);

    void setFace(Constraint<rational>* face);

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    virtual OptimizationOracle<rational>::Response maximizeDouble(
      std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const double* objectiveVector, const OptimizationOracle<rational>::Query& query);

    /**
     * \brief Maximize a rational objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    virtual OptimizationOracle<rational>::Response maximize(
      std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const rational* objectiveVector,
      const OptimizationOracle<rational>::Query& query);

  protected:
    void* _implementation;
  };

  class RationalMIPExtendedOptimizationOracle: public OptimizationOracle<rational>
  {
  public:
    RationalMIPExtendedOptimizationOracle(RationalMIPExtender* extender,
      std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const Constraint<rational>& face);

    virtual ~RationalMIPExtendedOptimizationOracle();

    OptimizationOracle<rational>::Response maximize(const rational* objectiveVector,
      const OptimizationOracle<rational>::Query& query) override;

    OptimizationOracle<rational>::Response maximizeDouble(const double* objectiveVector,
      const OptimizationOracle<rational>::Query& query) override;

  protected:
    RationalMIPExtender* _extender;
    std::shared_ptr<OptimizationOracle<double>> _approximateOracle;
    Constraint<rational> _face;
  };

  class RationalMIPExtendedSeparationOracle: public SeparationOracle<rational>
  {
  public:
    RationalMIPExtendedSeparationOracle(std::shared_ptr<SeparationOracle<double>> approximateOracle,
      const Constraint<rational>& face);

    virtual ~RationalMIPExtendedSeparationOracle();

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    virtual Response getInitial(const Query& query);

    /**
     * \brief Separates a point/ray with floating-point coordinates.
     *
     * Separates a point/ray with floating-point coordinates. This default implementation converts
     * the double objective to one of type \ref T.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param query Structure for query.
     * \param isPoint Whether a point shall be separated.
     * \param result Structure for returning the result.
     *
     * \returns \c true if and only if the point/ray was separated.
     **/

    virtual SeparationResponse<rational> separateDouble(const double* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery()) override;

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

    virtual SeparationResponse<rational> separate(const rational* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery()) override;

  protected:
    std::shared_ptr<SeparationOracle<double>> _approximateOracle;
    Constraint<rational> _face;
  };

} /* namespace ipo */

#endif /* IPO_RATIONAL_MIP && IPO_RATIONAL_LP */
