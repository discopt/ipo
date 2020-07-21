#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>

#include <memory>

#if defined(IPO_WITH_GMP) && defined(IPO_WITH_SOPLEX)

#define SOPLEX_WITH_GMP
#include <soplex.h>
#include <ipo/oracles.hpp>

namespace ipo
{
  class RationalMIPExtender
  {
  public:

    IPO_EXPORT
    RationalMIPExtender(const std::vector<bool>& integrality,
      const std::vector<std::pair<double, double>>& bounds);

    IPO_EXPORT
    virtual ~RationalMIPExtender();

    IPO_EXPORT
    void addConstraint(const Constraint<mpq_class>& constraint);

    IPO_EXPORT
    void addConstraint(const Constraint<double>& constraint);

    IPO_EXPORT
    void setFace(std::shared_ptr<Constraint<mpq_class>> face);

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<mpq_class>::Response maximizeDouble(
      std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const double* objectiveVector,
      const OptimizationOracle<mpq_class>::Query& query);

    /**
     * \brief Maximize a rational objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual OptimizationOracle<mpq_class>::Response maximize(
      std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      const mpq_class* objectiveVector,
      const OptimizationOracle<mpq_class>::Query& query);

  private:
    void setZeroObjective();

    void prepareRay();

    void extractRay(OptimizationOracle<mpq_class>::Response& result);

    void preparePoint(
      const OptimizationOracle<double>::Response::Point& approximatePoint);

    void extractPoint(OptimizationOracle<mpq_class>::Response& result,
      const mpq_class* objectiveVector);

  protected:
    soplex::SoPlex _spx;
    std::vector<bool> _integrality;
    int* _indices;
    mpq_t* _coefficients;
    mpq_t* _originalLowerBounds;
    mpq_t* _originalUpperBounds;
    std::shared_ptr<Constraint<mpq_class>> _currentFace;
  };

  class RationalMIPExtendedOptimizationOracle: public OptimizationOracle<mpq_class>
  {
  public:
    IPO_EXPORT
    RationalMIPExtendedOptimizationOracle(RationalMIPExtender* extender,
      std::shared_ptr<OptimizationOracle<double>> approximateOracle,
      std::shared_ptr<Constraint<mpq_class>> face);

    IPO_EXPORT
    virtual ~RationalMIPExtendedOptimizationOracle();

    IPO_EXPORT
    OptimizationOracle<mpq_class>::Response maximizeDouble(const double * objectiveVector,
      const OptimizationOracle<mpq_class>::Query& query) override;

    IPO_EXPORT
    OptimizationOracle<mpq_class>::Response maximize(const mpq_class* objectiveVector,
      const OptimizationOracle<mpq_class>::Query& query) override;

  protected:
    RationalMIPExtender* _extender;
    std::shared_ptr<OptimizationOracle<double>> _approximateOracle;
    std::shared_ptr<Constraint<mpq_class>> _face;
  };

  class RationalMIPExtendedSeparationOracle: public SeparationOracle<mpq_class>
  {
  public:
    IPO_EXPORT
    RationalMIPExtendedSeparationOracle(std::shared_ptr<SeparationOracle<double>> approximateOracle,
      std::shared_ptr<Constraint<mpq_class>> face);

    IPO_EXPORT
    virtual ~RationalMIPExtendedSeparationOracle();

    /**
     * \brief Returns initially known inequalities.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    IPO_EXPORT
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

    IPO_EXPORT
    virtual SeparationOracle<mpq_class>::Response separateDouble(const double* vector, bool isPoint,
      const SeparationOracle<mpq_class>::Query& query) override;

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
    virtual SeparationOracle<mpq_class>::Response separate(const mpq_class* vector, bool isPoint,
      const SeparationOracle<mpq_class>::Query& query) override;

  protected:
    std::shared_ptr<SeparationOracle<double>> _approximateOracle;
    std::shared_ptr<Constraint<mpq_class>> _face;
  };

} /* namespace ipo */

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
