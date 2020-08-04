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
    void setFace(Constraint<mpq_class>* face);

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual RationalOptimizationOracle::Response maximize(
      std::shared_ptr<RealOptimizationOracle> approximateOracle,
      const double* objectiveVector,
      const RationalOptimizationOracle::Query& query);

    /**
     * \brief Maximize a rational objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Parameters of query.
     * \return Optimization result.
     **/

    IPO_EXPORT
    virtual RationalOptimizationOracle::Response maximize(
      std::shared_ptr<RealOptimizationOracle> approximateOracle,
      const mpq_class* objectiveVector,
      const RationalOptimizationOracle::Query& query);

  private:
    void setZeroObjective();

    void prepareRay();

    void extractRay(RationalOptimizationOracle::Response& result);

    void preparePoint(
      const RealOptimizationOracle::Response::Point& approximatePoint);

    void extractPoint(RationalOptimizationOracle::Response& result,
      const mpq_class* objectiveVector);

  protected:
    soplex::SoPlex _spx;
    std::vector<bool> _integrality;
    int* _indices;
    mpq_t* _coefficients;
    mpq_t* _originalLowerBounds;
    mpq_t* _originalUpperBounds;
    Constraint<mpq_class>* _currentFace;
  };

  class RationalMIPExtendedOptimizationOracle: public RationalOptimizationOracle
  {
  public:
    IPO_EXPORT
    RationalMIPExtendedOptimizationOracle(RationalMIPExtender* extender,
      std::shared_ptr<RealOptimizationOracle> approximateOracle,
      const Constraint<mpq_class>& face);

    IPO_EXPORT
    virtual ~RationalMIPExtendedOptimizationOracle();

    IPO_EXPORT
    RationalOptimizationOracle::Response maximize(const mpq_class* objectiveVector,
      const RationalOptimizationOracle::Query& query) override;

    IPO_EXPORT
    RationalOptimizationOracle::Response maximize(const double* objectiveVector,
      const RationalOptimizationOracle::Query& query) override;

  protected:
    RationalMIPExtender* _extender;
    std::shared_ptr<RealOptimizationOracle> _approximateOracle;
    Constraint<mpq_class> _face;
  };

  class RationalMIPExtendedSeparationOracle: public RationalSeparationOracle
  {
  public:
    IPO_EXPORT
    RationalMIPExtendedSeparationOracle(std::shared_ptr<RealSeparationOracle> approximateOracle,
      const Constraint<mpq_class>& face);

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
    virtual RationalSeparationOracle::Response separate(const double* vector, bool isPoint,
      const RationalSeparationOracle::Query& query) override;

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
    virtual RationalSeparationOracle::Response separate(const mpq_class* vector, bool isPoint,
      const RationalSeparationOracle::Query& query) override;

  protected:
    std::shared_ptr<RealSeparationOracle> _approximateOracle;
    Constraint<mpq_class> _face;
  };

} /* namespace ipo */

#endif /* IPO_WITH_GMP && IPO_WITH_SOPLEX */
