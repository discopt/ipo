#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/oracles.hpp>

#include <memory>
#include <deque>

namespace ipo
{

  class Polyhedron : public std::enable_shared_from_this<Polyhedron>
  {
  public:
    Polyhedron(std::shared_ptr<OptimizationOracle> optimizationOracle);

    Polyhedron(std::shared_ptr<SeparationOracle> separationOracle);

    ~Polyhedron();

    /**
     * \brief Returns the polyhedron's ambient space.
     */

    IPO_EXPORT
    std::shared_ptr<Space> space() const;

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     * \param result Structure for returning the result.
     **/

    void maximize(const double* objectiveVector, const OptimizationOracle::Query& query,
      OptimizationOracle::Result& result);

#if defined(IPO_WITH_GMP)

    /**
     * \brief Maximize a rational objective vector.
     *
     * Maximize a rational objective vector. The default implementation just converts the
     * objective vector to a floating-point vector and calls \ref maximize.
     *
     * \param objectiveVector Array that maps coordinates to objective value coefficients.
     * \param query Structure for query.
     * \param result Structure for returning the result.
     **/

    IPO_EXPORT
    void maximize(const mpq_class* objectiveVector, const OptimizationOracle::Query& query,
      OptimizationOracle::Result& result);

#endif /* IPO_WITH_GMP */

  protected:
    void tuneOracles();
    
  protected:
    struct QueryStatistics
    {
      double runningTime;
      bool success;

      QueryStatistics(double runningTime, bool success);
    };

    template <class T>
    struct Data
    {
      double sumRunningTime; /// Sum of running times of history.
      unsigned int sumSuccess; /// Number of successful runs.
      double priority; /// \ref sumRunningTime / (\ref sumSuccess + 1)
      std::shared_ptr<T> oracle; /// Pointer to oracle.
      std::deque<QueryStatistics> history;
      bool isCache; /// \c true if this is a cache oracle.

      Data(std::shared_ptr<T> oracle, bool isCache = false);

      bool operator<(const Data<T>& other) const;

      void swap(Data<T>& other);
      
      void updateHistory(double runningTime, bool success, std::size_t historySize);
    };

    std::vector<Data<OptimizationOracle>> _optimization;
    std::vector<Data<SeparationOracle>> _separation;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
    std::size_t _historySize;
  };

} /* namespace ipo */

