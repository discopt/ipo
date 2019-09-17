#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>
#include <ipo/oracles.hpp>

#include <memory>
#include <vector>

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

  protected:
    struct QueryStatistics
    {
      double runningTime;
      bool success;
    };

    template <class T>
    struct Data
    {
      double expectedRunningTime;
      double expectedSuccess;
      std::shared_ptr<T> oracle; /// Pointer to oracle.
      std::vector<QueryStatistics> history;

      Data(std::shared_ptr<T> oracle);
    };

    std::vector<Data<OptimizationOracle>> _optimization;
    std::vector<Data<SeparationOracle>> _separation;
    double _normalizedRayEpsilon;
    double _normalizedPointEpsilon;
  };

} /* namespace ipo */

