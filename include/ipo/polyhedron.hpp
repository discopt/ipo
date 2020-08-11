#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>

#include <memory>

namespace ipo
{

  class RealPolyhedron: public std::enable_shared_from_this<RealPolyhedron>
  {
  public:
    typedef double Number;
    typedef RealOptimizationOracle OptimizationOracle;
    typedef RealSeparationOracle SeparationOracle;

    /// Deleted default constructor.

    RealPolyhedron() = delete;

    /**
     * \brief Creates a polyhedron in given ambient \p space without any oracles.
     */

    IPO_EXPORT
    RealPolyhedron(std::shared_ptr<Space> space);

    /**
     * \brief Creates the polyhedron defined by the \p optimizationOracle.
     */

    IPO_EXPORT
    RealPolyhedron(std::shared_ptr<RealOptimizationOracle> optimizationOracle);

    /**
     * \brief Creates the polyhedron defined by the \p separationOracle.
     */

    IPO_EXPORT
    RealPolyhedron(std::shared_ptr<RealSeparationOracle> separationOracle);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~RealPolyhedron();

    /**
     * \brief Returns the polyhedron's ambient space.
     */

    IPO_EXPORT
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    RealOptimizationOracle::Response maximize(const double* objectiveVector,
      const RealOptimizationOracle::Query& query);

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if point was not cached before.
     */

    IPO_EXPORT
    bool cachePoint(std::shared_ptr<sparse_vector<double>> point);

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if ray was not cached before.
     */

    IPO_EXPORT
    bool cacheRay(std::shared_ptr<sparse_vector<double>> ray);

    IPO_EXPORT
    std::size_t numCachedSolutions() const;

    IPO_EXPORT
    std::pair<bool, std::shared_ptr<sparse_vector<double>>> getCachedSolution(std::size_t index);

  protected:

    /// Pointer to private implementation.
    void* _implementation;
    /// Ambient space.
    std::shared_ptr<Space> _space;
  };

#if defined(IPO_WITH_GMP)

  class RationalPolyhedron: public std::enable_shared_from_this<RationalPolyhedron>
  {
  public:
    typedef double Number;
    typedef RationalOptimizationOracle OptimizationOracle;
    typedef RationalSeparationOracle SeparationOracle;

    /// Deleted default constructor.

    RationalPolyhedron() = delete;

    /**
     * \brief Creates a polyhedron in given ambient \p space without any oracles.
     */

    IPO_EXPORT
    RationalPolyhedron(std::shared_ptr<Space> space);

    /**
     * \brief Creates the polyhedron defined by the \p optimizationOracle.
     */

    IPO_EXPORT
    RationalPolyhedron(std::shared_ptr<RationalOptimizationOracle> optimizationOracle);

    /**
     * \brief Creates the polyhedron defined by the \p separationOracle.
     */

    IPO_EXPORT
    RationalPolyhedron(std::shared_ptr<RationalSeparationOracle> separationOracle);

    /**
     * \brief Destructor.
     */

    IPO_EXPORT
    virtual ~RationalPolyhedron();

    /**
     * \brief Returns the polyhedron's ambient space.
     */

    IPO_EXPORT
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

    /**
     * \brief Maximize a rational objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    RationalOptimizationOracle::Response maximize(const mpq_class* objectiveVector,
      const RationalOptimizationOracle::Query& query);

    /**
     * \brief Maximize a floating-point objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    RationalOptimizationOracle::Response maximize(const double* objectiveVector,
      const RationalOptimizationOracle::Query& query);

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if point was not cached before.
     */

    IPO_EXPORT
    bool cachePoint(std::shared_ptr<sparse_vector<mpq_class>> point);

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if ray was not cached before.
     */

    IPO_EXPORT
    bool cacheRay(std::shared_ptr<sparse_vector<mpq_class>> ray);

  protected:
    /// Pointer to private implementation.
    void* _implementation;
    /// Ambient space.
    std::shared_ptr<Space> _space;
  };

#endif /* IPO_WITH_GMP */

} /* namespace ipo */

