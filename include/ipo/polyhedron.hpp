#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>

#include <memory>

namespace ipo
{
  template <typename NumberType>
  class Polyhedron: public std::enable_shared_from_this<Polyhedron<NumberType>>,
    virtual public OptimizationOracle<NumberType>,
    virtual public SeparationOracle<NumberType>
  {
  public:
    typedef NumberType Number;
    typedef OptimizationOracle<Number> OptOracle;
    typedef SeparationOracle<Number> SepaOracle;

    /// Deleted default constructor.

    Polyhedron() = delete;

    /**
     * \brief Creates a polyhedron in given ambient \p space without any oracles.
     **/

    Polyhedron(std::shared_ptr<Space> space, const std::string& name = "");

    /**
     * \brief Creates the polyhedron defined by the \p optimizationOracle.
     **/

    Polyhedron(std::shared_ptr<OptOracle> optimizationOracle);

    /**
     * \brief Creates the polyhedron defined by the \p separationOracle.
     **/

    Polyhedron(std::shared_ptr<SepaOracle> separationOracle);

    /**
     * \brief Destructor.
     **/

    virtual ~Polyhedron();

    /**
     * \brief Maximize an objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    virtual OptimizationResponse<Number> maximize(const Number* objectiveVector,
      const OptimizationQuery<Number>& query) override;

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if point was not cached before.
     **/

    bool cachePoint(std::shared_ptr<sparse_vector<Number>> point);

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if ray was not cached before.
     **/

    bool cacheRay(std::shared_ptr<sparse_vector<Number>> ray);

    /**
     * \brief Returns the number of cached points plus rays.
     **/

    std::size_t numCachedSolutions() const;

    /**
     * \brief Returns a pair of an indicator and the solution \p index, where the indicator is \c true if it is a point.
     **/

    std::pair<bool, std::shared_ptr<sparse_vector<Number>>> getCachedSolution(std::size_t index);

    /**
     * \brief Tries to separate the point or ray \p vector from the polyhedron by a hyperplane (inequality).
     **/

    virtual SeparationResponse<Number> separate(const Number* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery()) override;

  protected:

    /// Pointer to private implementation.
    void* _implementation;
  };

} /* namespace ipo */

