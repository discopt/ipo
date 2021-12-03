#pragma once

#include <iostream>

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/oracles.hpp>

#include <memory>

namespace ipo
{
  template <typename NumberType>
  class Polyhedron: public std::enable_shared_from_this<Polyhedron<NumberType>>, OptimizationOracle<NumberType>,
    SeparationOracle<NumberType>
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

    IPO_EXPORT
    Polyhedron(std::shared_ptr<Space> space, const std::string& name = "");

    /**
     * \brief Creates the polyhedron defined by the \p optimizationOracle.
     **/

    IPO_EXPORT
    Polyhedron(std::shared_ptr<OptOracle> optimizationOracle);

    /**
     * \brief Creates the polyhedron defined by the \p separationOracle.
     **/

    IPO_EXPORT
    Polyhedron(std::shared_ptr<SepaOracle> separationOracle);

    /**
     * \brief Destructor.
     **/

    IPO_EXPORT
    virtual ~Polyhedron();

    /**
     * \brief Maximize an objective vector.
     *
     * \param objectiveVector Objective vector.
     * \param query Additional query information.
     * \return Optimization response.
     **/

    IPO_EXPORT
    virtual OptimizationResponse<Number> maximize(const Number* objectiveVector,
      const OptimizationQuery<Number>& query) override;

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if point was not cached before.
     **/

    IPO_EXPORT
    bool cachePoint(std::shared_ptr<sparse_vector<Number>> point);

    /**
     * \brief Adds given point to cache.
     *
     * Checks for duplicates.
     * 
     * \returns \c true if ray was not cached before.
     **/

    IPO_EXPORT
    bool cacheRay(std::shared_ptr<sparse_vector<Number>> ray);

    /**
     * \brief Returns the number of cached points plus rays.
     **/

    IPO_EXPORT
    std::size_t numCachedSolutions() const;

    /**
     * \brief Returns a pair of an indicator and the solution \p index, where the indicator is \c true if it is a point.
     **/

    IPO_EXPORT
    std::pair<bool, std::shared_ptr<sparse_vector<Number>>> getCachedSolution(std::size_t index);

    /**
     * \brief Tries to separate the point or ray \p vector from the polyhedron by a hyperplane (inequality).
     **/

    IPO_EXPORT
    virtual SeparationResponse<Number> separate(const Number* vector, bool isPoint,
      const SeparationQuery& query = SeparationQuery()) override;

    IPO_EXPORT
    inline std::shared_ptr<Space> space() const
    {
      return _space;
    }

  protected:

    /// Pointer to private implementation.
    void* _implementation;
    std::shared_ptr<Space> _space;
  };


//   class RealPolyhedron: public std::enable_shared_from_this<RealPolyhedron>, OptimizationOracle<double>,
//     SeparationOracle<double>
//   {
//   public:
//     typedef double Number;
//     typedef OptimizationOracle<double> OptOracle;
//     typedef SeparationOracle<double> SepaOracle;
// 
//     /// Deleted default constructor.
// 
//     RealPolyhedron() = delete;
// 
//     /**
//      * \brief Creates a polyhedron in given ambient \p space without any oracles.
//      **/
// 
//     IPO_EXPORT
//     RealPolyhedron(std::shared_ptr<Space> space, const std::string& name);
// 
//     /**
//      * \brief Creates the polyhedron defined by the \p optimizationOracle.
//      **/
// 
//     IPO_EXPORT
//     RealPolyhedron(std::shared_ptr<OptOracle> optimizationOracle);
// 
//     /**
//      * \brief Creates the polyhedron defined by the \p separationOracle.
//      **/
// 
//     IPO_EXPORT
//     RealPolyhedron(std::shared_ptr<SepaOracle> separationOracle);
// 
//     /**
//      * \brief Destructor.
//      **/
// 
//     IPO_EXPORT
//     virtual ~RealPolyhedron();
// 
//     /**
//      * \brief Returns the polyhedron's ambient space.
//      **/
// 
//     IPO_EXPORT
//     inline std::shared_ptr<Space> space() const
//     {
//       return _space;
//     }
// 
//     /**
//      * \brief Maximize a floating-point objective vector.
//      *
//      * \param objectiveVector Objective vector.
//      * \param query Additional query information.
//      * \return Optimization response.
//      **/
// 
//     IPO_EXPORT
//     OptOracle::Response maximize(const double* objectiveVector,
//       const OptOracle::Query& query);
// 
//     /**
//      * \brief Adds given point to cache.
//      *
//      * Checks for duplicates.
//      * 
//      * \returns \c true if point was not cached before.
//      **/
// 
//     IPO_EXPORT
//     bool cachePoint(std::shared_ptr<sparse_vector<double>> point);
// 
//     /**
//      * \brief Adds given point to cache.
//      *
//      * Checks for duplicates.
//      * 
//      * \returns \c true if ray was not cached before.
//      **/
// 
//     IPO_EXPORT
//     bool cacheRay(std::shared_ptr<sparse_vector<double>> ray);
// 
//     /**
//      * \brief Returns the number of cached points plus rays.
//      **/
// 
//     IPO_EXPORT
//     std::size_t numCachedSolutions() const;
// 
//     /**
//      * \brief Returns a pair of an indicator and the solution \p index, where the indicator is \c true if it is a point.
//      **/
// 
//     IPO_EXPORT
//     std::pair<bool, std::shared_ptr<sparse_vector<double>>> getCachedSolution(std::size_t index);
// 
//     /**
//      * \brief Tries to separate the point or ray \p vector from the polyhedron by a hyperplane (inequality).
//      **/
// 
//     IPO_EXPORT
//     SeparationResponse<double> separate(const double* vector, bool isPoint,
//       const SeparationQuery& query = SeparationQuery());
// 
//   protected:
// 
//     /// Pointer to private implementation.
//     void* _implementation;
//     /// Ambient space.
//     std::shared_ptr<Space> _space;
//   };
// 
// #if defined(IPO_WITH_GMP)
// 
//   class RationalPolyhedron: public std::enable_shared_from_this<RationalPolyhedron>
//   {
//   public:
//     typedef mpq_class Number;
//     typedef OptimizationOracle<mpq_class> OptOracle;
//     typedef SeparationOracle<mpq_class> SepaOracle;
// 
//     /// Deleted default constructor.
// 
//     RationalPolyhedron() = delete;
// 
//     /**
//      * \brief Creates a polyhedron in given ambient \p space without any oracles.
//      **/
// 
//     IPO_EXPORT
//     RationalPolyhedron(std::shared_ptr<Space> space);
// 
//     /**
//      * \brief Creates the polyhedron defined by the \p optimizationOracle.
//      **/
// 
//     IPO_EXPORT
//     RationalPolyhedron(std::shared_ptr<OptOracle> optimizationOracle);
// 
//     /**
//      * \brief Creates the polyhedron defined by the \p separationOracle.
//      **/
// 
//     IPO_EXPORT
//     RationalPolyhedron(std::shared_ptr<SepaOracle> separationOracle);
// 
//     /**
//      * \brief Destructor.
//      **/
// 
//     IPO_EXPORT
//     virtual ~RationalPolyhedron();
// 
//     /**
//      * \brief Returns the polyhedron's ambient space.
//      **/
// 
//     IPO_EXPORT
//     inline std::shared_ptr<Space> space() const
//     {
//       return _space;
//     }
// 
//     /**
//      * \brief Maximize a rational objective vector.
//      *
//      * \param objectiveVector Objective vector.
//      * \param query Additional query information.
//      * \return Optimization response.
//      **/
// 
//     IPO_EXPORT
//     OptOracle::Response maximize(const mpq_class* objectiveVector,
//       const OptOracle::Query& query);
// 
//     /**
//      * \brief Maximize a floating-point objective vector.
//      *
//      * \param objectiveVector Objective vector.
//      * \param query Additional query information.
//      * \return Optimization response.
//      **/
// 
//     IPO_EXPORT
//     OptOracle::Response maximize(const double* objectiveVector,
//       const OptOracle::Query& query);
// 
//     /**
//      * \brief Adds given point to cache.
//      *
//      * Checks for duplicates.
//      * 
//      * \returns \c true if point was not cached before.
//      **/
// 
//     IPO_EXPORT
//     bool cachePoint(std::shared_ptr<sparse_vector<mpq_class>> point);
// 
//     /**
//      * \brief Adds given point to cache.
//      *
//      * Checks for duplicates.
//      * 
//      * \returns \c true if ray was not cached before.
//      **/
// 
//     IPO_EXPORT
//     bool cacheRay(std::shared_ptr<sparse_vector<mpq_class>> ray);
// 
//     /**
//      * \brief Returns the number of cached points plus rays.
//      **/
// 
//     IPO_EXPORT
//     std::size_t numCachedSolutions() const;
// 
//     /**
//      * \brief Returns a pair of an indicator and the solution \p index, where the indicator is \c true if it is a point.
//      **/
// 
//     IPO_EXPORT
//     std::pair<bool, std::shared_ptr<sparse_vector<mpq_class>>> getCachedSolution(std::size_t index);
// 
//     /**
//      * \brief Tries to separate the point or ray \p vector from the polyhedron by a hyperplane (inequality).
//      **/
// 
//     IPO_EXPORT
//     SeparationResponse<mpq_class> separate(const mpq_class* vector, bool isPoint,
//       const SeparationQuery& query = SeparationQuery());
// 
//   protected:
//     /// Pointer to private implementation.
//     void* _implementation;
//     /// Ambient space.
//     std::shared_ptr<Space> _space;
//   };
// 
// #endif /* IPO_WITH_GMP */

} /* namespace ipo */

