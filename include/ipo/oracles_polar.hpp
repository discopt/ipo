#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/polyhedron.hpp>
#include <ipo/constraint.hpp>
#include <ipo/affine_hull.hpp>

#include <vector>
#include <ostream>

#if defined(IPO_DOUBLE_LP) || defined(IPO_RATIONAL_LP)

namespace ipo
{
  
  
//   struct SeparationQuery
//   {
//     /**
//      * \brief Decides when a constraint is satisfied with equality for double-precision
//      * computations.
//      */
//     double epsilonConstraints;
//     double timeLimit; /// Time limit in seconds.
// 
//     SeparationQuery();
// 
//     SeparationQuery& operator=(const SeparationQuery& other);
//   };

  const static int SEPARATION_ERROR_RUNNING = -2;
  const static int SEPARATION_ERROR_TIMEOUT = -3;
  const static int SEPARATION_ERROR_NUMERICS = -4;

  template <typename NumberType>
  class PolarSeparationOracle : public SeparationOracle<NumberType>
  {
  public:
    typedef NumberType Number;

    /// Query structure.
    typedef SeparationQuery Query;
    /// Response structure.
    typedef SeparationResponse<Number> Response;

    PolarSeparationOracle(std::shared_ptr<OptimizationOracle<Number>> optOracle,
      std::shared_ptr<SeparationOracle<Number>> sepaRelaxationOracle = NULL, const std::string& name = "");

    void setAffineHull(const AffineHull<Number>& affineHull);

    /**
     * \brief Returns initially known inequalities.
     *
     * The default implementation returns nothing.
     *
     * \param query Structure for query.
     * \return Separation result.
     **/

    virtual Response getInitial(const Query& query = Query());

    /**
     * \brief Separates a point/ray.
     *
     * \param vector Array that maps coordinates to point/ray coordinates.
     * \param isPoint Whether a point shall be separated.
     * \param query Additional query information.
     *
     * \returns Reponse structure.
     **/

    virtual Response separate(const Number* vector, bool isPoint, const Query& query = Query());

  protected:
    void* _implementation;
  };

//   template <typename NumberType>
//   struct Separated
//   {
//     typedef NumberType Number;
// 
//     Constraint<Number> constraint;                              /// The computed facet or equation.
//     double timeTotal;                                           /// Total time spent.
//     std::vector<std::shared_ptr<sparse_vector<Number>>> points; /// Points of inner description.
//     std::vector<std::shared_ptr<sparse_vector<Number>>> rays;   /// Rays of inner description.
// 
//     Separated();
// 
//     Separated(Separated<Number>&& other);
// 
//     Separated<Number>& operator=(const Separated<Number>& other);
// 
//     Separated<Number>& operator=(Separated<Number>&& other);
//   };
// 
//   template <typename Number>
//   Separated<Number> separatePoint(std::shared_ptr<Polyhedron<Number>> polyhedron,
//     const AffineHullResult<Number>& affineHullResult, const Number* point, const FacetQuery& query = FacetQuery());
// 
//   template <typename Number>
//   Separated<Number> separateRay(std::shared_ptr<Polyhedron<Number>> polyhedron,
//     const AffineHullResult<Number>& affineHullResult, const Number* ray, const FacetQuery& query = FacetQuery());

} /* namespace ipo */

#endif /* IPO_DOUBLE_LP || IPO_RATIONAL_LP */
