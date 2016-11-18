#ifndef IPO_POLYHEDRON_H_
#define IPO_POLYHEDRON_H_

#include <cmath>

#include "common.h"
#include "oracles.h"
#include "unique_vectors.h"

namespace ipo {

  class Polyhedron
  {
  protected:
    class CollectOracle : public OracleBase
    {
    public:
      /**
      * \brief Constructs an oracle that saves computed points, rays and valid inequalities.
      *
      * Constructs an oracle that saves computed points, rays and valid inequalities.
      */

      CollectOracle(const std::shared_ptr<OracleBase>& nextOracle);

      /**
      * \brief Destructor.
      *
      * Destructor.
      */

      virtual ~CollectOracle();

      /**
      * \brief Restricts the oracle to the face defined by \p newFace.
      *
      * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
      * For \p newFace equal to \c NULL we define \f$ F := P \f$.
      *
      * This implementation stores \p newFace such that currentFace() works properly and calls
      * setFace() for the next oracle.
      */

      virtual void setFace(const LinearConstraint& newFace = completeFaceConstraint());

      /**
      * \brief Wrapper method that calls the oracle's implementation.
      *
      * This method is called by maximize(), forwards the call to the next oracle if requested, calls the
      * maximizeImplementation() method, and finally forwards the call to the next oracle if necessary.
      *
      * \param result         After the call, contains the oracle's answer.
      * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
      * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
      * \param sort           Set this variable to true if points must be sorted.
      * \param checkDups      Set this variable to true if points or rays must be checked for duplicates.
      *
      * This implementation calls the standard implementation from \ref OracleBase and adds new points and rays to the cache.
      */

      virtual HeuristicLevel maximizeController(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
        bool& checkDups);

      /**
      * \brief Oracle's implementation to maximize the dense rational \p objective.
      *
      * This method is called by maximizeController() and contains the implementation of the oracle.
      *
      * \param result         After the call, contains the oracle's answer.
      * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
      * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
      * \param sort           Set this variable to true if points must be sorted.
      * \param checkDups      Set this variable to true if points or rays must be checked for duplicates.
      *
      * This implementation searches in the point / ray storage for suitable results.
      */

      virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
        bool& checkDups);

      friend class Polyhedron;

    protected:
      UniqueVectors _points; // Points found so far.
      UniqueVectors _rays; // Normalized rays found so far.
      UniqueVectors _normals; // Normals of normalized inequalities.
      std::vector<LinearConstraint> _inequalities; // Normalized inequalities found so far.
    };

    class Face
    {
    public:
      virtual ~Face();

      inline const LinearConstraint& inequality() const
      {
        return _inequality;
      }

      inline const AffineOuterDescription& outerDescription() const
      {
        return _outerDescription;
      }

      inline const InnerDescription& innerDescription() const
      {
        return _innerDescription;
      }

      inline const int dimension() const
      {
        return int(_innerDescription.points.size() + _innerDescription.rays.size()) - 1;
      }
      
      friend class Polyhedron;

    protected:

      Face(LinearConstraint& inequality);

      LinearConstraint _inequality;
      AffineOuterDescription _outerDescription;
      InnerDescription _innerDescription;
    };

  public:

    Polyhedron(const std::shared_ptr<OracleBase>& oracle);

    virtual ~Polyhedron();

    inline Space space() const
    {
      return _collectOracle->space();
    }

    inline std::size_t numPoints() const
    {
      return _collectOracle->_points.size();
    }

    inline std::size_t numRays() const
    {
      return _collectOracle->_rays.size();
    }

    inline std::size_t numInequalities() const
    {
      return _collectOracle->_inequalities.size();
    }

    void affineHull();

    inline int dimension()
    {
      affineHull();
      return _faces.front().dimension();
    }

    inline const AffineOuterDescription& affineHullOuterDescription()
    {
      affineHull();
      return _faces.front().outerDescription();
    }

    inline const InnerDescription& affineHullInnerDescription()
    {
      affineHull();
      return _faces.front().innerDescription();
    }

  protected:
    std::shared_ptr<CollectOracle> _collectOracle;

    VectorMap<Face> _faces;

    HeuristicLevel _affineHullLastCheapHeuristic;
    HeuristicLevel _affineHullLastModerateHeuristic;
    bool _affineHullApproximateDirections;
  };


} /* namespace ipo */

#endif /* IPO_POLYHEDRON_H_ */
