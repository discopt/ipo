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

      virtual void setFace(const LinearConstraint& newFace = completeFace());

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
    };

  public:
    Polyhedron(const std::shared_ptr<OracleBase>& oracle);
    virtual ~Polyhedron();

    void affineHull();

    inline int dimension()
    {
      affineHull();
      return int(_affineHullInner.points.size() + _affineHullInner.rays.size()) - 1;
    }

    inline const AffineOuterDescription& affineHullOuterDescription()
    {
      affineHull();
      return _affineHullOuter;
    }

    inline const InnerDescription& affineHullInnerDescription()
    {
      affineHull();
      return _affineHullInner;
    }

  protected:
    std::shared_ptr<CollectOracle> _collectOracle;

    bool _affineHullComputed;
    AffineOuterDescription _affineHullOuter;
    InnerDescription _affineHullInner;
    HeuristicLevel _affineHullLastCheapHeuristic;
    HeuristicLevel _affineHullLastModerateHeuristic;
    bool _affineHullApproximateDirections;
  };


} /* namespace ipo */

#endif /* IPO_POLYHEDRON_H_ */
