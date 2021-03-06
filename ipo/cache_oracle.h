#ifndef IPO_CACHE_ORACLE_H_
#define IPO_CACHE_ORACLE_H_

#include <cmath>

#include "common.h"
#include "oracles.h"
#include "unique_vectors.h"

namespace ipo {

  class CacheOracle : public OracleBase
  {
  public:
    enum Behavior
    {
      DISABLED,
      CACHE_ONLY,
      CACHE_AND_SEARCH
    };

    /**
     * \brief Constructs an oracle that stores points and rays explicitly.
     *
     * Constructs an oracle with given \p name that stores points and rays returned by the associated oracle \p nextOracle
     * explicitly.
     *
     * \param nextOracle    Next oracle to be called if cache search was not satisfactory.
     * \param outerBehavior Whether to cache/search points and rays.
     * \param innerBehavior Whether to cache/search valid inequalities.
     */

    CacheOracle(const std::shared_ptr<OracleBase>& nextOracle, Behavior outerBehavior = CACHE_ONLY,
      Behavior innerBehavior = CACHE_AND_SEARCH);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~CacheOracle();

    /**
     * \brief Enables or disables caching/searching of valid inequalities.
     *
     * Enables or disables caching/searching of valid inequalities. Disabling also frees the current cache. Note that searching
     * for valid inequalities requires solving an LP and only effects the returned heuristicLevel (\sa OracleBase).
     */

    void setOuterBehavior(Behavior outerBehavior);

    /**
     * \brief Returns the current caching/searching behavior for points and rays.
     *
     * Returns the current caching/searching behavior for points and rays.
     */

    inline Behavior outerBehavior() const
    {
      return _outerBehavior;
    }

    /**
     * \brief Enables or disables caching/searching of points and rays.
     *
     * Enables or disables caching/searching of points and rays. Disabling also frees the current cache.
     */

    void setInnerBehavior(Behavior outerBehavior);

    /**
     * \brief Returns the current caching/searching behavior for points and rays.
     *
     * Returns the current caching/searching behavior for points and rays.
     */

    inline Behavior innerBehavior() const
    {
      return _innerBehavior;
    }

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

    bool addPoint(const Vector& point);

    bool addRay(const Vector& ray);

    inline std::size_t numPoints() const
    {
      return _uniquePoints.size();
    }

    inline std::size_t numRays() const
    {
      return _uniqueRays.size();
    }

  protected:
    struct Data
    {
      Vector vector;

      int valueExponent;
      double valueMantissa;

      Data(Vector& vector);
      ~Data();

      void updateObjective(const soplex::VectorReal& approximateObjective, double approximateObjectiveBound);
      bool operator<(const Data& other) const;
    };

    void search(std::vector<Data>& vectors, const soplex::VectorReal& approximateObjective,
      double approximateObjectiveBound, bool handlingPoints, std::vector<Vector>& result);

  protected:
    Behavior _outerBehavior;
    Behavior _innerBehavior;

    // Structures for inner caching.

    UniqueVectors _uniquePoints;
    UniqueVectors _uniqueRays;
    std::vector<Data> _facePoints;
    std::vector<Data> _faceRays;

    // Structures for outer caching.

    soplex::SoPlex _inequalities;
  };

} /* namespace ipo */

#endif /* IPO_CACHE_ORACLE_H_ */
