#ifndef IPO_CACHE_ORACLE_H_
#define IPO_CACHE_ORACLE_H_

#include <cmath>

#include "ipo.h"
#include "oracles.h"
#include "unique_rational_vectors.h"


namespace ipo {

  class CacheOracle : public OracleBase
  {
  public:
    /**
     * \brief Constructs an oracle that stores points and directions explicitly.
     *
     * Constructs an oracle with given \p name in given \p space that stores points and directions
     * explicitly. The storage is external via the const references \p points and \p directions.
     */

    CacheOracle(const std::string& name, const Space& space,
      const UniqueRationalVectorsBase& points, const UniqueRationalVectorsBase& directions);

    /**
     * \brief Constructs an oracle that stores points and directions explicitly.
     *
     * Constructs an oracle with given \p name in given \p space that stores points and directions
     * explicitly. It is associated to \p nextOracle. The storage is external via the const
     * references \p points and \p directions.
     */

    CacheOracle(const std::string& name, OracleBase* nextOracle,
      const UniqueRationalVectorsBase& points, const UniqueRationalVectorsBase& directions);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~CacheOracle();


    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation stores \p newFace such that currentFace() works properly and calls
     * setFace() for the next oracle.
     */

    virtual void setFace(Face* newFace = NULL);

    /**
     * \brief Runs the oracle to maximize the dense rational \p objective.
     *
     * Runs the optimization oracle to maximize the given dense rational \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p objectiveBound is not exceeded, then the call must be forwarded to the
     * next oracle.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * This implementation searches in the associated point / direction storage.
     */

    virtual void maximize(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

  protected:
    struct VectorStats
    {
      int valueExponent;
      double valueMantissa;
      std::size_t sparsity;
      std::size_t index;

      double value; // TODO: remove

      VectorStats();
      VectorStats(double theObjectiveValue, std::size_t theSparsity, std::size_t theIndex);
      VectorStats& operator=(const VectorStats& other);
      bool operator<(const VectorStats& other) const;
    };

    typedef std::vector<std::size_t> FaceIndices;

    void updateFaceIndices(const UniqueRationalVectorsBase& vectors, FaceIndices& faceIndices,
      std::size_t& end, bool handlingPoints);

    void search(const UniqueRationalVectorsBase& vectors, const FaceIndices& faceIndices,
      const soplex::VectorReal& approxObjective, double approxObjectiveBound, bool handlingPoints,
      std::vector<std::size_t>& result);

    const UniqueRationalVectorsBase& _points;
    FaceIndices _facePoints;
    std::size_t _endFacePoints;

    const UniqueRationalVectorsBase& _directions;
    FaceIndices _faceDirections;
    std::size_t _endFaceDirections;

    std::vector<VectorStats> _vectorStats;
  };

} /* namespace ipo */

#endif /* IPO_CACHE_ORACLE_H_ */
