#ifndef IPO_ORACLES_H_
#define IPO_ORACLES_H_

#include <vector>
#include <limits>

#include "ipo.h"
#include "rows.h"
#include "space.h"
#include "unique_rational_vectors.h"

namespace ipo {

  class SparseVectorDuplicateChecker
  {
  public:
    SparseVectorDuplicateChecker(std::size_t dimension);
    ~SparseVectorDuplicateChecker();

    bool operator()(const soplex::SVectorRational* first, const soplex::SVectorRational* second);

  protected:
    std::vector<soplex::Rational const*> _firstValues; // Dense first vector as entry-pointers.
  };

  /**
   * \brief Defines a face of a polyhedron by a set of inequalities.
   *
   * Defines a face \f$F\f$ of a polyhedron \f$P\f$ by a set of inequalities.
   * It is used to control over which face an \ref OracleBase optimizes.
   **/

  class Face
  {
  public:
    /**
     * \brief Constructs the trivial face.
     *
     * Constructs the trivial face in \c space defined by \f$\left<\mathbbm{O},x\right> \leq 0\f$.
     * Use add() methods to add further inequalities.
     **/

    Face(const Space& space);

    /**
     * \brief Constructs a face definde by \c inequality.
     *
     * Constructs the face in \c space defined by \c inequality.
     * Use add() methods to add further inequalities,
     * that is, to intersect with other faces.
     **/

    Face(const Space& space, const soplex::LPRowRational& inequality);

    /**
     * \brief Constructs the interesection of the faces defined by \c inequalities.
     *
     * Constructs the intersection of the faces defined by \c inequalities in \c space.
     * Use add() methods to add further inequalities,
     * that is, to intersect with other faces.
     **/

    Face(const Space& space, const soplex::LPRowSetRational& inequalities);

    /**
     * \brief Destructor.
     */

    virtual ~Face();

    /**
     * \brief Adds \c inequality to the inequalities that define this face.
     *
     * Adds \c inequality to the inequality currently defining this face.
     * Geometrically this means to intersect the current face with the
     * one defined by \c inequality.
     */

    void add(const soplex::LPRowRational& inequality);

    /**
     * \brief Adds all \c inequalities to the ones that define this face.
     *
     * Adds all \c inequalities to the inequality currently defining this face.
     * Geometrically this means to intersect the current face with all
     * faces defined by \c inequalities.
     */

    void add(const soplex::LPRowSetRational& inequalities);

    /**
     * \brief Returns all defining inequalities.
     *
     * Returns a const reference to all inequalities added so far.
     */

    inline const soplex::LPRowSetRational& inequalities() const
    {
      return _inequalities;
    }

    /**
     * \brief Returns a dense normal vector of an inequality that defines this face.
     *
     * Returns a const reference to the dense normal vector of an inequality
     * that defines this face. The corresponding right-hand side can be obtained using rhs().
     *
     * \sa sparseNormal()
     */

    inline const soplex::VectorRational& denseNormal()
    {
      ensureSync();
      return _denseNormal;
    }

    /**
     * \brief Returns a sparse normal vector of an inequality that defines this face.
     *
     * Returns a const reference to the sparse normal vector of an inequality that defines this
     * face. The corresponding right-hand side can be obtained using rhs().
     *
     * \sa denseNormal()
     */

    inline const soplex::SVectorRational& sparseNormal()
    {
      ensureSync();
      return _sparseNormal;
    }

    /**
     * \brief Returns the right-hande side of an inequality that defines this face.
     *
     * Returns the right-hande side of an inequality that defines this face. The corresponding
     * normal vector can be obtained using sparseNormal() or denseNormal().
     */

    inline const soplex::Rational& rhs()
    {
      ensureSync();
      return _rhs;
    }

    /**
     * \brief Returns the maximum norm of the stored normal vector that defines this face.
     *
     * Returns the maximum norm of the vector returned by denseNormal() or sparseNormal().
     */

    inline const soplex::Rational& maxNorm()
    {
      ensureSync();
      return _maximumNorm;
    }

    /**
     * \brief Checks whether the given \p point lies in the face.
     *
     * Returns true iff the \p point lies in the face.
     */

    bool containsPoint(const soplex::SVectorRational& point);

    /**
     * \brief Checks whether the given \p direcvtion lies in the recession cone of the face.
     *
     * Returns true iff the \p direction lies in the recession cone of the face.
     */

    bool containsDirection(const soplex::SVectorRational& direction);

  protected:

    /**
     * \brief Ensures that all variables are in sync.
     *
     * Ensures that results of sparseNormal(), denseNormal() rhs() and maxNorm() are in sync with
     * inequalities().
     */

    void ensureSync();

    soplex::LPRowSetRational _inequalities; // Set of inequalities defining this face.
    bool _synced; // Whether \c _inequalities is in sync with the remaining variables.
    soplex::DVectorRational _denseNormal; // Dense normal vector of representing inequality.
    soplex::DSVectorRational _sparseNormal; // Sparse normal vector of representing inequality.
    soplex::Rational _maximumNorm; // Largest absolute number occuring in _normal.
    soplex::Rational _rhs; // Right-hand side of representing inequality.
  };

  /**
   * \brief Results of a call to an oracle.
   *
   * Stores the data that an optimization oracle can return, in particular the returned unbounded
   * directions or the returned points along with their objective values.
   *
   * The caller of the oracle is responsible for freeing the points and directions returned.
   */

  class OracleResult
  {
  public:
    struct Point
    {
      soplex::Rational objectiveValue;
      soplex::DSVectorRational const* point;
      std::size_t index;

      inline bool operator<(const Point& other) const
      {
        return objectiveValue > other.objectiveValue;
      }
    };

    struct Direction
    {
      soplex::DSVectorRational const* direction;
      std::size_t index;
    };

    /**
     * \brief Constructs a result capable of storing oracle answers.
     *
     * Constructs a result capable of storing oracle answers.
     */

    OracleResult();

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    ~OracleResult();

    // TODO: Call it heuristicLevel instead.

    /**
     * \brief Returns the heuristic level of the answer.
     *
     * Returns the heuristic level of the answer. If it is equal to zero and points are returned,
     * then an optimum must be among them.
     */

    inline std::size_t heuristic() const
    {
      return _heuristic;
    }

    /**
     * \brief Returns \c true iff neither points nor directions are returned.
     *
     * Returns \c true iff neither points nor directions are returned. In this case, heuristic()
     * must return 0, since a heuristic must forward a call if nothing was found.
     */

    inline bool isInfeasible()
    {
      return points.empty() && directions.empty();
    }

    /**
     * \brief Returns \c true iff points are returned.
     *
     * Returns \c true iff points (and hence no directions) are returned. If heuristic() returns 0,
     * then an optimum is among them.
     */

    inline bool isFeasible()
    {
      return !points.empty();
    }

    /**
     * \brief Returns \c true iff directions are returned.
     *
     * Returns \c true iff directions (and hence no points) are returned. In this case, heuristic()
     * must return 0, since an exact oracle cannot do better.
     */

    inline bool isUnbounded()
    {
      return !directions.empty();
    }

    /**
     * \brief Resets the status and removes all points and directions.
     *
     * Resets the status and removes all points and directions without freeing them.
     *
     * \param objective Objective vector of the current oracle call.
     */

    void buildStart(const soplex::VectorRational& objective);

    /**
     * \brief Adds the given \p point to the set of \c points.
     *
     * Adds the given \p point to the set of \c points without copying it.
     */

    void buildAddPoint(soplex::DSVectorRational const* point);

    /**
     * \brief Adds the given \p direction to the set of \c directions.
     *
     * Adds the given \p direction to the set of \c directions without copying it.
     */

    void buildAddDirection(soplex::DSVectorRational const* direction);

    /**
     * \brief Finishes construction of oracle answer.
     *
     * Finishes construction of oracle answer.
     *
     * \param heuristic              Sets the heuristic level of the answer.
     * \param computeObjectiveValues Compute the objective values of all points.
     * \param sort                   Sorts the points in descending order of objective value.
     * \param removeDuplicates       Removes duplicate points and directions.
     */

    void buildFinish(std::size_t heuristic, bool computeObjectiveValues, bool sort,
      bool removeDuplicates);

    /**
     * \brief Adds all stored points and directions to the respective containers.
     *
     * Adds all stored points and directions to the respective containers. Updates the indices
     * for the solutions that needed to be added. After this call, none of the solutions must
     * be freed by the user.
     */

    void addToContainers(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions);

    /**
     * \brief Checks the state for consistency.
     *
     * Checks the state for consistency, including result status, ordering of points and
     * presence of duplicates.
     */

    void checkConsistent();

  protected:

    /**
     * \brief Remove duplicate points and directions.
     *
     * Remove duplicate points and directions. If \p abort is \c true, then it aborts as soon as a
     * duplicate is detected.
     */

    bool removeDuplicates(bool abort);

  public:

    /**
     * \brief Array of returned feasible points.
     *
     * Array of returned feasible points. When an oracle call returns, this array must be sorted
     * in descending order of objective values and must not contain duplicate points.
     */

    std::vector<Point> points;

    /**
     * \brief Array of returned unbounded directions.
     *
     * Array of returned unbounded directions. When an oracle call returns, this array must not
     * contain duplicate directions.
     */

    std::vector<Direction> directions;

  protected:
    soplex::VectorRational const* _objective;
    std::size_t _heuristic;
  };

  struct ObjectiveBound
  {
    soplex::Rational value;
    bool strict;

    inline ObjectiveBound() : value(-soplex::infinity), strict(false)
    {

    }

    inline ObjectiveBound(const soplex::Rational& val, bool strct) : value(val), strict(strct)
    {

    }

    inline bool satisfiedBy(const soplex::Rational& other) const
    {
      if (strict)
        return other > value;
      else
        return other >= value;
    }
  };

  // TODO: Change interface such that default heuristicLevel pre- and postprocessing during the call
  // does not need to be implemented for each oracle.

  /**
   * \brief Base class for an optimization oracle.
   *
   * Base class for an optimization oracle for a polyhedron \f$ P \f$ and all its faces. The oracle
   * is either <b>exact</b> or it is associated to another optimization oracle, called the
   * <b>next oracle</b> and called <b>heuristic</b>. The next oracle in turn need not be exact.
   * All oracles reachable by going to the next one are called <b>associated</b>.
   * The \c isExact() method indicates whether the oracle is exact or not and the \c thisHeuristic()
   * method indicates the number of associated oracles (being zero if exact). The goal of this
   * concept is to allow the implementation of quick heuristcs or approximation algorithms
   * that forward a query to the next oracle if its own answer is not satisfactory.
   *
   * An instance has a reference to a \ref Space that manages the ambient space of \f$ P \f$,
   * and that must agree for all associated oracles. An actual implementation of an oracle
   * (or heuristic) must inherit from this class and implement the maximize() and setFace()
   * methods. If there is no direct way to implement the optimization over arbitrary faces,
   * consider inheriting from \ref FaceOracleBase instead.
   *
   * A call to any of the maximize() methods to optimize objective \f$ c \in \mathbb{Q}^n \f$
   * over the current face \f$ F \f$ (\c setFace()) with \c improveValue \f$ \gamma \f$ and
   * requested \c maxHeuristic and \c minHeuristic must obey the following rules, where we denote by
   * \f$ S \subseteq F \f$ and \f$ R \subseteq \text{recc}(P) \f$ the sets of returned \c points
   * and \c directions, respectively. The returned \ref OracleResult is denoted by \c result.
   *
   * \li \c result.heuristic must be at least \c minHeuristic and at most \c maxHeuristic.
   * \li One of \f$ S \f$ and \f$ R \f$ (or both) must be empty.
   * \li If \f$ S = R = \emptyset \f$, then \c result.heuristic must be equal to \c minHeuristic.
   *     In particular, if \c minHeuristic is 0, then \f$P = \emptyset\f$ must hold.
   * \li Every \f$ r \in R\f$ satisfies \f$\left<c,r\right> > 0\f$.
   * \li If \c result.heuristic is equal to zero, then \f$ R \neq \emptyset \f$ holds if and only
   *     if there exists a direction \f$r \in \text{recc}(F)\f$ with \f$\left<c,r\right> > 0\f$
   *     (given \f$P \neq \emptyset\f$).
   * \li If \c result.heuristic is equal to zero and \f$P \neq \emptyset\f$ holds,
   *     then \f$ P \neq \emptyset\f$ or \f$ R \neq \emptyset \f$ must hold.
   * \li If \f$ S \neq \emptyset \f$ and \c result.heuristic is equal to zero, then
   *     \f$ \max\{ \left<c,x\right> \mid x \in P\} = \{ \left<c,s\right> \mid s \in S\} \f$
   *     must hold.
   * \li If \f$ S \neq \emptyset \f$ and \c result.heuristic is positive, then
   *     \f$ \max\{ \left<c,s\right> \mid s \in S\} > \gamma \f$ must hold.
   **/

  class OracleBase
  {
  public:
    // Destructor.

    virtual ~OracleBase();

    /**
     * \brief Returns the name of the oracle.
     *
     * Returns the name of the oracle.
     */

    inline const std::string& name() const
    {
      return _name;
    }

    /**
     * \brief Returns the ambient \c space.
     *
     * Returns a reference to the ambient \c space.
     */

    inline const Space& space() const
    {
      return _space;
    }

    /**
     * \brief Returns the heuristic level of this oracle.
     *
     * Returns the heuristic level of this oracle which is equal to the number of heuristics it is
     * associated to. This number is equal to zero if the oracle is not associated to another
     * oracle. Otherwise, it is equal to the value of thisHeuristic() for the associated oracle
     * plus 1.
     */

    inline std::size_t thisHeuristic() const
    {
      return _thisHeuristic;
    }

    /**
     * \brief Returns the current face.
     *
     * Returns the current face.
     * \sa setFace()
     */

    inline Face* currentFace()
    {
      return _currentFace;
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
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual void maximize(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0) = 0;

    /**
     * \brief Runs the oracle to maximize the dense real \p objective.
     *
     * Runs the optimization oracle to maximize the given dense real \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p improveValue is not exceeded, then the call must be forwarded to the next
     * oracle.
     *
     * This default implementation directly forwards to the next oracle or calls the maximize()
     * method for dense rational objectives.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{F}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual void maximize(OracleResult& result, const soplex::VectorReal& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

    /**
     * \brief Runs the oracle to maximize the sparse rational \p objective.
     *
     * Runs the optimization oracle to maximize the given dense real \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p improveValue is not exceeded, then the call must be forwarded to the next
     * oracle.
     *
     * This default implementation directly forwards to the next oracle or calls the maximize()
     * method for dense rational objectives.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual void maximize(OracleResult& result, const soplex::SVectorRational& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

    /**
     * \brief Runs the oracle to maximize the sparse real \p objective.
     *
     * Runs the optimization oracle to maximize the given dense real \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p improveValue is not exceeded, then the call must be forwarded to the next
     * oracle.
     *
     * This default implementation directly forwards to the next oracle or calls the maximize()
     * method for dense rational objectives.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{F}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual void maximize(OracleResult& result, const soplex::SVectorReal& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

  protected:

    /**
     * \brief Constructs an oracle with given \p name in given \p space.
     *
     * Constructs an oracle with given \p name in given \p space.
     */

    OracleBase(const std::string& name, const Space& space);

    /**
     * \brief Constructs a heuristic with given \p name associated to \p nextOracle.
     *
     * Constructs a heuristic oracle with given \p name that is associated to \p nextOracle. The
     * ambient space is equal to that of \p nextOracle.
     */

    OracleBase(const std::string& name, OracleBase* nextOracle);

    /**
     * \brief Initializes datastructures that require the space.
     *
     * Initializes datastructures that require the space. This method should be called at the end
     * of the constructor of a subclass.
     */

    void initializedSpace();

  private:
    /**
     * \brief Default constructor is disabled.
     *
     * Default constructor is disabled.
     */

    OracleBase();

  protected:
    std::string _name; // Name of the oracle.
    const Space& _space; // Ambient space of the oracle.
    OracleBase* _nextOracle; // Next optimization oracle (or NULL if exact).
    std::size_t _thisHeuristic; // Number of associated oracles.
    Face* _currentFace; // Currently active face.
    soplex::DVectorRational _tempObjective; // Dense rational version of the current objective.
  };

  /**
   * \brief Base class for an optimization oracle without direct handling of faces.
   *
   * This class adds the ability to optimize over faces of \f$ P \f$ by tilting objective vectors.
   * If \f$ \left<a,x\right> \leq \beta \f$ defines the face, then this implementation
   * attempts to optimize a given objective vector \f$ \left<c,x\right> \f$ over \f$ F \f$ by
   * optimizing
   * \f$ \left<c + M \dfrac{ ||a||_{\infty} }{ ||c||_{\infty} } a,x\right> \f$ over \f$ P \f$
   * for sufficiently large \f$ M \in \mathbb{Q} \f$, increasing \f$ M \f$ as long as necessary.
   */

  class FaceOracleBase: public OracleBase
  {
  public:
    /**
     * Destructor.
     */

    virtual ~FaceOracleBase();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation also calls OracleBase::setFace(), which in turn calls
     * setFace() for the next oracle.
     */

    virtual void setFace(Face* newFace = NULL);

    /**
     * \brief Runs the oracle to maximize the dense rational \p objective.
     *
     * Runs the face optimization oracle to maximize the given dense rational \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p improveValue is not exceeded, then the call must be forwarded to the next
     * oracle.
     *
     * This method directly forwards to the next oracle or calls the unrestrictedMaximize()
     * method with a modified objective vector.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual void maximize(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveValue = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

  protected:
    /**
     * \brief Constructs an oracle with given \p name in given \p space.
     *
     * Constructs an exact face optimization oracle with given \p name in given \p space.
     *
     * \param name               Name of the oracle.
     * \param space              Ambient space.
     * \param numBlindIterations Number of times, \f$ M \f$ is increased, before testing whether
     *                           \f$ F = \emptyset \f$ holds.
     * \param initialM           Initial value of \f$ M \f$.
     */

    FaceOracleBase(const std::string& name, const Space& space,
      std::size_t numBlindIterations = 2, double initialM = 16);

    /**
     * \brief Constructs a heuristic with given \p name associated to \p nextOracle.
     *
     * Constructs a heuristic optimization oracle with given \p name that is associated to
     * \p nextOracle. The ambient space is equal to that of \p nextOracle.
     *
     * \param name               Name of the oracle.
     * \param nextOracle         Next oracle to forward calls to.
     * \param numBlindIterations Number of times, \f$ M \f$ is increased, before testing whether
     *                           \f$ F = \emptyset \f$ holds.
     * \param initialM           Initial value of \f$ M \f$.
     */

    FaceOracleBase(const std::string& name, OracleBase* nextOracle,
      std::size_t numBlindIterations = 2, double initialM = 16);

    /**
     * \brief Initializes datastructures that require the space.
     *
     * Initializes datastructures that require the space. This method should be called at the end
     * of the constructor of a subclass.
     */

    void initializedSpace();

    /**
     * \brief Implementation of the oracle.
     *
     * Implements the optimization oracle to maximize the given dense rational \p objective over
     * the whole polyhedron \f$ P \f$ (in contrast to the maximize() methods) and returns \p result.
     * If the objective value requested by \p improveValue is not exceeded, then the call must be
     * forwarded to the next oracle using \p originalObjective and \p originalImproveValue to allow
     * the next oracle to use (potentially more efficient) handling of face optimization.
     *
     * \param result               After the call, contains the oracle's answer.
     * \param objective            Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized over
     *                             \f$ P \f$.
     * \param objectiveBound       Objective value \f$ \gamma \in \mathbb{Q} \f$ that should be
     *                             exceeded.
     * \param originalObjective    Objective vector \f$ c' \in \mathbb{Q}^n \f$ to be maximized
     *                             over \f$ F \f$.
     * \param originalImproveValue Objective value \f$ \gamma' \in \mathbb{Q} \f$ that should be
     *                             exceeded when optimizing \f$ c' \f$ over \f$ F \f$.
     * \param maxHeuristic         Requested maximum heuristic level.
     * \param minHeuristic         Requested minimum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual void unrestrictedMaximize(OracleResult& result,
      const soplex::VectorRational& objective, const ObjectiveBound& improveValue,
      const soplex::VectorRational& originalObjective, const ObjectiveBound& orginalObjectiveBound,
      std::size_t maxHeuristic, std::size_t minHeuristic) = 0;

  protected:
    std::size_t _numBlindIterations; // Number of iterations assuming that the face is non-empty.
    soplex::Rational _M; // Dynamic constant for modifying objective.
    soplex::Rational _factor; // Cachable part of the constant for modifying objective.
    soplex::DVectorRational _modifiedObjective; // Modified objective.
  };

} /* namespace ipo */

#endif /* IPO_ORACLE_H_ */
