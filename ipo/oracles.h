#ifndef IPO_ORACLES_H_
#define IPO_ORACLES_H_

#include <vector>
#include <limits>

#include "common.h"
#include "rows.h"
#include "space.h"
#include "vectors.h"

namespace ipo {

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

    inline const DenseVector& denseNormal()
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

    inline const SparseVector& sparseNormal()
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

    inline const Rational& rhs()
    {
      ensureSync();
      return _rhs;
    }

    /**
     * \brief Returns the maximum norm of the stored normal vector that defines this face.
     *
     * Returns the maximum norm of the vector returned by denseNormal() or sparseNormal().
     */

    inline const Rational& maxNorm()
    {
      ensureSync();
      return _maximumNorm;
    }

    /**
     * \brief Checks whether the given \p point lies in the face.
     *
     * Returns true iff the \p point lies in the face.
     */

    bool containsPoint(const SparseVector& point);

    /**
     * \brief Checks whether the given \p direcvtion lies in the recession cone of the face.
     *
     * Returns true iff the \p direction lies in the recession cone of the face.
     */

    bool containsDirection(const SparseVector& direction);

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
    DenseVector _denseNormal; // Dense normal vector of representing inequality.
    SparseVector _sparseNormal; // Sparse normal vector of representing inequality.
    Rational _maximumNorm; // Largest absolute number occuring in _normal.
    Rational _rhs; // Right-hand side of representing inequality.
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
      SparseVector vector;

      Point(SparseVector& vector);
      Point(SparseVector& vector, const soplex::Rational& value);

      inline bool operator<(const Point& other) const
      {
        return objectiveValue > other.objectiveValue;
      }
    };

    struct Direction
    {
      SparseVector vector;

      Direction(SparseVector& vector);
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

    /**
     * \brief Returns the heuristic level of the answer.
     *
     * Returns the heuristic level of the answer. If it is equal to zero and points are returned,
     * then an optimum must be among them.
     */

    inline std::size_t heuristicLevel() const
    {
      return _heuristicLevel;
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
     * \brief Checks the state for consistency.
     *
     * Checks the state for consistency, including result status, ordering of points and
     * presence of duplicates.
     */

    void checkConsistent();
    
    /**
     * \brief Computes all points' objective values.
     * 
     * Computes all points' objective values.
     */
    
    void computeMissingObjectiveValues();

  protected:

    /**
     * \brief Remove duplicate points and directions.
     *
     * Remove duplicate points and directions.
     */

    void removeDuplicates();
    
    friend class OracleBase;

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
    DenseVector const* _objective;
    std::size_t _heuristicLevel;
  };

  struct ObjectiveBound
  {
    Rational value;
    bool strict;

    inline ObjectiveBound() : value(minusInfinity), strict(false)
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
   * \li \c result.heuristic must be at least \c minHeuristic and at most \c maxHeuristic, or equal to 0.
   * \li One of \f$ S \f$ and \f$ R \f$ (or both) must be empty.
   * \li If \f$ S = R = \emptyset \f$, then \c result.heuristic must be equal to 0 and \f$F = \emptyset\f$ must hold.
   * \li Every \f$ r \in R\f$ satisfies \f$\left<c,r\right> > 0\f$.
   * \li If \c result.heuristic is equal to 0, then \f$ R \neq \emptyset \f$ holds if and only
   *     if there exists a direction \f$r \in \text{recc}(F)\f$ with \f$\left<c,r\right> > 0\f$
   *     (given \f$F \neq \emptyset\f$).
   * \li If \c result.heuristic is equal to 0 and \f$F \neq \emptyset\f$ holds,
   *     then \f$S \neq \emptyset\f$ or \f$ R \neq \emptyset \f$ must hold.
   * \li If \f$ S \neq \emptyset \f$ and \c result.heuristic is equal to zero, then
   *     \f$ \max\{ \left<c,x\right> \mid x \in F\} = \{ \left<c,s\right> \mid s \in S\} \f$
   *     must hold.
   * \li If \f$ S \neq \emptyset \f$ and \c result.heuristic is positive, then
   *     \f$ \max\{ \left<c,s\right> \mid s \in S\} > \gamma \f$ must hold or \f$S\f$ must contain at most one element
   *     (to certify feasibility).
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

    inline std::size_t heuristicLevel() const
    {
      return _heuristicLevel;
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
     * Runs the optimization oracle to maximize the given dense rational \p objective over the current face \f$ F \f$ (see 
     * setFace()) and returns \p result. If \p maxHeuristic is less than thisHeuristic() or if the objective value requested by 
     * \p objectiveBound is not exceeded, then the call must be forwarded to the next oracle.
     * 
     * This implementation initializes the \p result data structure, calls the maximizeController() method (with a dense 
     * objective vector), and finally postprocesses the result object.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param minHeuristic   Requested minimum heuristic level.
     * \param maxHeuristic   Requested maximum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    void maximize(OracleResult& result, const DenseVector& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(), std::size_t minHeuristic = 0,
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max());

    /**
     * \brief Runs the oracle to maximize the sparse rational \p objective.
     *
     * Runs the optimization oracle to maximize the given sparse rational \p objective over the current face \f$ F \f$ (see 
     * setFace()) and returns \p result. If \p maxHeuristic is less than thisHeuristic() or if the objective value requested by 
     * \p objectiveBound is not exceeded, then the call must be forwarded to the next oracle.
     * 
     * This implementation initializes the \p result data structure, calls the maximizeController() method (with a dense 
     * objective vector), and finally postprocesses the result object.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param minHeuristic   Requested minimum heuristic level.
     * \param maxHeuristic   Requested maximum heuristic level.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    void maximize(OracleResult& result, const SparseVector& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(), std::size_t minHeuristic = 0,
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max());

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
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual std::size_t maximizeController(OracleResult& result, const DenseVector& objective,
      const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic, bool& sort, bool& checkDups);

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
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual std::size_t maximizeImplementation(OracleResult& result, const DenseVector& objective,
      const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups) = 0;

    /**
     * \brief Returns the current face.
     *
     * Returns the current face. \sa setFace()
     */

    virtual Face* currentFace();

  private:
    /**
     * \brief Default constructor is disabled.
     *
     * Default constructor is disabled.
     */

    OracleBase();
    
    Face* _currentFace; // Currently active face.

  protected:
    friend class FaceOracleBase;
    
    std::string _name; // Name of the oracle.
    const Space& _space; // Ambient space of the oracle.
    OracleBase* _nextOracle; // Next optimization oracle (or NULL if exact).
    std::size_t _heuristicLevel; // Number of associated oracles.
    DenseVector _tempObjective; // Dense rational versi::DVectorRational on of the current objective.
  };

  /**
   * \brief Base class for an optimization oracle without direct handling of faces.
   *
   * This class adds the ability to optimize over faces of \f$ P \f$ by tilting objective vectors.
   * If \f$ \left<a,x\right> \leq \beta \f$ defines the face, then this implementation attempts to optimize a given objective 
   * vector \f$ \left<c,x\right> \f$ over \f$ F \f$ by optimizing
   * \f$ \left<c + M \dfrac{ ||a||_{\infty} }{ ||c||_{\infty} } a,x\right> \f$ over \f$ P \f$
   * for sufficiently large \f$ M \in \mathbb{Q} \f$, doubling \f$ M \f$ a certain nummber of times before giving up.
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
    
  protected:

    /**
     * \brief Constructs an oracle with given \p name in given \p space.
     *
     * Constructs an exact face optimization oracle with given \p name in given \p space.
     *
     * \param name                      Name of the oracle.
     * \param space                     Ambient space.
     * \param maxInfeasibleIterations   Maximum number of iterations before (heuristically) checking if the face is empty.
     * \param initialM                  Initial value of \f$ M \f$.
     */

    FaceOracleBase(const std::string& name, const Space& space, std::size_t maxInfeasibleIterations = 4,
      double initialM = 16);

    /**
     * \brief Constructs a heuristic with given \p name associated to \p nextOracle.
     *
     * Constructs a heuristic optimization oracle with given \p name that is associated to
     * \p nextOracle. The ambient space is equal to that of \p nextOracle.
     *
     * \param name               Name of the oracle.
     * \param nextOracle         Next oracle to forward calls to.
     * \param maxInfeasibleIterations   Maximum number of iterations before (heuristically) checking if the face is empty.
     * \param initialM           Initial value of \f$ M \f$.
     */

    FaceOracleBase(const std::string& name, OracleBase* nextOracle, std::size_t maxInfeasibleIterations = 4,
      double initialM = 16);

    /**
     * \brief Initializes datastructures that require the space.
     *
     * Initializes datastructures that require the space. This method should be called at the end
     * of the constructor of a subclass.
     */

    void initializedSpace();
    
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
     * The maximizeImplementation() function is called repeatedly with titled objective vectors (see \ref FaceOracleBase).
     */

    virtual std::size_t maximizeController(OracleResult& result, const DenseVector& objective,
      const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic, bool& sort, bool& checkDups);

    /**
     * \brief Returns the current face.
     *
     * Returns the current face. \sa setFace(). Since handling of faces is done via modifications of the objective vector, this 
     * implementation always returns NULL in order to remove face-handling from maximizeImplementation()'s requirements. To 
     * actually access the current face, call OracleBase::currentFace().
     */

    virtual Face* currentFace();
    
      
  protected:
    std::size_t _maxInfeasibleIterations; // Maximum number of iterations before (heuristically) checking if the face is empty.
    soplex::Rational _M; // Dynamic constant for modifying objective.
    soplex::Rational _factor; // Cachable part of the constant for modifying objective.
    soplex::DVectorRational _modifiedObjective; // Modified objective.
  };

} /* namespace ipo */

#endif /* IPO_ORACLE_H_ */
