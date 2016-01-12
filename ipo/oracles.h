#ifndef IPO_ORACLES_H_
#define IPO_ORACLES_H_

#include <vector>
#include <limits>

#include "ipo.h"
#include "spx_gmp.h"
#include "rows.h"

namespace ipo {

  /**
   * \brief Defines a face of a polyhedron by a set of inequalities.
   * 
   * Defines a face \f$F\f$ of a polyhedron \f$P\f$ by a set of inequalities.
   * It is used to create an optimization oracle for \f$F\f$.
   * If the optimization oracle class for \f$P\f$ inherits from
   * \ref FaceOptimizationOracleBase, the face can be controlled directly.
   * In any case one can construct an instance of \ref \FaceOptimizationOracle
   * that calls the optimization oracle for \f$P\f$ (maybe multiple times per call).
   **/

  class Face
  {
  public:
    /**
     * Creates the trivial face for ambient dimension \c numVariables defined by
     * \f$\left<\mathbbm{O},x\right> \leq 0\f$.
     * Use \ref add() methods to add further inequalities.
     **/

    Face(std::size_t numVariables);

    /**
     * Creates the face for ambient dimension \c numVariables defined by
     * \c inequality.
     * Use \ref add() methods to add further inequalities,
     * that is, to intersect with other faces.
     **/

    Face(std::size_t numVariables, const soplex::LPRowRational& inequality);

    /**
     * Creates the intersection of the faces defined by
     * \c inequalities.
     * The ambient dimension is \c numVariables.
     * Use \ref add() methods to add further inequalities,
     * that is, to intersect with other faces.
     **/

    Face(std::size_t numVariables, const soplex::LPRowSetRational& inequalities);

    // Destructor.

    virtual ~Face();

    /**
     * Adds \c inequality to the inequality currently defining this face.
     * Geometrically this means to intersect the current face with the
     * one defined by \c inequality.
     */

    void add(const soplex::LPRowRational& inequality);

    /**
     * Adds all \c inequalities to the inequality currently defining this face.
     * Geometrically this means to intersect the current face with all
     * faces defined by \c inequalities.
     */

    void add(const soplex::LPRowSetRational& inequalities);

    /**
     * Returns a const reference to all inequalities added so far.
     */

    inline const soplex::LPRowSetRational& inequalities() const
    {
      return _inequalities;
    }

    /**
     * Returns a const reference to the normal vector of an inequality
     * that defines this face.
     *
     * \sa rhs()
     */

    inline const soplex::SVectorRational& normal()
    {
      ensureSync();
      return _normal;
    }

    /**
     * Returns the right-hande side of an inequality
     * that defines this face.
     *
     * \sa normal()
     */

    inline const soplex::Rational& rhs()
    {
      ensureSync();
      return _rhs;
    }

    /**
     * Returns the largest absolute value of a coefficient of
     * the vector returned by \ref normal().
     */

    inline const soplex::Rational& largestAbsCoefficient()
    {
      ensureSync();
      return _largestAbsCoefficient;
    }

  protected:

    // Ensures that results of normal(), rhs() and
    // largestAbsCoefficient() are in sync with inequalities().

    void ensureSync();

    soplex::LPRowSetRational _inequalities; // Set of inequalities defining this face.
    bool _synced; // Whether _worker, _normal, _largestAbsCoefficient and _rhs are in sync with _inequalities.
    soplex::DVectorRational _worker; // Dense normal vector of representing inequality.
    soplex::DSVectorRational _normal; // Sparse normal vector of representing inequality.
    soplex::Rational _largestAbsCoefficient; // Largest absolute number occuring in _normal.
    soplex::Rational _rhs; // Right-hand side of representing inequality.
  };

  /**
   * \brief Results of a call to an optimization oracle.
   *
   * Stores the data that an oracle can return,
   * in particular the returned \c points (with their \c objectives) or unbounded \c directions.
   * The latter objects must be freed by the user.
   * It also contains information about the optimality of the result.
   */

  struct OptimizationResult
  {
    /**
     *
     * \brief If \c true, then optimality is guaranteed.
     * If points are returned, then it is \c true iff an optimum is among the returned points.
     * If neither points nor directions are returned,
     * then it is \c true iff the oracle guarantees emptyness of \f$P\f$.
     */

    bool optimal;

    /**
     * Index of the best returned point.
     */

    std::size_t bestIndex;

    /**
     * Objective value of best returned point.
     */

    soplex::Rational bestValue;

    /**
     * Points returned by the optimization oracle.
     */

    std::vector<soplex::DSVectorRational*> points;

    /**
     * Objective values corresponding to \points.
     */

    std::vector<soplex::Rational> objectives;

    /**
     * Unbounded directions returned.
     */

    std::vector<soplex::DSVectorRational*> directions;

    /**
     * Returns \c true iff the optimization oracle claims \f$P = \emptyset\f$.
     */

    inline bool isInfeasible() const
    {
      return objectives.empty() && directions.empty();
    }

    /**
     * Returns \c true iff the optimization oracle returned
     * unbounded \c directions.
     */

    inline bool isUnbounded() const
    {
      return !directions.empty();
    }

    /**
     * Returns \c true iff the optimization oracle found
     * no unbounded directions, but \points.
     */

    inline bool isFeasible() const
    {
      return !objectives.empty();
    }

    /**
     * \brief Initializes the structure during oracle call.
     *
     * Initializes the structure to return results in ambient dimension \numVariables.
     * This method should be called in an oracle implementation before
     * actually adding \c points or \c directions.
     */

    void reset(std::size_t numVariables);

    /**
     * \brief Creates a new point and adds it to the set of \c points.
     *
     * Creates a new point and adds it to the set of \c points.
     *
     * \returns Reference to sparse vector representing this point.
     */

    soplex::DSVectorRational& newPoint();

    /**
     * \brief Marks result as feasible during oracle call.
     *
     * Marks result as feasible, computes all objective values,
     * setting index and value of best returned point.
     * This method should be called in an oracle implementation after
     * adding all \c points or \c directions.
     *
     * \param objective Original objective vector passed to the oracle.
     */

    void setFeasible(const soplex::VectorRational& objective);

    /**
     * \brief Updates index and value of best returned point.
     *
     * Updates index and value of best returned point.
     *
     * \param value New best objective value.
     * \param index New index of best returned point.
     */

    void setBest(const soplex::Rational& value, std::size_t index = std::numeric_limits<std::size_t>::max());

    /**
     * \brief Marks result as infeasible during oracle call.
     *
     * Marks result as infeasible.
     * This method should be called in an oracle implementation.
     */

    void setInfeasible();

    /**
     * \brief Marks result as unbounded during oracle call.
     *
     * Marks result as unbounded.
     * This method should be called in an oracle implementation.
     */

    void setUnbounded();

    /**
     * \brief Inspects all points and rays for duplicates and removes them.
     *
     * Inspects all points and rays for duplicates and removes them.
     */

    void filterDuplicates();

#ifdef IPO_DEBUG
  public:
    void checkConsistent() const;
    bool hasDuplicates() const;
#endif

  protected:

    /*
     * Removes duplicates from a given set of vectors.
     */

    void filterDuplicates(std::vector<soplex::DSVectorRational*>& vectors);

  protected:
    std::size_t _numVariables; /// Ambient dimension.

  };

  /**
   * \brief Base class for an oracle.
   *
   * Base class that every optimization oracle must inherit from,
   * implementing the \c maximize() and \c improve() methods or the \c run() method.
   * This class also manages the names of all variables,
   * and hence the ambient dimension.
   *
   * Every optimization oracle implementation for a polyhedron \f$P \subseteq \mathbb{R}^n\f$
   * must obey the following rules when called with an objective vector \f$c \in \mathbb{R}^n\f$:
   *
   * \li Return \c points \f$S \subseteq P\f$
   *     or \c directions \f$R \subseteq \text{recc}(P)\f$, but not both.
   *     If \f$P = \emptyset \f$ holds, it must not return anything.
   * \li Every returned direction \f$r \in R\f$ must satisfy \f$\left<c,r\right> > 0\f$.
   * \li If \c forceOptimal is \c true and there exists a direction \f$r \in \text{recc}(P)\f$
   *     with \f$\left<c,r\right> > 0\f$
   *     (given \f$P \neq \emptyset\f$), then it must return such a direction.
   * \li If \c forceOptimal is \c true and \f$P \neq \emptyset\f$ holds,
   *     then it <b>must</b> return \c points or \c directions.
   * \li If \c forceOptimal is \c true and it returns \c points,
   *     then \f$ \max\{ \left<c,x\right> \mid x \in P\} = \{ \left<c,s\right> \mid s \in S\} \f$
   *     must hold.
   **/

  class OptimizationOracleBase
  {
  public:
    // Destructor.

    virtual ~OptimizationOracleBase();

    /**
     * \brief Returns the ambient dimension of \f$ P \f$.
     *
     * Returns the ambient dimension of \f$ P \f$.
     */

    inline std::size_t numVariables() const
    {
#ifdef IPO_DEBUG
      if (!_initialized)
        throw std::runtime_error("OptimizationOracleBase not initialized!");
#endif
      return _variableNames.size();
    }

    /**
     * \brief Returns the name of the variable of index \c var.
     *
     * Returns the name of the variable of index \c var,
     * where \c var must nonnegative and less than \ref numVariables().
     */

    inline const std::string& variableName(std::size_t var) const
    {
#ifdef IPO_DEBUG
      if (!_initialized)
        throw std::runtime_error("OptimizationOracleBase not initialized!");
#endif
      return _variableNames[var];
    }

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
     * \brief Prints the \c row (inequality / equation) to \c stream using the variable names of the oracle.
     *
     * Prints the \c row (inequality / equation) to \c stream using variable names of the oracle.
     * Does not emit a newline character.
     */

    void printRow(std::ostream& stream, const soplex::LPRowRational& row) const;

    /**
     * \brief Prints the row (inequality / equation) \c index from the set of \c rows
     * to \c stream using the variable names of the oracle.
     *
     * Prints the row (inequality / equation) \c index from the set of \c rows
     * to \c stream using the variable names of the oracle.
     * Does not emit a newline character.
     */

    void printRow(std::ostream& stream, const soplex::LPRowSetRational& rows, std::size_t index) const;

    /**
     * \brief Prints the \c rows (inequality / equation)
     * to \c stream using the variable names of the oracle.
     *
     * Prints the \c rows (inequality / equation)
     * to \c stream using the variable names of the oracle.
     * Emit a newline character after each row.
     */

    void printRows(std::ostream& stream, const soplex::LPRowSetRational& rows) const;

    /**
     * \brief Prints \c vector to \c stream using the variable names of the oracle.
     *
     * Prints \c vector to \c stream using the variable names of the oracle.
     * Does not emit a newline character.
     */

    void printVector(std::ostream& stream, const soplex::SVectorRational* vector) const;

    /**
     * \brief Runs the optimization oracle for the given dense rational \c objective, returning \c result.
     *
     * Runs the optimization oracle for the given dense rational \c objective, returning \c result.
     * The default implementation calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void maximize(OptimizationResult& result, const soplex::VectorRational& objective,
        bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given dense floating-point \c objective, returning \c result.
     *
     * Runs the optimization oracle for the given dense floating-point \c objective, returning \c result.
     * The default implementation converts it to a rational objective and calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void maximize(OptimizationResult& result, const soplex::VectorReal& objective, bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given sparse rational \c objective, returning \c result.
     *
     * Runs the optimization oracle for the given sparse rational \c objective, returning \c result.
     * The default implementation converts it to a dense objective and calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void maximize(OptimizationResult& result, const soplex::SVectorRational& objective,
        bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given sparse floating-point \c objective, returning \c result.
     *
     * Runs the optimization oracle for the given sparse floating-point \c objective, returning \c result.
     * The default implementation converts it to a dense rational objective and calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void maximize(OptimizationResult& result, const soplex::SVectorReal& objective, bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given dense rational \c objective
     * in order to obtain a better objective than \c value, returning \c result.
     *
     * Runs the optimization oracle for the given dense rational \c objective, returning \c result.
     * Even if \c forceOptimal is \c true, it may abort with non-optimal points
     * if it found one whose objective is greater than \c value.
     * The default implementation calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void improve(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational& value, bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given dense floating-point \c objective
     * in order to obtain a better objective than \c value, returning \c result.
     *
     * Runs the optimization oracle for the given dense floating-point \c objective, returning \c result.
     * Even if \c forceOptimal is \c true, it may abort with non-optimal points
     * if it found one whose objective is greater than \c value.
     * The default implementation converts it to a rational objective and calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void improve(OptimizationResult& result, const soplex::VectorReal& objective, const soplex::Rational& value,
        bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given sparse rational \c objective
     * in order to obtain a better objective than \c value, returning \c result.
     *
     * Runs the optimization oracle for the given sparse rational \c objective, returning \c result.
     * Even if \c forceOptimal is \c true, it may abort with non-optimal points
     * if it found one whose objective is greater than \c value.
     * The default implementation converts it to a sparse objective and calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void improve(OptimizationResult& result, const soplex::SVectorRational& objective,
        const soplex::Rational& value, bool forceOptimal = true);

    /**
     * \brief Runs the optimization oracle for the given sparse floating-point \c objective
     * in order to obtain a better objective than \c value, returning \c result.
     *
     * Runs the optimization oracle for the given sparse floating-point \c objective, returning \c result.
     * Even if \c forceOptimal is \c true, it may abort with non-optimal points
     * if it found one whose objective is greater than \c value.
     * The default implementation converts it to a sparse rational objective and calls the \ref run() method.
     *
     * \params forceOptimal Controls optimality requirements
     *         (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void improve(OptimizationResult& result, const soplex::SVectorReal& objective,
        const soplex::Rational& value, bool forceOptimal = true);

  protected:
    /**
     * \brief Main method for an actual implementation of an optimization oracle.
     *
     * Main method for an actual implementation of an optimization oracle.
     * Runs the optimization oracle for the given dense rational \c objective, returning \c result.
     * If \c improveValue is not \c NULL and \c forceOptimal is \c true,
     * then it may abort with non-optimal points
     * if it found one whose objective is greater than \c *improveValue.
     *
     * \params improveValue
     *   If not NULL, allows early termination if objective \c *improveValue is exceeded.
     * \params forceOptimal
     *   Controls optimality requirements (see Detailed Description of \ref OptimizationOracleBase).
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    /**
     * \brief Constructs an optimization oracle with given \c name.
     *
     * Constructs an optimization oracle with given \c name.
     * An actual implementation must call one of the \ref initialize() methods
     * from its own constructor.
     */

    OptimizationOracleBase(const std::string& name);

    /**
     * \brief Initializes the oracle with the given set of variables.
     *
     * Initializes the oracle with the given set of variables.
     * This method should be called from the constructor of a subclass.
     *
     * \param variableNames Vector of variable names.
     */

    void initialize(const std::vector<std::string>& variableNames);

    /**
     * \brief Initializes the oracle using the variables of another oracle.
     *
     * Initializes the oracle using the variables of another \c oracle.
     * This method should be called from the constructor of a subclass.
     */

    void initialize(const OptimizationOracleBase* oracle);

  private:
    /*
     * This should not be used.
     */

    OptimizationOracleBase();

    /*
     * \brief Implementation of row printing.
     *
     * Implementation of row printing.
     */

    void printRow(std::ostream& stream, const soplex::Rational* lhs, const soplex::Rational* rhs,
        const soplex::SVectorRational& vector) const;

    std::string _name; // Name of the oracle.
    std::vector<std::string> _variableNames; // Variables names.
    soplex::DVectorRational _objective; // Variable holding the dense rational version of the current objective.
#ifdef IPO_DEBUG
    bool _initialized;
#endif
  };

  /**
   * \brief An oracle that calls a given first oracle, and, if not satisfied, a second one.
   *
   * An optimization oracle that calls a given first oracle, and, if not satisfied, a second one.
   */

  class ChainedOptimizationOracle: public OptimizationOracleBase
  {
  public:

    /**
     * \brief Constructor.
     *
     * Constructor.
     *
     * \param first Oracle to be called first.
     * \param second Oracle to be called if the answer of the first is not satisfactory.
     */

    ChainedOptimizationOracle(OptimizationOracleBase* first, OptimizationOracleBase* second);

    /**
     * Destructor.
     */

    virtual ~ChainedOptimizationOracle();

  protected:
    /**
     * \brief Implementation of the oracle.
     *
     * Implementation of the oracle.
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

  protected:
    OptimizationOracleBase* _first; // Pointer to the first oracle.
    OptimizationOracleBase* _second; // Pointer to the second oracle.
  };

  /**
   * \brief Base class for an oracle that can optimize over arbitrary faces.
   *
   * An optimization oracle for a polyhedron \f$ P \f$
   * that can optimize over arbitrary faces.
   * The current face is controlled using the \ref setFace() methods.
   */

  class FaceOptimizationOracleBase: public OptimizationOracleBase
  {
  public:
    /**
     * Destructor.
     */

    virtual ~FaceOptimizationOracleBase();

    /**
     * Sets the face over which the oracle should optimize.
     *
     * \param face Pointer to a \ref Face structure, or \c NULL for the whole polyhedron.
     */

    virtual Face* setFace(Face* face = NULL);

  protected:
    /**
     * \brief Constructs a face optimization oracle with given \c name.
     *
     * Constructs a face optimization oracle with given \c name.
     * An actual implementation must call one of the \ref initialize() methods
     * from its own constructor.
     */

    FaceOptimizationOracleBase(const std::string& name);

    /**
     * \brief Method that is called when a new face is activated.
     *
     * Method that is called when a new \c face is activated.
     * If a (non-trivial) face was active before,
     * it is ensured that \ref faceDisabled() is called before.
     */

    virtual void faceEnabled(Face* face) = 0;

    /**
     * \brief Method that is called when a face is deactivated.
     *
     * Method that is called when a face is deactivated.
     */

    virtual void faceDisabled(Face* face) = 0;

  protected:
    Face* _face; // Pointer to currently active face.
  };

  /**
   * \brief Defines an affine projection.
   *
   * Defines an affine projection \f$ \pi(x) = Ax + b \f$
   * including variable names of the image space.
   * Its main purpose is to create a \ref ProjectedOptimizationOracle
   * from another optimization oracle.
   */

  class Projection
  {
  public:
    /**
     * \brief Construct a projection without image space variables.
     *
     * Construct the projection without image space variables.
     * Use \c addVariable() methods to add some.
     */

    Projection();

    /**
     * \brief Construct an orthogonal projection for a given oracle.
     *
     * Construct an orthogonal projection of the ambient space of a given \c oracle
     * onto a subset of the variables.
     *
     * \param oracle
     *   Given oracle, used to define the source space and the variable names.
     * \param variableSubset
     *   Indices of the subset of variables to project onto.
     */

    Projection(const OptimizationOracleBase* oracle, const std::vector<std::size_t>& variableSubset);

    /**
     * \brief Destructor.
     */

    virtual ~Projection();

    /**
     * \brief Add a variable to the image space.
     *
     * Add a variable to the image space by
     * adding a given row to the projection matrix and to the constant term.
     *
     * \param name
     *   Name of the new variable.
     * \param variableMap
     *   Sparse rational row of the projection matrix \f$ A \f$.
     * \param shift
     *   Rational entry for shift vector \f$ b \f$.
     */

    void addVariable(const std::string& name, const soplex::SVectorRational& variableMap,
        const soplex::Rational& shift = soplex::Rational(0));

    /**
     * \brief Copy a source variable to the image space.
     *
     * Copy a source variable to the image space.
     *
     * \param oracle
     *   Oracle to copy the name from.
     * \param originalVariable.
     *   Index of the original variable.
     * \param shift
     *   Rational entry for shift vector \f$ b \f$.
     */

    void addVariable(const OptimizationOracleBase* oracle, std::size_t originalVariable, const soplex::Rational& shift =
        soplex::Rational(0));

    /**
     * \brief Project a point.
     *
     * Project a \c point.
     */

    void projectPoint(const soplex::DVectorRational& point, soplex::DSVectorRational& image) const;

    /**
     * \brief Returns the number of variables of the image space.
     *
     * Returns the number of variables of the image space.
     */

    inline std::size_t numVariables() const
    {
      return _map.size();
    }

    /**
     * \brief Returns a const reference to the vector of variable names.
     *
     * Returns a const reference to the vector of variable names.
     */

    inline const std::vector<std::string>& names() const
    {
      return _names;
    }

    /**
     * \brief Returns a row of the projection matrix.
     *
     * Returns a row of the projection matrix \f$ A \f$ for a given \c variable
     * as a const reference to a sparse rational vector.
     */

    inline const soplex::SVectorRational& map(std::size_t variable) const
    {
      return _map[variable];
    }

    /**
     * \brief Returns an entry of the projection vector.
     *
     * Returns an entry of the projection vector \f$ b \f$ for a given \c variable.
     */

    inline const soplex::Rational& shift(std::size_t var) const
    {
      return _shift[var];
    }

  protected:
    std::vector<std::string> _names; // Image space variable names.
    std::vector<soplex::DSVectorRational> _map; // Rows of the projection matrix.
    std::vector<soplex::Rational> _shift; // Entries of the projection vector.
  };

  /**
   * \brief An oracle for the projection of a polyhedron defined by another oracle.
   *
   * Defines an optimization oracle for the projection \f$ \pi(P) \f$
   * of a polyhedron \f$ P \f$ defined by another oracle
   * for given projection map \f$ \pi \f$
   *
   * \sa Projection
   */

  class ProjectedOptimizationOracle: public OptimizationOracleBase
  {
  public:
    /**
     * \brief Constructor for given projection and source oracle.
     *
     * Constructor for given \c projection and source \c oracle.
     *
     * \param name
     *   Name of the new oracle.
     * \param projection
     *   Projection map.
     * \param oracle
     *   Source oracle.
     */

    ProjectedOptimizationOracle(const std::string& name, const Projection& projection, OptimizationOracleBase* oracle);

    /**
     * \brief Destructor.
     */

    virtual ~ProjectedOptimizationOracle();

  protected:

    /**
     * \brief Implementation of the oracle.
     *
     * Implementation of the oracle.
     * For a projection \f$ \pi(x) = Ax + b \f$ it
     * calls the oracle with objective \f$ A^{\intercal} c \f$
     * where \f$ c \f$ is the given objective vector.
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

  protected:
    const Projection& _projection; // Projection.
    OptimizationOracleBase* _oracle; // Source oracle.
    soplex::DVectorRational _liftedObjective; // Dense rational lifted objective A^t c
    OptimizationResult _result; // Result structure for result of source oracle.
  };

} /* namespace ipo */

#endif /* IPO_ORACLE_H_ */
