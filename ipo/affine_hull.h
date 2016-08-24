#ifndef IPO_AFFINE_HULL_H_
#define IPO_AFFINE_HULL_H_

#include <map>

#include "common.h"
#include "cpu_timer.h"
#include "oracles.h"

namespace ipo {

  /**
   * \brief Interface for information-retrieval during an affine-hull computation.
   * 
   * Interface for information-retrieval during an affine-hull computation.
   */
  
  class AffineHullState
  {
  public:
    /**
     * \brief Constructor.
     * 
     */

    AffineHullState();

    /**
     * \brief Destructor.
     */

    virtual ~AffineHullState();

    /**
     * \brief Returns the space of the oracle.
     * 
     * Returns a const-reference to the space of the oracle.
     */

    virtual const Space& space() const = 0;

    /**
     * \brief Returns the number of spanning points found so far.
     * 
     * Returns the number of spanning points found so far.
     */

    virtual std::size_t numPoints() const = 0;

    /**
     * \brief Returns the number of spanning rays found so far.
     * 
     * Returns the number of spanning rays found so far.
     */

    virtual std::size_t numRays() const = 0;

    /**
     * \brief Returns the number of valid independent equations found so far.
     * 
     * Returns the number of valid independent equations found so far.
     */

    virtual std::size_t numEquations() const = 0;

    /**
     * \brief Returns the number of candidate equations.
     * 
     * Returns the number of candidate equations. Together with the equations (see numEquations()) they are linearly independent.
     */

    virtual std::size_t numCandidateEquations() const = 0;

    /**
     * \brief Returns whether an oracle returned infeasible.
     * 
     * Returns whether an oracle returned infeasible.
     */
    
    virtual bool infeasible() const = 0;

    /**
     * \brief Returns the best lower bound on the dimension known so far.
     * 
     * Returns the best lower bound on the dimension known so far.
     */
    
    inline long lowerBound() const
    {
      return long(numPoints() + numRays()) - 1;
    }

    /**
     * \brief Returns the best exact upper bound on the dimension known so far.
     * 
     * Returns the best exact upper bound on the dimension known so far.
     */

    inline long upperBound() const
    {
      if (infeasible())
        return -1;
      else
        return long(space().dimension() - numEquations());
    }

    /**
     * \brief Returns the best approximate upper bound on the dimension known so far.
     * 
     * Returns the best approximate upper bound on the dimension known so far. Is correct if the candidate equations are indeed
     * valid.
     */

    inline long candidateUpperBound() const
    {
      if (infeasible())
        return -1;
      else
        return long(space().dimension() - numEquations() - numCandidateEquations());
    }

    /**
     * \brief Returns the index of the candidate equation that is being verified.
     * 
     * Returns the index of the candidate equation that is being verified.
     */

    virtual std::size_t verifyIndex() const = 0;

    /**
     * \brief Returns true iff the verification was successful so far.
     * 
     * Returns true iff the verification was successful so far.
     */

    virtual bool verifySuccess() const = 0;

    /**
     * \brief Returns the number of approximate LU solves in the last iteratior.
     * 
     * Returns the number of approximate LU solves in the last iteratior.
     */

    virtual std::size_t directionApproximateSolves() const = 0;
    
    /**
     * \brief Returns the number of exact LU solves in the last iteratior.
     * 
     * Returns the number of exact LU solves in the last iteratior.
     */

    virtual std::size_t directionExactSolves() const = 0;
    
    /**
     * \brief Returns the number of nonzeros of the direction vector.
     * 
     * Returns the number of nonzeros of the direction vector.
     */

    virtual std::size_t directionNonzeros() = 0;
    
    /**
     * \brief Returns the number of bits of the direction vector.
     * 
     * Returns the number of bits of the direction vector.
     */

    virtual std::size_t directionBitsize() = 0;
    
    /**
     * \brief Returns the direction vector.
     * 
     * Returns the direction vector.
     */

    virtual const soplex::VectorRational& directionVector() = 0;
    

    /**
     * \brief Returns the maximum allowed heuristic level of the current oracle call.
     * 
     * Returns the maximum allowed heuristic level of the current oracle call.
     */

    virtual std::size_t oracleMaxHeuristicLevel() const = 0;

    /**
     * \brief Returns the minimum allowed heuristic level of the current oracle call.
     * 
     * Returns the minimum allowed heuristic level of the current oracle call.
     */

    virtual std::size_t oracleMinHeuristicLevel() const = 0;

    /**
     * \brief Returns the heuristic level of the oracle's last answer.
     * 
     * Returns the heuristic level of the oracle's last answer.
     */

    virtual std::size_t oracleResultHeuristicLevel() const = 0;

    /**
     * \brief Returns the number of points returned by the last oracle call.
     * 
     * Returns the number of points returned by the last oracle call.
     */

    virtual std::size_t oracleNumPoints() const = 0;

    /**
     * \brief Returns the number of rays returned by the last oracle call.
     * 
     * Returns the number of rays returned by the last oracle call.
     */

    virtual std::size_t oracleNumRays() const = 0;

    /**
     * \brief Returns the current set of spanning points.
     * 
     * Returns a const-reference to the current set of spanning points.
     */

    virtual const std::vector<Vector>& spanningPoints() const = 0;

    /**
     * \brief Returns the current set of spanning rays.
     * 
     * Returns a const-reference to the current set of spanning rays.
     */

    virtual const std::vector<Vector>& spanningRays() const = 0;

    /**
     * \brief Returns the current set of valid independent equations.
     * 
     * Returns a const-reference to the current set of valid independent equations.
     */

    virtual const std::vector<LinearConstraint>& equations() const = 0;

    /**
     * \brief Returns the current set of candidate equations.
     * 
     * Returns a const-reference to the current set of candidate equations.
     */

    virtual const std::vector<LinearConstraint>& candidateEquations() const = 0;

  };
  
  /**
   * \brief Base class for an observer for affine-hull computations.
   * 
   * Base class for an observer for affine-hull computations.
   */

  class AffineHullHandler
  {
  public:
    enum Event
    {
      BEGIN,
      EQUATIONS_INITIALIZED,
      LOOP,
      DIRECTIONS_APPROXIMATE_BEGIN,
      DIRECTIONS_APPROXIMATE_END,
      DIRECTIONS_EXACT_BEGIN,
      DIRECTIONS_EXACT_END,
      ORACLE_ZERO_BEGIN,
      ORACLE_ZERO_END = ORACLE_ZERO_BEGIN + 1,
      ORACLE_MAXIMIZE_BEGIN,
      ORACLE_MAXIMIZE_END = ORACLE_MAXIMIZE_BEGIN + 1,
      ORACLE_MINIMIZE_BEGIN,
      ORACLE_MINIMIZE_END = ORACLE_MINIMIZE_BEGIN + 1,
      ORACLE_VERIFY_BEGIN,
      ORACLE_VERIFY_END = ORACLE_VERIFY_BEGIN + 1,
      POINT_BEGIN,
      POINT_END,
      RAY_BEGIN,
      RAY_END,
      EQUATION_CANDIDATE,
      EQUATION_FINAL,
      VERIFY_BEGIN,
      VERIFY_END,
      END,
    };

    /**
     * \brief Default constructor.
     * 
     * Default constructor.
     */

    AffineHullHandler();

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    virtual ~AffineHullHandler();

    /**
     * \brief This method is called by the algorithm.
     * 
     * This method is called by the algorithm in certain steps.
     */

    virtual void notify(Event event, AffineHullState& state) = 0;
  };
  
  /**
   * \brief Computes the affine hull of a polyhedron defined by the given \p oracle.
   * 
   * Computes the affine hull of a polyhedron defined by the given \p oracle.
   * 
   * \param handlers              Set of \ref AffineHullHandler objects whose notify() methods are called appropriately.
   * \param oracle                Oracle to be used.
   * \param givenEquations        Set of equations whose validity assumed without checking.
   * \param resultPoints          Set of affinely independent points that, together with \p resultRays, spans the affine hull.
   * \param resultRays            Set of linearly independent rays that, together with \p resultPoints, spans the affine hull.
   * \param resultEquations       Set of linearly independent equations describing the affine hull.
   * \param lastCheapHeuristic    Smallest heuristic level for which max-min alternating is performed.
   * \param lastModerateHeuristic Smallest heuristic level to be used for collecting equation candidates.
   * \param approximateDirections If true, computes directions approximately and estimates their sparsity to decide which
   *                              direction vectors should be computed exactly.
   */

  void affineHull(std::vector<AffineHullHandler*>& handlers, const std::shared_ptr<OracleBase>& oracle,
    const std::vector<LinearConstraint>& givenEquations, std::vector<Vector>& resultPoints, std::vector<Vector>& resultRays, 
    std::vector<LinearConstraint>& resultEquations, std::size_t lastCheapHeuristic = std::numeric_limits<std::size_t>::max(), 
    std::size_t lastModerateHeuristic = 0, bool approximateDirections = true);

  class DebugAffineHullHandler: public AffineHullHandler
  {
  public:
    /**
     * \brief Constructor.
     * 
     * Constructor.
     * 
     * \param stream Where the output will be written to.
     */

    DebugAffineHullHandler(std::ostream& stream, bool printPointsAndRays = false, bool printEquations = false,
      bool printDirections = false);

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    virtual ~DebugAffineHullHandler();

    /**
     * \brief Notification method, see \ref AffineHullHandler.
     * 
     * Notification method, see \ref AffineHullHandler.
     */

    virtual void notify(Event event, AffineHullState& state);

  protected:
    std::ostream& _stream;
    bool _printPointsAndRays;
    bool _printEquations;
    bool _printDirections;
  };

  /**
   * \brief Computation of the affine hull of a polyhedron.
   *
   * This namespace contains classes that allow to compute the affine hull of a polyhedron.
   * Use it as follows:
   *
   * \code
   * // ..
   * // Create an oracle.
   * // ...
   *
   * std::vector<LinearConstraint> equations;
   * AffineHull::ProgressOutput hullOutput;
   * AffineHull::Result hull;
   * hull.run(equations, oracle, hullOutput);
   *
   * std::cout << "Dimension: " << hull.dimension() << std::endl;
   *
   * \endcode
   *
   * There are different output classes, namely
   * \li
   *   \ref QuietOutput
   * \li
   *   \ref ProgressOutput
   * \li
   *   \ref DebugOutput
   */

  namespace AffineHull {

    class Implementation;

    class InformationBase
    {
    public:
      InformationBase();
      virtual ~InformationBase();

      const std::string& oracleName() const;
      std::size_t numVariables() const;
      int dimensionLowerBound() const;
      int dimensionUnsafeUpperBound() const;
      int dimensionSafeUpperBound() const;
      const Vector& spanningPoint(std::size_t i) const;
      std::size_t numSpanningPoints() const;
      const Vector& spanningRay(std::size_t i) const;
      std::size_t numSpanningRays() const;
      const std::vector<std::size_t>& basicColumns() const;
      const std::vector<LinearConstraint>& irredundantEquations() const;
      const std::vector<LinearConstraint>& redundantEquations() const;
      const std::vector<LinearConstraint>& potentialEquations() const;
      std::size_t numIrredundantEquations() const;
      std::size_t numRedundantEquations() const;
      std::size_t numPotentialEquations() const;
      std::size_t numCacheQueries() const;
      std::size_t numCacheHits() const;
      std::size_t numHeuristicCalls() const;
      std::size_t numOracleCalls() const;
      std::size_t numApproximateDirectionSolves() const;
      std::size_t numExactDirectionSolves() const;
      std::size_t lastDirectionBitsize() const;
      std::size_t maxDirectionBitsize() const;

      friend class Implementation;
      friend class Result;

    protected:
      bool hasImplementation() const;

    private:
      void ensureImplementation() const;

      Implementation* _implementation;
    };

    class OutputBase: protected InformationBase
    {
    public:
      OutputBase();
      virtual ~OutputBase();

    protected:
      bool isRunning() const;

      virtual void onStart();
      virtual void onAddedInitialEquations(std::size_t numAllEquations);
      virtual void onBeforeApproximateDirections();
      virtual void onAfterApproximateDirections(std::size_t numComputed);
      virtual void onBeforeExactDirections();
      virtual void onAfterExactDirections(std::size_t numComputed);
      virtual void onBeforeCache();
      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleZero();
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize();
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleMinimize();
      virtual void onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleVerify(std::size_t verifyIndex);
      virtual void onAfterOracleVerify(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforePoint(bool twoPoints);
      virtual void onAfterPoint(bool twoPoints);
      virtual void onBeforeRay();
      virtual void onAfterRay();
      virtual void onPotentialEquation();
      virtual void onEquation();
      virtual void onBeforeVerifyImmediate();
      virtual void onAfterVerifyImmediate(bool success);
      virtual void onBeforeVerifyDelayed(std::size_t numVerifications);
      virtual void onAfterVerifyDelayed(std::size_t numVerifications);
      virtual void onRemovedRedundantEquations(std::size_t numRemoved);
      virtual void onEnd();

      friend class Implementation;
    };

    /**
     * \brief Quiet output class for affine hull computation.
     */

    class QuietOutput: public OutputBase
    {
    public:
      QuietOutput();
      virtual ~QuietOutput();
    };

    /**
     * \brief Pretty output class for affine hull computation.
     */

    class ProgressOutput: public OutputBase
    {
    public:
      ProgressOutput(std::size_t indent = 0);
      virtual ~ProgressOutput();

    protected:
      double timeStamp();
      void onProgress();

      virtual void onStart();
      virtual void onAddedInitialEquations(std::size_t numAllEquations);
      virtual void onBeforeApproximateDirections();
      virtual void onAfterApproximateDirections(std::size_t numComputed);
      virtual void onBeforeExactDirections();
      virtual void onAfterExactDirections(std::size_t numComputed);
      virtual void onBeforeCache();
      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleZero();
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize();
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleMinimize();
      virtual void onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleVerify(std::size_t verifyIndex);
      virtual void onAfterOracleVerify(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforePoint(bool twoPoints);
      virtual void onAfterPoint(bool twoPoints);
      virtual void onBeforeRay();
      virtual void onAfterRay();
      virtual void onPotentialEquation();
      virtual void onEquation();
      virtual void onBeforeVerifyImmediate();
      virtual void onBeforeVerifyDelayed(std::size_t numVerifications);
      virtual void onAfterVerifyDelayed(std::size_t numVerifications);
      virtual void onEnd();

    protected:
      std::string _indent;
      std::size_t _numVerificationCalls;
      CPUTimer _timer;
      bool _lastHeuristic;
      double _lastTime;
      double _timeStarted;
      double _timeApproxDirections;
      double _timeExactDirections;
      double _timeCache;
      double _timeHeuristics;
      double _timeOracles;
      double _timeFactorization;
    };

    /**
     * \brief Very verbose output class for affine hull computation.
     */

    class DebugOutput: public OutputBase
    {
    public:
      DebugOutput();
      virtual ~DebugOutput();

    protected:
      void printStatus();
      void timeStamp(const std::string& event = "", bool printNewline = true);

      virtual void onStart();
      virtual void onAddedInitialEquations(std::size_t numAllEquations);
      virtual void onBeforeApproximateDirections();
      virtual void onAfterApproximateDirections(std::size_t numComputed);
      virtual void onBeforeExactDirections();
      virtual void onAfterExactDirections(std::size_t numComputed);
      virtual void onBeforeCache();
      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleZero();
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize();
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleMinimize();
      virtual void onAfterOracleMinimize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleVerify(std::size_t verifyIndex);
      virtual void onAfterOracleVerify(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforePoint(bool twoPoints);
      virtual void onAfterPoint(bool twoPoints);
      virtual void onBeforeRay();
      virtual void onAfterRay();
      virtual void onPotentialEquation();
      virtual void onEquation();
      virtual void onBeforeVerifyImmediate();
      virtual void onAfterVerifyImmediate(bool success);
      virtual void onBeforeVerifyDelayed(std::size_t numVerifications);
      virtual void onAfterVerifyDelayed(std::size_t numVerifications);
      virtual void onRemovedRedundantEquations(std::size_t numRemoved);
      virtual void onEnd();

    protected:
      double _lastTime;
      std::map<std::string, double> _times;
      CPUTimer _timer;
    };

    /**
     * \brief Actual results of an affine hull computation.
     */

    class Result: public InformationBase
    {
    public:
      Result();
      virtual ~Result();

      int dimension() const;

      int run(std::vector<LinearConstraint>& equations, const std::shared_ptr<OracleBase>& oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification = 0, bool removeRedundantEquations = false);

      friend class Implementation;
    };

    int run(std::vector<LinearConstraint>& equations, const std::shared_ptr<OracleBase>& oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification = 0, bool removeRedundantEquations = false);

  } /* namespace AffineHull */

} /* namespace ipo */

#endif /* IPO_AFFINE_HULL_H_ */
