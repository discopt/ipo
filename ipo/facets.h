#ifndef IPO_FACETS_H_
#define IPO_FACETS_H_

#include <set>

#include "common.h"
#include "spx_gmp.h"
#include "timer.h"
#include "oracles.h"

namespace ipo {

  class FacetSeparationState
  {
  public:
    /**
     * \brief Constructor.
     *
     */

    FacetSeparationState();

    /**
     * \brief Destructor.
     */

    virtual ~FacetSeparationState();

    /**
     * \brief Returns the space of the oracle.
     *
     * Returns a const-reference to the space of the oracle.
     */

    virtual const Space& oracleSpace() const = 0;

    /**
     * \brief Returns the space of the polar LP.
     *
     * Returns a const-reference to the space of the polar LP.
     */

    virtual const Space& polarSpace() const = 0;

    /**
     * \brief Returns the feasibility / optimality tolerance of the LP.
     *
     * Returns the feasibility / optimality tolerance of the LP.
     */

    virtual double polarTolerance() const = 0;

    /**
     * \brief Returns the number of spanning points that are currently in the LP.
     *
     * Returns the number of spanning points that are currently in the LP.
     */

    virtual std::size_t polarNumPoints() const = 0;

    /**
     * \brief Returns the number of spanning rays that are currently in the LP.
     *
     * Returns the number of spanning rays found that are currently in the LP.
     */

    virtual std::size_t polarNumRays() const = 0;

    /**
     * \brief Returns the number of points that were added in this iteration.
     *
     * Returns the number of points that were added in this iteration.
     */

    virtual std::size_t numAddedPoints() const = 0;

    /**
     * \brief Returns the number of rays that were added in this iteration.
     *
     * Returns the number of rays that were added in this iteration.
     */

    virtual std::size_t numAddedRays() const = 0;

    /**
     * \brief Returns the number of rows of the current LP.
     *
     * Returns the number of rows of the current LP.
     */

    virtual std::size_t polarNumRowsLP() const = 0;

    /**
     * \brief Returns the number of columns of the current LP.
     *
     * Returns the number of columns of the current LP.
     */

    virtual std::size_t polarNumColumnsLP() const = 0;

    /**
     * \brief Returns the number of nonzeros of the current LP.
     *
     * Returns the number of nonzeros of the current LP.
     */

    virtual std::size_t polarNumNonzerosLP() const = 0;

    /**
     * \brief Returns the last solution of the polar LP as an inequality.
     *
     * Returns the last solution of the polar LP as an inequality.
     */

    virtual LinearConstraint currentInequality() const = 0;

    /**
     * \brief Returns the objective value of the last oracle call.
     *
     * Returns the objective value of the last oracle call, i.e., the smallest right-hand side for currentInequality()
     * that is valid.
     */

    virtual soplex::Rational oracleObjectiveValue() const = 0;

    /**
     * \brief Returns true iff we are currently solving an approximate LP.
     *
     * Returns true iff we are currently solving an approximate LP.
     */

    virtual bool approximateSolve() const = 0;

    /**
     * \brief Returns true iff we are currently solving an exact LP.
     *
     * Returns true iff we are currently solving an exact LP.
     */

    virtual bool exactSolve() const = 0;

    /**
     * \brief Returns the maximum allowed heuristic level of the current oracle call.
     *
     * Returns the maximum allowed heuristic level of the current oracle call.
     */

    virtual HeuristicLevel oracleMaxHeuristicLevel() const = 0;

    /**
     * \brief Returns the minimum allowed heuristic level of the current oracle call.
     *
     * Returns the minimum allowed heuristic level of the current oracle call.
     */

    virtual HeuristicLevel oracleMinHeuristicLevel() const = 0;

    /**
     * \brief Returns the heuristic level of the oracle's last answer.
     *
     * Returns the heuristic level of the oracle's last answer.
     */

    virtual HeuristicLevel oracleResultHeuristicLevel() const = 0;

    /**
     * \brief Returns the point or ray that is currently being added.
     *
     * Returns the point or ray that is currently being added.
     */

    virtual Vector currentVector() const = 0;

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
     * \brief Returns the point or ray to be separated.
     *
     * Returns the point or ray to be separated.
     */

    virtual Vector separationTarget() const = 0;

    /**
     * \brief Returns true iff we are separating a ray.
     *
     * Returns true iff we are separating a ray.
     */

    virtual bool separatingRay() const = 0;

    /**
     * \brief Returns true iff we just separated by a facet.
     *
     * Returns true iff we just separated by a facet (if at all).
     */

    virtual bool separatedByFacet() const = 0;

    /**
     * \brief Returns true iff we just separated by an equation.
     *
     * Returns true iff we just separated by an equation (if at all).
     */

    virtual bool separatedByEquation() const = 0;
  };

  /**
   * \brief Base class for an observer for facet computations.
   *
   * Base class for an observer for facet computations.
   */

  class FacetSeparationHandler
  {
  public:
    enum Event
    {
      LP_BEGIN, // Before (re)solving the current polar LP.
      LP_END, // After (re)solving the current polar LP.
      ORACLE_BEGIN, // Before an oracle call.
      ORACLE_END = ORACLE_BEGIN + 1, // After an oracle call.
      POINTS_BEGIN, // Before adding point-constraints to the LP.
      POINT, // Every point-constraint added.
      POINTS_END, // After adding point-constraints to the LP.
      RAYS_BEGIN, // Before adding ray-constraints to the LP.
      RAY, // Every ray-constraint added.
      RAYS_END, // After adding ray-constraints to the LP.
      BEGIN, // Algorithm started.
      INITIALIZED, // LP initialized.
      APPROXIMATE_SOLVE_BEGIN, // Before solving approximate polar LP to optimality.
      APPROXIMATE_SOLVE_END, // After solving approximate polar LP to optimality.
      EXACT_SOLVE_BEGIN, // Before solving exact polar LP to optimality.
      EXACT_SOLVE_END, // After solving exact polar LP to optimality.
      END // Algorithm finished.
    };

    /**
     * \brief Default constructor.
     *
     * Default constructor.
     */

    FacetSeparationHandler();

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~FacetSeparationHandler();

    /**
     * \brief This method is called by the algorithm.
     *
     * This method is called by the algorithm in certain steps.
     */

    virtual void notify(Event event, FacetSeparationState& state) = 0;
  };


  /**
   * \brief Class for debugging facet computations.
   *
   * Class for debugging facet computations. It prints at least one line for each notification.
   */

  class DebugFacetSeparationHandler: public FacetSeparationHandler
  {
  public:
    /**
     * \brief Constructor.
     *
     * Constructor.
     *
     * \param stream Where the output will be written to.
     */

    DebugFacetSeparationHandler(std::ostream& stream, bool printPointsAndRays = false, bool printInequalities = false);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~DebugFacetSeparationHandler();

    /**
     * \brief Notification method, see \ref FacetSeparationHandler.
     *
     * Notification method, see \ref FacetSeparationHandler.
     */

    virtual void notify(Event event, FacetSeparationState& state);

  protected:
    std::ostream& _stream;
    bool _printPointsAndRays;
    bool _printInequalities;
  };

  /**
    * \brief Class for collecting statistics of a facet computation.
    *
    * Class for collecting statistics of a facet computation.
    */

  class StatisticsFacetSeparationHandler: public FacetSeparationHandler
  {
  public:
    /**
      * \brief Constructor.
      *
      * Constructor.
      */

    StatisticsFacetSeparationHandler();

    /**
      * \brief Destructor.
      *
      * Destructor.
      */

    virtual ~StatisticsFacetSeparationHandler();

    /**
      * \brief Notification method, see \ref FacetSeparationHandler.
      *
      * Notification method, see \ref FacetSeparationHandler.
      */

    virtual void notify(Event event, FacetSeparationState& state);

    /**
     * \brief Resets all statistics.
     *
     * Resets all statistics.
     */

    void reset();

    /**
     * \brief Returns the number of oracle queries.
     *
     * Returns the number of oracle queries.
     */

    inline std::size_t numOracleQueries() const
    {
      return _numOracleQueries;
    }

    /**
     * \brief Returns the vector storing how often an oracle returned with each heuristicLevel.
     *
     * Returns the vector storing how often an oracle returned with each heuristicLevel.
     */

    inline const std::vector<std::size_t>& numHeuristicLevelAnswers() const
    {
      return _numHeuristicLevelAnswers;
    }

    /**
     * \brief Returns how many approximate LPs were solved.
     *
     * Returns how many approximate LPs were solved.
     */

    inline std::size_t numApproximateLPs() const
    {
      return _numApproximateLPs;
    }

    /**
     * \brief Returns how many exact LPs were solved.
     *
     * Returns how many exact LPs were solved.
     */

    inline std::size_t numExactLPs() const
    {
      return _numExactLPs;
    }

    /**
     * \brief Returns the total running time for solving approximate LPs.
     *
     * Returns the total running time for solving approximate LPs.
     */

    inline double timeApproximateLPs() const
    {
      return _timeApproximateLPs;
    }

    /**
     * \brief Returns the total running time for solving exact LPs.
     *
     * Returns the total running time for solving exact LPs.
     */

    inline double timeExactLPs() const
    {
      return _timeExactLPs;
    }

    /**
     * \brief Returns the total running time for oracles.
     *
     * Returns the total running time for oracles.
     */

    inline double timeOracles() const
    {
      return _timeOracles;
    }

    /**
     * \brief Returns the total running time of the algorithm.
     *
     * Returns the total running time of the algorithm.
     */

    inline double timeAll() const
    {
      return _timeAll;
    }

  protected:
    Timer _timer;
    std::size_t _numOracleQueries;
    std::vector<std::size_t> _numHeuristicLevelAnswers;
    std::size_t _numApproximateLPs;
    std::size_t _numExactLPs;
    double _timeApproximateLPs;
    double _timeExactLPs;
    double _timeOracles;
    double _timeAll;
    double _timeLastBegin;
    double _timeLastEvent;
  };




  /**
   * \brief Separates a given \p point by a facet or an equation.
   *
   * Separates a given \p point by a facet or an equation.
   *
   * \param oracle      Oracle defining the polyhedron P.
   * \param point       Point to be separated.
   * \param spanning    Points and rays spanning P's affine hull.
   * \param handlers    Set of \ref FacetSeparationHandler objects whose notify() methods are called appropriately.
   * \param constraint  If the \p point lies outside P, this variable will contain the separating inequality.
   * \param certificate If not \c NULL, it will contain points and rays that span the facet or equation.
   *
   * \return True iff \p point can be separated.
   */

  bool separatePoint(const std::shared_ptr<OracleBase>& oracle, const Vector& point, const InnerDescription& spanning,
    std::vector<FacetSeparationHandler*>& handlers, LinearConstraint& constraint, InnerDescription* certificate = NULL);

  /**
   * \brief Separates a given \p ray by a facet or an equation.
   *
   * Separates a given \p ray by a facet or an equation.
   *
   * \param oracle      Oracle defining the polyhedron P.
   * \param ray         Ray to be separated.
   * \param spanning    Points and rays spanning P's affine hull.
   * \param handlers    Set of \ref FacetSeparationHandler objects whose notify() methods are called appropriately.
   * \param constraint  If the \p ray lies outside P's recession cone, this variable will contain the separating inequality.
   * \param certificate If not \c NULL, it will contain points and rays that span the facet or equation.
   *
   * \return True iff \p ray can be separated.
   */

  bool separateRay(const std::shared_ptr<OracleBase>& oracle, const Vector& ray, const InnerDescription& spanning,
    std::vector<FacetSeparationHandler*>& handlers, LinearConstraint& constraint, InnerDescription* certificate = NULL);

  ///////////////// OLD ////////////////////

  /**
   * \brief Facet-separation for a polyhedron.
   *
   * This namespace contains classes that allow the separation of facets of a polyhedron.
   * If the given point is not in the affine hull, it may also yield a violated equation.
   * Use it as follows:
   *
   * \code
   * // ..
   * // Create an oracle and a point to be separated.
   * // ...
   *
   * // We first compute the affine hull.
   *
   * UniqueRationalVectors points(oracle->numVariables());
   * UniqueRationalVectors rays(oracle->numVariables());
   * soplex::LPRowSetRational equations;
   * AffineHull::QuietOutput hullOutput;
   * AffineHull::Result hull;
   * hull.run(points, rays, equations, oracle, hullOutput);
   *
   * // We now create the separator...
   *
   * Separation::ProgressOutput separateOutput;
   * Separation::Result separator(points, rays, hull.spanningPoints(),
   *   hull.spanningRays(), hull.basicColumns(), oracle);
   *
   * // ... and separate a point.
   *
   * separator.separatePoint(point, separateOutput);
   * if (separate.violation() <= 0)
   *   std::cout << "The ";
   *
   * // We extract the inequality and print it.
   *
   * soplex::LPRowRational inequality;
   * separator.inequality(inequality);
   * oracle->printRow(std::cout, inequality);
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

  namespace Separation {

    /**
     * \brief A certificate structure for faces.
     *
     * Structure containing points and rays that span a face.
     */

    struct Certificate
    {
      std::vector<Vector> points;
      std::vector<Vector> rays;
    };

    class Implementation;

    class InformationBase
    {
    public:
      InformationBase();
      virtual ~InformationBase();

      bool separatingPoint() const;

      inline bool separatingRay() const
      {
        return !separatingPoint();
      }

      bool separatedFacet() const;
      bool separatedEquation() const;

      const std::string& oracleName() const;
      std::size_t numVariables() const;

      std::size_t numRowsLP() const;
      std::size_t numColumnsLP() const;
      std::size_t numNonzerosLP() const;

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
      virtual void onEnd(bool separated);
      virtual void onBeforeSolve(bool stabilizing, double penalty);
      virtual void onAfterSolve(bool stabilizing, double penalty, double mainObjective, double penaltyCosts);
      virtual void onPenaltyDecrease(double penalty);
      virtual void onBeforeCache();
      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleCall(bool forceOptimal);
      virtual void onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
          bool lastIteration);
      virtual void onBeforeAddPoint();
      virtual void onAfterAddPoint();
      virtual void onBeforeAddRay();
      virtual void onAfterAddRay();

      friend class Implementation;
    };

    /**
     * \brief Quiet output class for the facet separation.
     */

    class QuietOutput: public OutputBase
    {
    public:
      QuietOutput();
      virtual ~QuietOutput();
    };

    /**
     * \brief Very verbose output class for the facet separation.
     */

    class DebugOutput: public OutputBase
    {
    public:
      DebugOutput();
      virtual ~DebugOutput();

    protected:
      virtual void onStart();
      virtual void onEnd(bool separated);
      virtual void onBeforeSolve(bool stabilizing, double penalty);
      virtual void onAfterSolve(bool stabilizing, double penalty, double mainObjective, double penaltyCosts);
      virtual void onPenaltyDecrease(double penalty);
      virtual void onBeforeCache();
      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleCall(bool forceOptimal);
      virtual void onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
          bool lastIteration);
      virtual void onBeforeAddPoint();
      virtual void onAfterAddPoint();
      virtual void onBeforeAddRay();
      virtual void onAfterAddRay();
    };

    /**
     * \brief Pretty output class for the facet separation.
     */

    class ProgressOutput: public OutputBase
    {
    public:
      ProgressOutput(std::size_t indent = 0);
      virtual ~ProgressOutput();

      void printStatistics();

    protected:
      virtual void onStart();
      virtual void onEnd(bool separated);
      virtual void onBeforeSolve(bool stabilizing, double penalty);
      virtual void onAfterSolve(bool stabilizing, double penalty, double mainObjective, double penaltyCosts);
      virtual void onPenaltyDecrease(double penalty);
      virtual void onBeforeCache();
      virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleCall(bool forceOptimal);
      virtual void onAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays,
          bool lastIteration);

    protected:
      double timeStamp();

      std::string _indent;
      Timer _timer;
      std::size_t _numStabilizationLP;
      std::size_t _numMainLP;
      std::size_t _numCache;
      std::size_t _numHeuristics;
      std::size_t _numOracles;
      double _lastTime;
      double _timeOverall;
      double _timeStarted;
      double _timeStabilizationLP;
      double _timeMainLP;
      double _timeCache;
      double _timeHeuristics;
      double _timeOracles;
    };

    /**
     * \brief Actual results of a facet separation.
     */

    class Result: public InformationBase
    {
    public:
      Result(const std::vector<Vector>& spanningPoints, const std::vector<Vector>& spanningRays,
        const std::vector<std::size_t>& columnBasis, const std::shared_ptr<OracleBase>& oracle);
      virtual ~Result();

      LinearConstraint inequality() const;
      void certificate(Certificate& certificate) const;
      const Rational& violation() const;

      bool separatePoint(const Vector& targetPoint, OutputBase& output);
      bool separateRay(const Vector& targetRay, OutputBase& output);

      friend class Implementation;
    };

  } /* namespace Separation */

} /* namespace ipo */

#endif /* IPO_FACETS_H_ */
