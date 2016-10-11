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
     * \brief Returns true iff we are separating a ray.
     * 
     * Returns true iff we are separating a ray.
     */

    virtual bool separatingRay() const = 0;   
  };

  /**
   * \brief Base class for an observer for affine-hull computations.
   * 
   * Base class for an observer for affine-hull computations.
   */

  class FacetSeparationHandler
  {
  public:
    enum Event
    {
      BEGIN,
      EQUATIONS_INITIALIZED,
      LOOP,
      APPROXIMATE_SOLVE_BEGIN,
      APPROXIMATE_SOLVE_END,
      EXACT_SOLVE_BEGIN,
      EXACT_SOLVE_END,
      ORACLE_BEGIN,
      ORACLE_END = ORACLE_BEGIN + 1,
      POINT_BEGIN,
      POINT_END,
      RAY_BEGIN,
      RAY_END,
      END,
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
