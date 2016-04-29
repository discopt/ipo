#ifndef IPO_AFFINE_HULL_H_
#define IPO_AFFINE_HULL_H_

#include "ipo.h"
#include "spx_gmp.h"
#include "cpu_timer.h"
#include "unique_rational_vectors.h"
#include "oracles.h"

namespace ipo {

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
   * UniqueRationalVectors points(oracle->numVariables());
   * UniqueRationalVectors directions(oracle->numVariables());
   * soplex::LPRowSetRational equations;
   * AffineHull::ProgressOutput hullOutput;
   * AffineHull::Result hull;
   * hull.run(points, directions, equations, oracle, hullOutput);
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
      const VectorSubset& spanningPoints() const;
      std::size_t numSpanningPoints() const;
      const VectorSubset& spanningDirections() const;
      std::size_t numSpanningDirections() const;
      const VectorSubset& basicColumns() const;
      const VectorSubset& irredundantEquations() const;
      std::size_t numIrredundantEquations() const;
      const VectorSubset& potentialEquations() const;
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
      virtual void onAfterCache(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleZero();
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize();
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleMinimize();
      virtual void onAfterOracleMinimize(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleVerify(std::size_t verifyIndex);
      virtual void onAfterOracleVerify(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforePoint(bool twoPoints);
      virtual void onAfterPoint(bool twoPoints);
      virtual void onBeforeDirection();
      virtual void onAfterDirection();
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
      virtual void onAfterCache(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleZero();
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize();
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleMinimize();
      virtual void onAfterOracleMinimize(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleVerify(std::size_t verifyIndex);
      virtual void onAfterOracleVerify(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforePoint(bool twoPoints);
      virtual void onAfterPoint(bool twoPoints);
      virtual void onBeforeDirection();
      virtual void onAfterDirection();
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
      virtual void onAfterCache(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleZero();
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize();
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleMinimize();
      virtual void onAfterOracleMinimize(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforeOracleVerify(std::size_t verifyIndex);
      virtual void onAfterOracleVerify(std::size_t numPoints, std::size_t numDirections);
      virtual void onBeforePoint(bool twoPoints);
      virtual void onAfterPoint(bool twoPoints);
      virtual void onBeforeDirection();
      virtual void onAfterDirection();
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

      int run(UniqueRationalVectorsBase& points,
        UniqueRationalVectorsBase& directions,soplex::LPRowSetRational& equations,
        OracleBase* oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification = 0, bool removeRedundantEquations = false);

      friend class Implementation;
    };

    int run(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions, soplex::LPRowSetRational& equations,
        OracleBase* oracle, OutputBase& output,
        std::size_t minHeuristicBeforeVerification = 0, bool removeRedundantEquations = false);

  } /* namespace AffineHull */

} /* namespace ipo */

#endif /* IPO_AFFINE_HULL_H_ */
