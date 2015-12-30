#ifndef IPO_AFFINE_HULL_H_
#define IPO_AFFINE_HULL_H_

#include "oracles.h"
#include "spx_gmp.h"
#include "unique_rational_vectors.h"
#include "cpu_timer.h"

namespace ipo {

  namespace AffineHull {

    class Implementation;

    enum Options
    {
      REDUNDANT_EQUATIONS_KEEP = 0, // Given equations are kept.
      REDUNDANT_EQUATIONS_REMOVE = 1, // Only irredundant given equations are kept.
      MASK_REDUNDANT_EQUATIONS = 1,
      CACHE_USE = 0, // Use point / ray cache.
      CACHE_SKIP = 2, // Don't use point / ray cache.
      MASK_CACHE = 2,
      ORACLE_DELAYED = 0, // Use heuristic calls and verify equations at the end.
      ORACLE_ONLY = 4, // Always use oracle calls.
      ORACLE_IMMEDIATE = 8, // Use heuristic calls and verify equation immediately.
      ORACLE_NEVER = 12, // Always use heuristic oracle calls and trust.
      MASK_ORACLE = 12
    };

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
      const VectorSubset& spanningRays() const;
      std::size_t numSpanningRays() const;
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
      bool optionRemoveRedundantEquations() const;
      bool optionUseCache() const;
      int optionOracleUsage() const;

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
      virtual void onBeforeOracleZero(bool forceOptimal);
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize(bool forceOptimal);
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleMinimize(bool forceOptimal);
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

    class QuietOutput: public OutputBase
    {
    public:
      QuietOutput();
      virtual ~QuietOutput();
    };

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
      virtual void onBeforeOracleZero(bool forceOptimal);
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize(bool forceOptimal);
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleMinimize(bool forceOptimal);
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
      virtual void onBeforeOracleZero(bool forceOptimal);
      virtual void onAfterOracleZero(std::size_t numPoints);
      virtual void onBeforeOracleMaximize(bool forceOptimal);
      virtual void onAfterOracleMaximize(std::size_t numPoints, std::size_t numRays);
      virtual void onBeforeOracleMinimize(bool forceOptimal);
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

    class Result: public InformationBase
    {
    public:
      Result();
      virtual ~Result();

      int dimension() const;

      int run(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, soplex::LPRowSetRational& equations,
          OptimizationOracleBase* oracle, OutputBase& output,
          int options = REDUNDANT_EQUATIONS_KEEP | CACHE_USE | ORACLE_DELAYED);

      friend class Implementation;
    };

    int run(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, soplex::LPRowSetRational& equations,
        OptimizationOracleBase* oracle, OutputBase& output,
        int options = REDUNDANT_EQUATIONS_KEEP | CACHE_USE | ORACLE_DELAYED);

  }

}

#endif /* IPO_AFFINE_HULL_H_ */
