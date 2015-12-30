#ifndef IPO_FACETS_H_
#define IPO_FACETS_H_

#include <set>

#include "oracles.h"
#include "cpu_timer.h"
#include "unique_rational_vectors.h"
#include "spx_gmp.h"

namespace ipo {

  namespace Separation {

    struct Certificate
    {
      VectorSubset pointIndices;
      VectorSubset rayIndices;
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

    class QuietOutput: public OutputBase
    {
    public:
      QuietOutput();
      virtual ~QuietOutput();
    };

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
      CPUTimer _timer;
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

    class Result: public InformationBase
    {
    public:
      Result(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, const VectorSubset& spanningPoints,
          const VectorSubset& spanningRays, const VectorSubset& columnBasis, OptimizationOracleBase* oracle);
      virtual ~Result();

      void inequality(soplex::LPRowRational& inequality) const;
      void inequality(soplex::LPRowSetRational& inequalities) const;
      void certificate(Certificate& certificate) const;
      const soplex::Rational& violation() const;

      bool separatePoint(const Point* targetPoint, OutputBase& output);
      bool separateRay(const Ray* targetRay, OutputBase& output);

      friend class Implementation;
    };

  }

}

#endif /* IPO_FACETS_H_ */
