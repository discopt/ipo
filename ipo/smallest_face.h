#ifndef IPO_SMALLEST_FACE_H_
#define IPO_SMALLEST_FACE_H_

#include <set>

#include "ipo.h"
#include "spx_gmp.h"
#include "unique_rational_vectors.h"
#include "oracles.h"
#include "affine_hull.h"

namespace ipo {

  /**
   * \brief Computation of the smallest face containing a given point.
   *
   * This namespace contains classes that allow to compute the smallest face
   * that contains a given point.
   * Use it as follows:
   *
   * \code
   * // ..
   * // Create an oracle and a point.
   * // ...
   *
   * UniqueRationalVectors points(oracle->numVariables());
   * UniqueRationalVectors rays(oracle->numVariables());
   *
   * SmallestFace::ProgressOutput smallestFaceOutput;
   * SmallestFace::Result smallestFace(points, rays, oracle);
   *
   * smallestFace.run(point, smallestFaceOutput);
   * std::cout << "Dimension of smallest face: " << smallestFace.dimension() << std::endl;
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

  namespace SmallestFace {

    class Implementation;
    class NormalConeOptimizationOracle;

    class InformationBase
    {
    public:
      InformationBase();
      virtual ~InformationBase();

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
      virtual void onBeforeVerifyElements(std::size_t numVerifications);
      virtual void onAfterVerifyElements(std::size_t numVerifications);
      virtual void onBeforeOracleVerifyElement(std::size_t verifyIndex);
      virtual void onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool verified);
      virtual void onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t numEquations);
      virtual void onConeBeforeSolve(bool stabilizing);
      virtual void onConeAfterSolve(bool stabilizing);
      virtual void onConePenaltyDecrease();
      virtual void onConeBeforeCache();
      virtual void onConeAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onConeBeforeOracleCall(bool forceOptimal);
      virtual void onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays);
      virtual void onConeBeforeAddPoint();
      virtual void onConeAfterAddPoint();
      virtual void onConeBeforeAddRay();
      virtual void onConeAfterAddRay();
      virtual void onEnd();

      virtual AffineHull::OutputBase& normalConeHullOutput() = 0;

      friend class NormalConeOptimizationOracle;
      friend class Implementation;
    };

    /**
     * \brief Quiet output class for smallest-face computation.
     */

    class QuietOutput: public OutputBase
    {
    public:
      QuietOutput();
      virtual ~QuietOutput();

    protected:
      virtual AffineHull::OutputBase& normalConeHullOutput();

      AffineHull::QuietOutput _normalConeHullOutput;
    };

    /**
     * \brief Very verbose output class for smallest-face computation.
     */

    class DebugOutput: public OutputBase
    {
    public:
      DebugOutput();
      virtual ~DebugOutput();

    protected:
      virtual void onStart();
      virtual void onBeforeVerifyElements(std::size_t numVerifications);
      virtual void onAfterVerifyElements(std::size_t numVerifications);
      virtual void onBeforeOracleVerifyElement(std::size_t verifyIndex);
      virtual void onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool verified);
      virtual void onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t numEquations);
      virtual void onConeBeforeSolve(bool stabilizing);
      virtual void onConeAfterSolve(bool stabilizing);
      virtual void onConeBeforeCache();
      virtual void onConeAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onConeBeforeOracleCall(bool forceOptimal);
      virtual void onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays);
      virtual void onEnd();

    protected:
      virtual AffineHull::OutputBase& normalConeHullOutput();

      AffineHull::DebugOutput _normalConeHullOutput;
    };

    /**
     * \brief Pretty output class for smallest-face computation.
     */

    class ProgressOutput: public OutputBase
    {
      class AffineHullOutput: public AffineHull::OutputBase
      {
      public:
        AffineHullOutput(std::size_t indent = 0);
        virtual ~AffineHullOutput();

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
        virtual void onBeforeDirection();
        virtual void onAfterDirection();
        virtual void onPotentialEquation();
        virtual void onEquation();
        virtual void onBeforeVerifyImmediate();
        virtual void onBeforeVerifyDelayed(std::size_t numVerifications);
        virtual void onAfterVerifyDelayed(std::size_t numVerifications);
        virtual void onEnd();

        friend class ProgressOutput;

      protected:
        std::string _indent;
        std::size_t _numVerificationCalls;
        CPUTimer _timer;
        bool _lastApprox;
        double _lastTime;
        double _timeStarted;
        double _timeApproxDirections;
        double _timeExactDirections;
        double _timeCache;
        double _timeFactorization;
        double _timeLP;
        double _timeMainCache;
        double _timeMainHeuristics;
        double _timeMainOracles;
        std::size_t _numMainCache;
        std::size_t _numMainHeuristics;
        std::size_t _numMainOracles;
      };

    public:
      ProgressOutput(std::size_t indent = 0);
      virtual ~ProgressOutput();

    protected:
      virtual void onBeforeVerifyElements(std::size_t numVerifications);
      virtual void onBeforeOracleVerifyElement(std::size_t verifyIndex);
      virtual void onAfterOracleVerifyElement(std::size_t numPoints, std::size_t numRays, bool verified);
      virtual void onAddedInitials(std::size_t numPoints, std::size_t numRays, std::size_t numEquations);
      virtual void onConeBeforeSolve(bool stabilizing);
      virtual void onConeAfterSolve(bool stabilizing);
      virtual void onConePenaltyDecrease();
      virtual void onConeBeforeCache();
      virtual void onConeAfterCache(std::size_t numPoints, std::size_t numRays);
      virtual void onConeBeforeOracleCall(bool forceOptimal);
      virtual void onConeAfterOracleCall(bool forceOptimal, bool feasible, std::size_t numPoints, std::size_t numRays);

    protected:
      virtual AffineHull::OutputBase& normalConeHullOutput();

      AffineHullOutput _normalConeHullOutput;
      std::string _indent;
      std::size_t _numVerifications;
    };

    /**
     * \brief Actual results of a smallest-face computation.
     */

    class Result: public InformationBase
    {
    public:
      Result(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, OptimizationOracleBase* oracle);
      virtual ~Result();

      void getMaximizingObjective(soplex::DSVectorRational& maximizingObjective) const;

      int dimension() const;

      int run(const Point* target, OutputBase& output,
          const std::vector<const soplex::SVectorRational*>& normalConeElements, bool verifyElements);
      int run(const Point* target, OutputBase& output, const soplex::LPRowSetRational& affineHullEquations);
      int run(const Point* target, OutputBase& output);

      bool isVertex(const Point* target, OutputBase& output,
          const std::vector<const soplex::SVectorRational*>& normalConeElements, bool verifyElements);
      bool isVertex(const Point* target, OutputBase& output, const soplex::LPRowSetRational& affineHullEquations);
      bool isVertex(const Point* target, OutputBase& output);

      bool areAdjacent(const Point* first, const Point* second, OutputBase& output,
          const std::vector<const soplex::SVectorRational*>& normalConeElements, bool verifyElements);
      bool areAdjacent(const Point* first, const Point* second, OutputBase& output,
          const soplex::LPRowSetRational& affineHullEquations);
      bool areAdjacent(const Point* first, const Point* second, OutputBase& output);

      friend class Implementation;
    };

  } /* namespace SmallestFace */

} /* namespace ipo */

#endif /* IPO_SMALLEST_FACE_H_ */
