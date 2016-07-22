#ifndef IPO_FACETS_H_
#define IPO_FACETS_H_

#include <set>

#include "common.h"
#include "spx_gmp.h"
#include "cpu_timer.h"
#include "unique_rational_vectors.h"
#include "oracles.h"

namespace ipo {

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
     * Structure containing indices of points and directions that span a face.
     */

    struct Certificate
    {
      VectorSubset pointIndices;
      VectorSubset directionIndices;
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

    /**
     * \brief Actual results of a facet separation.
     */

    class Result: public InformationBase
    {
    public:
      Result(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& rays, const VectorSubset& spanningPoints,
          const VectorSubset& spanningRays, const VectorSubset& columnBasis, OracleBase* oracle);
      virtual ~Result();

      void inequality(soplex::LPRowRational& inequality) const;
      void inequality(soplex::LPRowSetRational& inequalities) const;
      void certificate(Certificate& certificate) const;
      const soplex::Rational& violation() const;

      bool separatePoint(const Point* targetPoint, OutputBase& output);
      bool separateRay(const Direction* targetRay, OutputBase& output);

      friend class Implementation;
    };

  } /* namespace Separation */

} /* namespace ipo */

#endif /* IPO_FACETS_H_ */
