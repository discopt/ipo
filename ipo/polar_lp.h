#ifndef IPO_POLAR_LP_H_
#define IPO_POLAR_LP_H_

#include <vector>
#include <set>
#include <map>

#include "common.h"
#include "unique_rational_vectors.h"
#include "oracles.h"

#define IPO_POLAR_STABIL_HIST_SIZE 2

namespace ipo {

  /// Class for solving polar LPs using lazy separation, cut aging and stabilization.

  class PolarLP
  {
  public:
    struct Basis
    {
      std::vector<soplex::SPxSolver::VarStatus> columnStatus;
      std::vector<soplex::SPxSolver::VarStatus> constraintStatus;
      std::set<std::size_t> tightPoints;
      std::set<std::size_t> tightRays;
    };

    PolarLP(UniqueRationalVectorsBase& points, UniqueRationalVectorsBase& directions,
      OracleBase* oracle, double initialPenalty = 1024.0, int maxAge = 30);
    virtual ~PolarLP();

  protected:
    inline std::size_t n() const
    {
      return _n;
    }

    inline std::size_t numConstraints() const
    {
      return _constraintsToRows.size();
    }

    inline std::size_t numRows() const
    {
      return _stabilizing ? _stabLP->numRowsReal() : _mainLP->numRowsRational();
    }

    inline std::size_t numColumns() const
    {
      return _stabilizing ? _stabLP->numColsReal() : _mainLP->numColsRational();
    }

    inline std::size_t numNonzeros() const
    {
      return _stabilizing ? _stabLP->numNonzerosReal() : _mainLP->numRowsRational();
    }

    inline double lastMainObjective() const
    {
      return _lastMainObjective;
    }

    inline double lastStabilizationPenaltyCosts() const
    {
      return _lastPenaltyCosts;
    }

    inline double currentStabilizationPenalty() const
    {
      return _stabPenalty;
    }

    void setBounds(std::size_t column, const soplex::Rational& lower, const soplex::Rational& upper);
    std::size_t addConstraint(const soplex::Rational& lhs, const soplex::SVectorRational& row,
        const soplex::Rational& rhs);
    void updateConstraint(std::size_t index, const soplex::LPRowRational& row);
    void updateConstraint(std::size_t index, const soplex::Rational& lhs, const soplex::SVectorRational& normal,
        const soplex::Rational& rhs);
    void updateConstraint(std::size_t index, const soplex::Rational& lhs, const soplex::VectorRational& normal,
        const soplex::Rational& rhs);
    void updateObjective(const soplex::VectorRational& objective);
    std::size_t addPointContraint(std::size_t index);
    std::size_t addRayContraint(std::size_t index);

    void setBasis(const Basis& basis);
    void getBasis(Basis& basis);
    soplex::Rational getObjectiveValue();
    void getPrimalSolution(soplex::VectorRational& solution);
    void getPrimalSolution(soplex::DSVectorRational& solution);
    void stabilizedPresolve();
    void optimize(bool perturbeObjective);
    void reoptimizePerturbed();

    /// Callbacks for the output.

    virtual void onBeforeSolve(bool stabilizing);
    virtual void onAfterSolve(bool stabilizing);
    virtual void onPenaltyDecrease();
    virtual void onBeforeCache();
    virtual void onAfterCache(std::size_t numPoints, std::size_t numRays);
    virtual void onBeforeOracleCall();
    virtual void onAfterOracleCall(bool feasible, std::size_t numPoints, std::size_t numRays,
      bool lastIteration);
    virtual void onBeforeAddPoint();
    virtual void onAfterAddPoint();
    virtual void onBeforeAddRay();
    virtual void onAfterAddRay();

  private:
    bool addPointsAndRays(VectorSubset& pointIndices, VectorSubset& rayIndices, bool stabLP);
    void addPointRow(std::size_t index, bool stabLP);
    void addRayRow(std::size_t index, bool stabLP);
    void searchCacheApproximate(VectorSubset& indices, UniqueRationalVectorsBase& objects, bool points,
        const soplex::VectorReal& approxObjective, double approxRhs, std::size_t maxAdd, double epsilon);
    void searchCache(VectorSubset& indices, UniqueRationalVectorsBase& objects, bool points,
        const soplex::VectorReal& approxObjective, const soplex::VectorRational& exactObjective,
        const soplex::Rational& rhs, std::size_t maxAdd);
    void maximizeOracle(VectorSubset& pointIndices, VectorSubset& directionIndices,
        const soplex::VectorRational& exactObjective, const soplex::Rational& rhs);

  protected:
    UniqueRationalVectorsBase& _points;
    UniqueRationalVectorsBase& _directions;
    OracleBase* _oracle;

//    int iteration;

  protected:
//   private: TODO: for debugging only!
    struct RowInfo
    {
      char type;
      std::size_t index;
      int age;
    };

    const std::size_t _n;
    const std::size_t _d;
    const std::size_t _offsetLower;
    const std::size_t _offsetUpper;
    std::size_t _offsetStabilizationRows;
    bool _stabilizing;

    /// Main LP

    soplex::SoPlex* _mainLP;
    int _maxAge;
    std::vector<RowInfo> _rowInfos;
    std::vector<std::size_t> _constraintsToRows;
    soplex::DVectorRational _lastPrimalSolution;
    soplex::DVectorRational _currentPrimalSolution;

    /// Stabilization

    struct StabilizationInfo
    {
      double history[IPO_POLAR_STABIL_HIST_SIZE]; // History of last entries.
      bool active; // True if history contains distinct values.
      double trustLower; // Lower bound of trust region.
      double trustUpper; // Upper bound of trust region.
      std::size_t lowerRow;
      std::size_t upperRow;
    };
    double _initialPenalty;
    double _stabPenalty;

    soplex::SoPlex* _stabLP;
    std::vector<RowInfo> _stabRowInfos;
    std::vector<StabilizationInfo> _stabColInfos;

    /// Oracle

    OracleResult _result;

    /// Query data

    double _lastMainObjective;
    double _lastPenaltyCosts;
  };

} /* namespace ipo */

#endif /* IPO_POLAR_LP_H_ */

