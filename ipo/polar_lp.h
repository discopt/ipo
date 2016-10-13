#ifndef IPO_POLAR_LP_H_
#define IPO_POLAR_LP_H_

#include <vector>
#include <set>
#include <map>

#include "common.h"
#include "vectors.h"
#include "oracles.h"

#define IPO_POLAR_STABIL_HIST_SIZE 2

#include "soplex_reproduce.h"

namespace ipo {

  class XPolarLP;

  /**
   * \brief Base class for an observer for polar-LP computations.
   *
   * Base class for an observer for polar-LP computations.
   */

  class PolarLPHandler
  {
  public:
    enum Event
    {
      LP_BEGIN,
      LP_END,
      ORACLE_BEGIN,
      ORACLE_END = ORACLE_BEGIN + 1,
      POINTS_BEGIN,
      POINT,
      POINTS_END,
      RAYS_BEGIN,
      RAY,
      RAYS_END,
      OBJECTIVE_SET,
      ROW_ADDED,
      ROW_UPDATED,
      SOLVE_BEGIN,
      SOLVE_END,
    };

    /**
     * \brief Default constructor.
     *
     * Default constructor.
     */

    PolarLPHandler();

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~PolarLPHandler();

    /**
     * \brief This method is called by the algorithm.
     *
     * This method is called by the algorithm in certain steps.
     */

    virtual void notify(Event event, XPolarLP& polarLP) = 0;
  };

  class XPolarLP
  {
  protected:

    struct RowInfo
    {
      bool dynamic;
      char type;
      Vector vector;
      int age;

      RowInfo(bool dynamic);
//       RowInfo(char type, int age = -1);
      RowInfo(bool dynamic, char type, const Vector& vector, int age = -1);
    };

  public:

    XPolarLP(const std::shared_ptr<OracleBase>& oracle, PolarLPHandler& handler, bool approximate);

    ~XPolarLP();

    void clear();

    void setObjective(const soplex::VectorRational& objective);

    std::size_t addRow(const Rational& lhs, const soplex::SVectorRational& normalVector, const Rational& rhs, bool dynamic);

    void updateRow(std::size_t row, const Rational& lhs, const soplex::SVectorRational& normalVector, const Rational& rhs);

    void addPointRow(const Vector& point, bool dynamic);

    void addRayRow(const Vector& ray, bool dynamic);

    void solve();

    inline Rational getObjectiveValue() const
    {
      return _spx->objValueRational();
    }

    LinearConstraint currentInequality() const;

    inline soplex::Rational oracleObjectiveValue() const
    {
      return _oracleObjectiveValue;
    }

    void getTightPointsRays(InnerDescription& tightPointsRays, bool dynamicOnly = false);

    /**
     * \brief Returns the space of the polar LP.
     *
     * Returns a const-reference to the space of the polar LP.
     */

    inline const Space& polarSpace() const
    {
      return *_space;
    }

    /**
      * \brief Returns the feasibility / optimality tolerance of the LP.
      *
      * Returns the feasibility / optimality tolerance of the LP.
      */

    inline double getTolerance() const
    {
      return _spx->realParam(soplex::SoPlex::FEASTOL);
    }

    /**
     * \brief Returns the number of points that are currently in the LP.
     *
     * Returns the number of points that are currently in the LP.
     */

    inline std::size_t numPointsLP() const
    {
      return _numPointsLP;
    }

    /**
     * \brief Returns the number of rays that are currently in the LP.
     *
     * Returns the number of rays found that are currently in the LP.
     */

    inline std::size_t numRaysLP() const
    {
      return _numRaysLP;
    }

    /**
     * \brief Returns the number of rows of the current LP.
     *
     * Returns the number of rows of the current LP.
     */

    inline std::size_t numRowsLP() const
    {
      return _spx->numRowsRational();
    }

    /**
     * \brief Returns the number of columns of the current LP.
     *
     * Returns the number of columns of the current LP.
     */

    inline std::size_t numColumnsLP() const
    {
      return _spx->numColsRational();
    }

    /**
     * \brief Returns the number of nonzeros of the current LP.
     *
     * Returns the number of nonzeros of the current LP.
     */

    inline std::size_t numNonzerosLP() const
    {
      return _spx->numNonzerosRational();
    }

    /**
     * \brief Returns the maximum allowed heuristic level of the current oracle call.
     *
     * Returns the maximum allowed heuristic level of the current oracle call.
     */

    inline HeuristicLevel oracleMaxHeuristicLevel() const
    {
      return _oracleMaxHeuristicLevel;
    }

    /**
     * \brief Returns the minimum allowed heuristic level of the current oracle call.
     *
     * Returns the minimum allowed heuristic level of the current oracle call.
     */

    inline HeuristicLevel oracleMinHeuristicLevel() const
    {
      return _oracleMinHeuristicLevel;
    }

    /**
     * \brief Returns the heuristic level of the oracle's last answer.
     *
     * Returns the heuristic level of the oracle's last answer.
     */

    inline HeuristicLevel oracleResultHeuristicLevel() const
    {
      return _result.heuristicLevel();
    }

    /**
     * \brief Returns the point or ray that is currently being added.
     *
     * Returns the point or ray that is currently being added.
     */

    inline Vector currentVector() const
    {
      return _currentVector;
    }

    /**
     * \brief Returns the number of points returned by the last oracle call.
     *
     * Returns the number of points returned by the last oracle call.
     */

    inline std::size_t oracleNumPoints() const
    {
      return _result.points.size();
    }

    /**
     * \brief Returns the number of rays returned by the last oracle call.
     *
     * Returns the number of rays returned by the last oracle call.
     */

    inline std::size_t oracleNumRays() const
    {
      return _result.rays.size();
    }

    inline bool isApproximate() const
    {
      return _spx->realParam(soplex::SoPlex::FEASTOL) != 0.0;
    }

  protected:

    void notify(PolarLPHandler::Event event);

  protected:

    bool _approximate;
    soplex::LPColSetRational _initialColumns;
    soplex::LPRowSetRational _initialRows;
    Space* _space;
    std::shared_ptr<OracleBase> _oracle;
    soplex::DVectorRational _oracleObjective;
    OracleResult _result;
    PolarLPHandler& _handler;
    soplex::SoPlex* _spx;
    soplex::DVectorRational _solution;
    soplex::DVectorRational _denseColumnVector;
    soplex::DSVectorRational _sparseColumnVector;
    std::vector<RowInfo> _rows;
    std::size_t _firstDynamicRow;

    Vector _currentVector;
    std::size_t _numPointsLP;
    std::size_t _numRaysLP;
    HeuristicLevel _oracleMaxHeuristicLevel;
    HeuristicLevel _oracleMinHeuristicLevel;
    soplex::Rational _oracleObjectiveValue;
  };


  /////////////////////// OLD INTERFACE ///////////////////////

  /// Class for solving polar LPs using lazy separation, cut aging and stabilization.

  class PolarLP
  {
  public:
    struct Basis
    {
      std::vector<soplex::SPxSolver::VarStatus> columnStatus;
      std::vector<soplex::SPxSolver::VarStatus> constraintStatus;
      std::vector<Vector> tightPoints;
      std::vector<Vector> tightRays;
    };

    PolarLP(const std::shared_ptr<OracleBase>& oracle, double initialPenalty = 1024.0, int maxAge = 30);
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
    std::size_t addPointConstraint(const Vector& point);
    std::size_t addRayConstraint(const Vector& ray);

    void setBasis(const Basis& basis);
    void getBasis(Basis& basis);
    soplex::Rational getObjectiveValue();
    void getPrimalSolution(soplex::VectorRational& solution);
    void getPrimalSolution(soplex::DSVectorRational& solution); // TODO: all versions useful?
    Vector getPrimalSolution();
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
    bool addPointsAndRays(std::vector<Vector>& pointIndices, std::vector<Vector>& rayIndices, bool stabLP);
    void addPointRow(Vector& point, bool stabLP);
    void addRayRow(Vector& ray, bool stabLP);

  protected:
    std::shared_ptr<OracleBase> _oracle;

//    int iteration;

  protected:
//   private: TODO: for debugging only!
    struct RowInfo
    {
      char type;
      Vector vector;
      int age;

      RowInfo();
      RowInfo(char type, int age = -1);
      RowInfo(char type, const Vector& vector, int age = -1);
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

