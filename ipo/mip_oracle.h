#ifndef IPO_MIP_ORACLE_H_
#define IPO_MIP_ORACLE_H_

#include "common.h"
#include "rational.h"

#include <string>
#include <vector>

#ifdef WITH_SCIP
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
#endif
#endif

#include "oracles.h"
#include "lp.h"

namespace ipo {

  /**
   * \brief Base oracle for oracles based on approximate mixed-integer-programming solvers.
   *
   * Base oracle for oracles based on approximate mixed-integer-programming solvers.
   * It postprocesses the returned floating-point solutions to get feasible rational ones.
   */

  class MIPOracleBase: public OracleBase
  {
  public:
    inline const std::shared_ptr<MixedIntegerLinearSet>& mixedIntegerLinearSet() const
    {
      return _mixedIntegerLinearSet;
    }

  protected:
    /**
     * \brief Constructs the oracle.
     *
     * Constructs the oracle based on a MixedIntegerSet that is passed to the initialize() method. The latter need not be
     * complete, i.e., there may be inequalities missing. In this case, the separate() method should be implemented by the
     * inheriting class, which is then queried with a potential solution and must produce additional inequalities. Note that this
     * is only required to complete the continuous part of a solution, i.e., the integer variables of the solution will remain
     * fixed until the completion is finished. An actual MIP solver oracle is provided by inheriting from this class.
     *
     * \param name            Name of the oracle.
     * \param nextOracle      Next associated oracle.
     */

    MIPOracleBase(const std::string& name, const std::shared_ptr<OracleBase>& nextOracle = NULL);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~MIPOracleBase();

    /**
     * \brief Initializes the oracle for the given \p mixedIntegerLinearSet.
     *
     * Initializes the oracle for the given \p mixedIntegerLinearSet. Also calls OracleBase::initializeSpace().
     */

    void initialize(const std::shared_ptr<MixedIntegerLinearSet>& mixedIntegerLinearSet);

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation adds a corresponding equation to the LP that is used for the
     * postprocessing.
     */

    virtual void setFace(const LinearConstraint& newFace = completeFaceConstraint());

  protected:

    /**
     * \brief Oracle's implementation to maximize the dense rational \p objective.
     *
     * This method is called by maximizeController() and contains the implementation of the oracle.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param sort           Set this variable to true if points must be sorted.
     * \param checkDups      Set this variable to true if points or rays must be checked for duplicates.
     *
     * This implementation calls the virtual method solverMaximize() which is supposed to return a set of floating-point
     * solutions, which are then postprocessed.
     */

    virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
      bool& checkDups);

    /**
     * \brief Solver-specific maximization method.
     *
     * Solver-specific maximization method.
     *
     * \param objective Objective vector to be maximized.
     * \param points    Returned set of points (each entry is an array of length \c n \c that is free'd by the caller.
     * \param rays      Returned set of rays (each entry is an array of length \c n \c that is free'd by the caller.
     * \param hitLimit  Set to true iff solver reached some limit, e.g., due to time or memory constraints.
     */

    virtual void solverMaximize(double* objective, double objectiveBound, std::vector<double*>& points,
      std::vector<double*>& rays, bool& hitLimit) = 0;

    /**
     * \brief Method that can add lazy inequalities that cut off a given point.
     *
     * Method that can add lazy inequalities that cut off a given \p point. If not implemented, it is assumed that the MIP is
     * exactly the one passed to the constructor of the class.
     */

    virtual void separatePoint(const Vector& point, std::vector<LinearConstraint>& cuts);

    /**
     * \brief Method that can add lazy inequalities that cut off a given ray.
     *
     * Method that can add lazy inequalities that cut off a given \p ray. If not implemented, it is assumed that the MIP's
     * unbounded rays are exactly the ones induced by the MIP that was passed to the constructor of the class.
     */

    virtual void separateRay(const Vector& ray, std::vector<LinearConstraint>& cuts);

  private:

    void prepareSolver(const soplex::VectorRational& objective);

    void restoreSolver();

    Vector extendPoint(double* point, soplex::Rational& objectiveValue, bool& isRay);

    Vector computeRay();

    std::shared_ptr<MixedIntegerLinearSet> _mixedIntegerLinearSet;
    std::shared_ptr<LinearProgram> _correctionLP;
    double* _objective; // Objective to be passed to the underlying floating-pointer based LP solver.
    std::vector<LinearConstraint> _cuts;
    std::vector<double*> _points; // Array of points returned by the underlying LP solver.
    std::vector<double*> _rays; // Array of rays returned by the underlying LP solver.
  };

} /* namespace ipo */

#endif /* IPO_MIP_ORACLE_H_ */
