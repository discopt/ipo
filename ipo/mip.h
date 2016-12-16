#ifndef IPO_MIP_H_
#define IPO_MIP_H_

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

namespace ipo {

  /**
   * \brief A mixed-integer set.
   *
   * A mixed-integer set. It is typically used to create oracles for optimizing over it, e.g., using \ref SCIPOracle.
   */

  class MixedIntegerSet
  {
  public:
    struct Variable
    {
      bool integral;
      Rational upperBound;
      Rational lowerBound;

      Variable();
      Variable(const Rational& lowerBound, const Rational& upperBound, bool integral);
    };

  public:
#ifdef IPO_WITH_SCIP
    /**
     * \brief Constructs a \c MixedIntegerSet from a \c SCIP instance.
     *
     * Constructs a \c MixedIntegerSet from a \c SCIP instance. Only explicitly stated linear constraints of the instance are
     * considered.
     */

    MixedIntegerSet(SCIP* scip);
#endif /* IPO_WITH_SCIP */

    /**
     * \brief Destructor.
     */

    virtual ~MixedIntegerSet();

    /**
     * \brief Returns the space.
     *
     * Returns a const-reference to the space.
     */

    inline const Space& space() const
    {
      return _space;
    }

    /**
     * \brief Returns the number of columns.
     *
     * Returns the number of columns.
     */

    inline std::size_t numVariables() const
    {
      return _space.dimension();
    }

    /**
     * \brief Returns the number of row constraints.
     *
     * Returns the number of row constraints.
     */

    inline std::size_t numRows() const
    {
      return _rowNames.size();
    }

    /**
     * \brief Returns a \p variable.
     *
     * Returns a const-reference to the specified \p variable.
     */

    inline const Variable& variable(std::size_t variable) const
    {
      return _variables[variable];
    }

    /**
     * \brief Returns a \p row constraint.
     *
     * Returns a const-reference to a \p row constraint.
     */

    inline const LinearConstraint& rowConstraint(std::size_t row) const
    {
      return _rowConstraints[row];
    }

    /**
     * \brief Returns all row constraints.
     *
     * Returns a const-reference to vector containing all row constraint.
     */

    inline const std::vector<LinearConstraint>& rowConstraints() const
    {
      return _rowConstraints;
    }

    /**
     * \brief Returns a \p row name.
     *
     * Returns a const-refrence to a \p row name.
     */

    inline const std::string& rowName(std::size_t row) const
    {
      return _rowNames[row];
    }

    /**
     * \brief Returns whether the given variable is integral.
     *
     * Returns whether the given \c variable has an integrality constraint.
     */

    inline bool isIntegral(std::size_t variable) const
    {
      return _variables[variable].integral;
    }

//     /**
//      * \brief Returns true if the given \c point satisfies all bound constraints.
//      *
//      * Returns true if the given \c point satisfies all bound constraints.
//      */
//
//     bool checkPointBounds(const soplex::SVectorRational* point) const;
//
//     /**
//      * \brief Returns true if the given \c point satisfies all row constraints.
//      *
//      * Returns true if the given \c point satisfies all row constraints.
//      */
//
//     bool checkPointRows(const soplex::SVectorRational* point);
//
//     /**
//      * \brief Returns true if the given \c point satisfies all integrality constraints.
//      *
//      * Returns true if the given \c point satisfies all integrality constraints.
//      */
//
//     bool checkPointIntegral(const soplex::SVectorRational* point) const;
//
//     /**
//      * \brief Returns true if the given \c point is feasible.
//      *
//      * Returns true if the given \c point is feasible.
//      */
//
//     bool checkPoint(const soplex::SVectorRational* point);
//
//     /**
//      * \brief Returns true if the given \c ray satisfies all bound constraints.
//      *
//      * Returns true if the given \c ray satisfies all bound constraints.
//      */
//
//     bool checkRayBounds(const soplex::SVectorRational* ray) const;
//
//     /**
//      * \brief Returns true if the given \c ray satisfies all row constraints.
//      *
//      * Returns true if the given \c ray satisfies all row constraints.
//      */
//
//     bool checkRayRows(const soplex::SVectorRational* ray);
//
//     /**
//      * \brief Returns true if the given \c ray is feasible.
//      *
//      * Returns true if the given \c ray is feasible.
//      */
//
//     bool checkRay(const soplex::SVectorRational* ray);
//
    /**
      * \brief Restricts the MixedIntegerSet to the given face.
      *
      * Restricts the MixedIntegerSet to the given face by adding an equation constraint.
      */

    void setFace(const LinearConstraint& newFace);

//     /**
//      * \brief Returns all row constraints.
//      *
//      * Returns all row constraints or only those that are inequalities or equations,
//      * respectively.
//      *
//      * \param rows
//      *   Rows structure to write to.
//      * \param inequalities
//      *   Whether to extract inequalities.
//      * \param equations
//      *   Whether to extract equations.
//      * \param names
//      *   If not \c NULL, writes corresponding row names.
//      */
//
//     void getConstraints(soplex::LPRowSetRational& rows, bool inequalities, bool equations,
//         std::vector<std::string>* names = NULL);
//
//     /**
//      * \brief Extracts equations that correspond to fixed variables.
//      *
//      * Extracts equations that correspond to fixed variables.
//      *
//      * \param rows
//      *   Rows structure to write to.
//      * \param names
//      *   If not \c NULL, writes corresponding row names.
//      */
//
//     void getFixedVariableEquations(soplex::LPRowSetRational& rows,
//       std::vector<std::string>* names = NULL);

  protected:
    Space _space; // Space with column names.
    std::vector<Variable> _variables;
    std::vector<LinearConstraint> _rowConstraints; // Row constraints.
    std::vector<std::string> _rowNames; // Row constraints' names.
    LinearConstraint _currentFace; // Currently active face constraint.
  };

  /**
   * \brief Base oracle for oracles based on approximate mixed-integer-programming solvers.
   *
   * Base oracle for oracles based on approximate mixed-integer-programming solvers.
   * It postprocesses the returned floating-point solutions to get feasible rational ones.
   */

  class MIPOracleBase: public OracleBase
  {
  public:
    inline const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet() const
    {
      return _mixedIntegerSet;
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
     * \brief Initializes the oracle for the given \p mixedIntegerSet.
     *
     * Initializes the oracle for the given \p mixedIntegerSet. Also calls OracleBase::initializeSpace().
     */

    void initialize(const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet);

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

//     /**
//      * \brief Initializes the LP and solver.
//      *
//      * Initializes the LP and solver.
//      */
//
//     void initializeLP();
//
//     /**
//      * \brief Corrects a given point.
//      *
//      * Corrects a given point by rounding integer variables
//      * and solving an LP for the continuous ones.
//      */
//
//     soplex::DSVectorRational* correctPoint(const soplex::SVectorRational* point,
//         const soplex::VectorRational& objective);
//
//     /**
//      * \brief Corrects a given ray.
//      *
//      * Corrects a given unbounded ray
//      * by solving an LP.
//      */
//
//     soplex::DSVectorRational* correctDirection(const soplex::SVectorRational* direction,
//       const soplex::VectorRational& objective);
//
//     /**
//      * \brief Oracle's implementation to maximize the dense rational \p objective.
//      *
//      * This method is called by maximizeController() and contains the implementation of the oracle.
//      *
//      *
//      * \param result         After the call, contains the oracle's answer.
//      * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
//      * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
//      * \param sort           Set this variable to true if points must be sorted.
//      * \param checkDups      Set this variable to true if points or rays must be checked for duplicates.
//      *
//      * This implementation
//      * For requirements on the behavior, see Detailed Description of \ref OracleBase.
//      */
//
//     virtual std::size_t maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
//       const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
// = 0;


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
     */

    virtual void solverMaximize(double* objective, double objectiveBound, std::vector<double*>& points,
      std::vector<double*>& rays) = 0;

    /**
     * \brief Method that can add lazy inequalities that cut off a given point.
     *
     * Method that can add lazy inequalities that cut off a given \p point. If not implemented, it is assumed that the MIP is
     * exactly the one passed to the constructor of the class.
     */

    virtual void separatePoint(const soplex::VectorRational& point, soplex::LPRowSetRational& cuts);

    /**
     * \brief Method that can add lazy inequalities that cut off a given ray.
     *
     * Method that can add lazy inequalities that cut off a given \p ray. If not implemented, it is assumed that the MIP's
     * unbounded rays are exactly the ones induced by the MIP that was passed to the constructor of the class.
     */

    virtual void separateRay(const soplex::VectorRational& ray, soplex::LPRowSetRational& cuts);

  private:

    void prepareSolver(const soplex::VectorRational& objective);

    void restoreSolver();

    Vector extendPoint(double* point, soplex::Rational& objectiveValue);

    Vector computeRay();

    std::shared_ptr<MixedIntegerSet> _mixedIntegerSet;
    soplex::SoPlex* _spx; // LP solver with the correction LP.
    std::size_t _numRows;
    double* _objective;
    soplex::DVectorRational _lpResult;
    std::vector<int> _lpRowPermutation;
    soplex::LPRowSetRational _separateResult;
    std::vector<double*> _points;
    std::vector<double*> _rays;
  };

} /* namespace ipo */

#endif /* IPO_MIP_H_ */
