#ifndef IPO_MIXED_INTEGER_PROGRAM_H_
#define IPO_MIXED_INTEGER_PROGRAM_H_

#include "common.h"

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
   * \brief A mixed-integer program.
   *
   * A mixed-integer program.
   * It is typically used to create oracles for optimizing over it,
   * e.g., using \ref SCIPOracle and \ref ExactSCIPOracle.
   */

  class MixedIntegerProgram
  {
  public:
#ifdef WITH_SCIP
    /**
     * \brief Constructs a \c MixedIntegerProgram from a \c SCIP instance.
     *
     * Constructs a \c MixedIntegerProgram from a \c SCIP instance.
     * Only explicitly stated linear constraints of the instance are considered.
     */

    MixedIntegerProgram(Space& space, SCIP* scip);
#endif

    /**
     * \brief Destructor.
     */

    virtual ~MixedIntegerProgram();

    /**
     * \brief Returns the space.
     *
     * Returns a const reference to the space.
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

    inline std::size_t numColumns() const
    {
      return _space.dimension();
    }

    /**
     * \brief Returns the number of rows.
     *
     * Returns the number of rows.
     */

    inline std::size_t numRows() const
    {
      return _rowNames.size();
    }

    /**
     * \brief Returns the columns.
     *
     * Returns a const reference to the columns.
     */

    inline const soplex::LPColSetRational& columns() const
    {
      return _columns;
    }

    /**
     * \brief Returns the rows.
     *
     * Returns a const reference to the rows.
     */

    inline const soplex::LPRowSetRational& rows() const
    {
      return _rows;
    }

    /**
     * \brief Returns whether the given variable is integral.
     *
     * Returns whether the given \c variable has an integrality constraint.
     */

    inline bool isIntegral(std::size_t variable) const
    {
      return _integrality[variable];
    }

    /**
     * \brief Returns the name of the given row.
     *
     * Returns the name of the given \c row.
     */

    inline const std::string& rowName(std::size_t row) const
    {
      return _rowNames[row];
    }

    /**
     * \brief Returns true if the given \c point satisfies all bound constraints.
     *
     * Returns true if the given \c point satisfies all bound constraints.
     */

    bool checkPointBounds(const soplex::SVectorRational* point) const;

    /**
     * \brief Returns true if the given \c point satisfies all row constraints.
     *
     * Returns true if the given \c point satisfies all row constraints.
     */

    bool checkPointRows(const soplex::SVectorRational* point);

    /**
     * \brief Returns true if the given \c point satisfies all integrality constraints.
     *
     * Returns true if the given \c point satisfies all integrality constraints.
     */

    bool checkPointIntegral(const soplex::SVectorRational* point) const;

    /**
     * \brief Returns true if the given \c point is feasible.
     *
     * Returns true if the given \c point is feasible.
     */

    bool checkPoint(const soplex::SVectorRational* point);

    /**
     * \brief Returns true if the given \c ray satisfies all bound constraints.
     *
     * Returns true if the given \c ray satisfies all bound constraints.
     */

    bool checkRayBounds(const soplex::SVectorRational* ray) const;

    /**
     * \brief Returns true if the given \c ray satisfies all row constraints.
     *
     * Returns true if the given \c ray satisfies all row constraints.
     */

    bool checkRayRows(const soplex::SVectorRational* ray);

    /**
     * \brief Returns true if the given \c ray is feasible.
     *
     * Returns true if the given \c ray is feasible.
     */

    bool checkRay(const soplex::SVectorRational* ray);

    /**
     * \brief Restricts the MIP to the given face.
     *
     * Restricts the MIP to the given face by adding an equation constraint.
     */

    void setFace(Face* newFace);

    /**
     * \brief Returns all row constraints.
     *
     * Returns all row constraints or only those that are inequalities or equations,
     * respectively.
     *
     * \param rows
     *   Rows structure to write to.
     * \param inequalities
     *   Whether to extract inequalities.
     * \param equations
     *   Whether to extract equations.
     * \param names
     *   If not \c NULL, writes corresponding row names.
     */

    void getConstraints(soplex::LPRowSetRational& rows, bool inequalities, bool equations,
        std::vector<std::string>* names = NULL);

    /**
     * \brief Extracts equations that correspond to fixed variables.
     *
     * Extracts equations that correspond to fixed variables.
     *
     * \param rows
     *   Rows structure to write to.
     * \param names
     *   If not \c NULL, writes corresponding row names.
     */

    void getFixedVariableEquations(soplex::LPRowSetRational& rows,
      std::vector<std::string>* names = NULL);

  protected:
    const Space& _space; // Space with column names.
    soplex::LPColSetRational _columns; // Columns
    soplex::LPRowSetRational _rows; // Rows
    std::vector<std::string> _rowNames; // Row names
    std::vector<bool> _integrality; // Integrality constraints
    soplex::DVectorRational _worker; // Temporary dense rational vector.
    Face* _currentFace; // Currently active face.
  };

  /**
   * \brief Base oracle for oracles based on approximate mixed-integer-programming solvers.
   *
   * Base oracle for oracles based on approximate mixed-integer-programming solvers.
   * It postprocesses the returned floating-point solutions to get feasible rational ones.
   */

  class MIPOracleBase: public OracleBase
  {
  protected:

    /**
     * \brief Constructs the oracle.
     *
     * Constructs the oracle based on the MIP that is passed. The latter need not be complete, i.e., there may be inequalities
     * missing. In this case, the separate() method should be implemented by the inheriting class, which is then queried with
     * a potential solution and must produce additional inequalities. Not that this is only required to complete the continuous
     * part of a solution, i.e., the integer variables of the solution will remain fixed until the completion is finished.
     * The actual MIP solver oracle is provided by inheriting from this class.
     *
     * \param name              Name of the oracle.
     * \param mip               Associated mixed-integer program.
     */

    MIPOracleBase(const std::string& name, const Space& space, const MixedIntegerProgram& mip);

    /**
     * \brief Constructs the oracle.
     *
     * Constructs the oracle based on the MIP that is passed. The latter need not be complete, i.e., there may be inequalities
     * missing. In this case, the separate() method should be implemented by the inheriting class, which is then queried with
     * a potential solution and must produce additional inequalities. Not that this is only required to complete the continuous
     * part of a solution, i.e., the integer variables of the solution will remain fixed until the completion is finished.
     * The actual MIP solver oracle is provided by inheriting from this class.
     *
     * \param name              Name of the oracle.
     * \param mip               Associated mixed-integer program.
     * \param nextOracle        Next associated oracle.
     */

    MIPOracleBase(const std::string& name, OracleBase* nextOracle, const MixedIntegerProgram& mip);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~MIPOracleBase();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation adds a corresponding equation to the LP that is used for the
     * postprocessing.
     */

    virtual void setFace(Face* newFace = NULL);

  protected:
    
    struct Column
    {
      bool integral;
      soplex::Rational upper;
      soplex::Rational lower;
    };

    void initializeLP(const MixedIntegerProgram& mip);
    
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
//      * Corrects a given ray (unbounded direction)
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

    virtual std::size_t maximizeImplementation(OracleResult& result, const DenseVector& objective,
      const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);

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
     * unbounded directions are exactly the ones induced by the MIP that was passed to the constructor of the class.
     */

    virtual void separateRay(const soplex::VectorRational& ray, soplex::LPRowSetRational& cuts);

  private:
    
    void prepareSolver(const soplex::VectorRational& objective);
    
    void restoreSolver();

    SparseVector extendPoint(double* point, soplex::Rational& objectiveValue);
    
    SparseVector computeRay();

    soplex::SoPlex* _spx; // LP solver with the correction LP.
    std::vector<Column> _columns;
    std::size_t _numRows;
    double* _objective;
    soplex::DVectorRational _lpResult;
    std::vector<int> _lpRowPermutation;
    soplex::LPRowSetRational _separateResult;
    std::vector<double*> _points;
    std::vector<double*> _rays;

//     MixedIntegerProgram& _mip; // Associated mixed-integer program.
//     OracleBase* _approximateOracle; // Approximate oracle.
//     soplex::DVectorRational _denseVector; // Temporary dense rational vector.
  };

} /* namespace ipo */

#endif /* IPO_MIXED_INTEGER_PROGRAM_H_ */
