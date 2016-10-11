#ifndef IPO_SCIP_ORACLES_H_
#define IPO_SCIP_ORACLES_H_

#include <set>
#include <map>
#include <limits>

#include "common.h"
#include "oracles.h"
#include "mip.h"

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
#endif

namespace ipo {

  class MIPOracleBase;

  /**
   * A map to associate SCIP variables with variable indices.
   */

  typedef std::map<SCIP_VAR*, std::size_t> SCIPvarToIndexMap;

  /**
   * \brief Returns a canonical mapping from SCIP variables to variable indices.
   *
   * Returns a canonical mapping from SCIP variables to variable indices.
   * It maps original variables in the order SCIP returns them.
   * If the SCIP instance was transformed already, then
   * it also maps the corresponding transformed variables to the same indices.
   *
   * \param originalSCIP
   *   SCIP instance
   * \param map
   *   Variable map the method writes into.
   */

  void getSCIPvarToIndexMap(SCIP* originalSCIP, SCIPvarToIndexMap& map);

  /**
   * \brief Extracts the objective from a SCIP instance.
   *
   * Extracts the objective from the given SCIP instance. The variables are ordered canonically (see \ref getSCIPvarToIndexMap()).
   * By default, the objective is scaled such that it corresponds to a maximization problem.
   *
   * \param scip
   *   SCIP instance
   * \param makeMaximization
   *   If \c true, the objective is scaled such that it
   *   corresponds to a maximization problem.
   */

  Vector getSCIPObjective(SCIP* scip, bool makeMaximization = true);

  /**
   * \brief An oracle based on the SCIP solver.
   *
   * An oracle for the convex hull of the solutions returned by the SCIP instance. The computed floating-point solutions are 
   * turned into rational ones by the underlying \ref MIPOracleBase.
   *
   * An instance is either constructed from a \c SCIP instance or from a \ref MixedIntegerSet.
   */

  class SCIPOracle: public MIPOracleBase
  {
  public:
    /**
     * \brief Constructs a SCIP oracle with given \p name, optionally associated to \p nextOracle.
     *
     * Constructs a SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is
     * defined via the \p originalSCIP instance (and must be equal to that of \p nextOracle). The oracle is implemented by 
     * calling SCIP on a copy of the given \p originalSCIP instance.
     */

    SCIPOracle(const std::string& name, SCIP* originalSCIP, const std::shared_ptr<OracleBase>& nextOracle = NULL);

    /**
     * \brief Constructs a SCIP oracle with given \p name in given \p space.
     *
     * Constructs a SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is equal to
     * that of \p nextOracle and of the space of the given \p mixedIntegerSet. The oracle is implemented by calling SCIP in order 
     * to solve mixed-integer programs over the \p mixedIntegerSet.
     */

    SCIPOracle(const std::string& name, const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet, 
      const std::shared_ptr<OracleBase>& nextOracle = NULL);

    /**
     * \brief Destructor.
     */

    virtual ~SCIPOracle();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation adds an equation to the underlying SCIP instance.
     */

    virtual void setFace(const LinearConstraint& newFace = completeFace());

  protected:

    std::shared_ptr<MixedIntegerSet> constructFromSCIP(SCIP* originalSCIP);

    void constructFromMixedIntegerSet(const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet);

    virtual void solverMaximize(double* objective, double objectiveBound, std::vector<double*>& points,
      std::vector<double*>& rays);

  protected:
    SCIP* _scip; // SCIP instance
    std::vector<SCIP_VAR*> _variables; // SCIP variables.
    SCIP_CONS* _faceConstraint; // Special equation constraint for optimizing over a face.
  };

  /**
   * \brief An oracle based on the exact IP solver scip-ex.
   *
   * An oracle for the given \ref MixedIntegerSet that uses scip-ex (via an external call) to optimize.
   */
  
  class ExactSCIPOracle : public OracleBase
  {
  public:
    /**
     * \brief Constructs an exact SCIP oracle with given \p name in given \p space.
     *
     * Constructs an exact SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is 
     * equal to that of \p nextOracle and of the space of the given \p mixedIntegerSet. The oracle is implemented by calling
     * scip-ex to solve mixed-integer programs over the \p mixedIntegerSet.
     */

    ExactSCIPOracle(const std::string& name, const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet, 
      const std::shared_ptr<OracleBase>& nextOracle = NULL);

    /**
     * \brief Destructor.
     */

    virtual ~ExactSCIPOracle();

    /**
     * \brief Sets the path of the scip-ex binary used.
     * 
     * Sets the path of the scip-ex binary used.
     */

    void setBinaryPath(const std::string& path);

    /**
     * \brief Returns the path of the scip-ex binary used.
     * 
     * Returns the path of the scip-ex binary used.
     */

    inline const std::string& binaryPath() const
    {
      return _binary;
    }

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation adds an equation to the underlying SCIP instance.
     */

    virtual void setFace(const LinearConstraint& newFace = completeFace());

  protected:
    
    /**
     * \brief Creates a temporary directory to work in.
     * 
     * Creates a temporary directory to work in.
     */
    
    void createWorkingDirectory();
    
    /**
     * \brief Removes the temporary directory.
     * 
     * Removes the temporary directory created by createWorkingDirectory().
     */
    
    void deleteWorkingDirectory();

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
     * This implementation creates a ZIMPL file containing the model and calls the scip-ex binary to solve it.
     */

    virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
      bool& checkDups);
    
    void writeModel(const soplex::VectorRational& objective);
    
    void solveModel();
    
    VectorData* parseOutput();

  protected:
    std::string _binary;
    std::string _workingDirectory;
    std::shared_ptr<MixedIntegerSet> _mixedIntegerSet;
  };

//   class ExactSCIPOptimizationOracle: public FaceOptimizationOracleBase
//   {
//   public:
//
//     /**
//      * \brief Constructor.
//      *
//      * Constructs an oracle for the given \ref MIP.
//      * It calls ExactSCIP externally and reads back the optimal
//      * solution from a file.
//      * For this the constructor creates a temporary directory.
//      *
//      * \param name
//      *   Name of the oracle.
//      * \param exactBinary
//      *   Binary of ExactSCIP.
//      * \param mip
//      *   Mixed-integer program to optimize over.
//      * \param heuristic
//      *   If not \c NULL, this other oracle is called if optimality is not required.
//      * \param timeLimit
//      *   Time limit (in seconds) for the solver.
//      */
//
//     ExactSCIPOptimizationOracle(const std::string& name, const std::string& exactBinary, MixedIntegerProgram& mip,
//         FaceOptimizationOracleBase* heuristic, double timeLimit = std::numeric_limits<double>::max());
//
//     /**
//      * \brief Destructor.
//      *
//      * Destructor.
//      * It deletes the temporary directory.
//      */
//
//     virtual ~ExactSCIPOptimizationOracle();
//
//     /**
//      * Whether to scale the objective before passing it to ExactSCIP.
//      */
//
//     bool scaleObjective;
//
//   protected:
//
//     /**
//      * \brief Actual implementation.
//      *
//      * The actual implementation.
//      * It writes the MIP with the current objective vector to
//      * a ZIMPL model and calls the ExactSCIP solver as an external program,
//      * redirecting the output into the temporary directory.
//      */
//
//     virtual void run(OracleResult& result, const soplex::VectorRational& objective,
//         const soplex::Rational* improveValue, bool forceOptimal);
//
//     /**
//      * \brief Method that is called when a new face is activated.
//      *
//      * Method that is called when a new \c face is activated.
//      * If a (non-trivial) face was active before,
//      * it is ensured that \ref faceDisabled() is called before.
//      *
//      * The implementation adds an equation constraint
//      * to the underlying mixed-integer program.
//      */
//
//     virtual void faceEnabled(Face* face);
//
//     /**
//      * \brief Method that is called when a face is deactivated.
//      *
//      * Method that is called when a face is deactivated.
//      *
//      * The implementation removes the constraint that was added in \ref faceEnabled().
//      */
//
//     virtual void faceDisabled(Face* face);
//
//     /**
//      * \brief Creates a temporary directory.
//      *
//      * Creates a temporary directory that contains
//      * the model file and logs.
//      */
//
//     void createTempDirectory();
//
//     /**
//      * \brief Removes the temporary directory.
//      *
//      * Removes the temporary directory created in \ref createTempDirectory().
//      */
//
//     void deleteTempDirectory();
//
//     /**
//      * \brief Writes a ZIMPL model into the temporary directory.
//      *
//      * Writes a ZIMPL model into the temporary directory.
//      */
//
//     void writeModel(const soplex::VectorRational& objective);
//
//     /**
//      * Calls the ExactSCIP solver.
//      */
//
//     void callSolver();
//
//     /**
//      * \brief Parses ExactSCIP's output.
//      *
//      * Parses ExactSCIP's output and extracts the solution vector.
//      */
//
//     soplex::DSVectorRational* parseOutput();
//
//   protected:
//
//     const std::string _binary; // The ExactSCIP binary.
//     MixedIntegerProgram& _mip; // Reference to the underlying MIP.
//     FaceOptimizationOracleBase* _heuristic; // Optional heuristic oracle.
//     std::string _path; // Path to the temporary directory.
//     double _timeLimit; // Time limit for the solver (in seconds).
//   };

} /* namespace ipo */

#endif /* IPO_SCIP_ORACLES_H_ */
