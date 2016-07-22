#ifndef IPO_SCIP_ORACLES_H_
#define IPO_SCIP_ORACLES_H_

#include <set>
#include <map>
#include <limits>

#include "common.h"
#include "oracles.h"
#include "mixed_integer_program.h"

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
#endif


namespace ipo {

  class MixedIntegerProgram;

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
   * Extracts the objective from the given SCIP instance
   * as a dense rational vector. The variables
   * are ordered canonically (see \ref getSCIPvarToIndexMap()).
   * By default, the objective is scaled such that
   * it corresponds to a maximization problem.
   *
   * \param originalSCIP
   *   SCIP instance
   * \param objective
   *   Objective vector for the objective. Its dimension is set appropriately.
   * \param makeMaximization
   *   If \c true, the objective is scaled such that it
   *   corresponds to a maximization problem.
   */

  void getSCIPObjective(SCIP* originalSCIP, soplex::DVectorRational& objective, bool makeMaximization = true);

  /**
   * \brief An oracle based on the SCIP solver.
   *
   * An oracle for the convex hull of the solutions returned by the SCIP instance. The computed
   * floating-point solutions are turned into rational ones using rational-reconstruction
   * techniques. This may in principle lead to infeasible solutions, in particular if continuous
   * solutions are present. To produce correct solutions one may use a
   * \ref MixedIntegerProgramCorrectorOracle that postprocesses the solutions.
   *
   * An instance is either constructed from a \c SCIP instance
   * or from a \ref MixedIntegerProgram.
   */

  class SCIPOracle: public OracleBase
  {
  public:
     /**
     * \brief Constructs a SCIP oracle with given \p name in given \p space.
     *
     * Constructs a SCIP oracle with given \p name in given \p space. The oracle is implemented by
     * calling SCIP on a copy of the given \p originalSCIP instance.
     */

    SCIPOracle(const std::string& name, const Space& space, SCIP* originalSCIP);

     /**
     * \brief Constructs a SCIP oracle with given \p name in given \p space.
     *
     * Constructs a SCIP oracle with given \p name in given \p space. The oracle is implemented by
     * calling SCIP on a copy of the given \p originalSCIP instance. If \p space is 0-dimensional,
     * it is augmented according to the variables of the SCIP instance.
     */

    SCIPOracle(const std::string& name, Space& space, SCIP* originalSCIP);

    /**
     * \brief Constructs a heuristic SCIP oracle with given \p name associated to \p nextOracle.
     *
     * Constructs a heuristic SCIP oracle with given \p name that is associated to \p nextOracle.
     * The ambient space is equal to that of \p nextOracle. The oracle is implemented by calling
     * SCIP on a copy of the given \p originalSCIP instance.
     */

    SCIPOracle(const std::string& name, OracleBase* nextOracle, SCIP* originalSCIP);

    /**
     * \brief Constructs a SCIP oracle with given \p name in given \p space.
     *
     * Constructs a SCIP oracle with given \p name. The oracle is implemented by calling SCIP in
     * order to solve the mixed-integer program \p mip. The ambient space is equal to that of
     * \p mip.
     */

    SCIPOracle(const std::string& name, const MixedIntegerProgram& mip);

    /**
     * \brief Constructs a SCIP oracle with given \p name in given \p space.
     *
     * Constructs a SCIP oracle with given \p name that is associated to \p nextOracle. The oracle
     * is implemented by calling SCIP in order to solve the mixed-integer program \p mip. The
     * ambient space is equal to that of \p mip and to that of \p nextOracle.
     */

    SCIPOracle(const std::string& name, OracleBase* nextOracle, const MixedIntegerProgram& mip);

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

    virtual void setFace(Face* newFace = NULL);

    /**
     * \brief Runs the oracle to maximize the dense rational \p objective.
     *
     * Runs the optimization oracle to maximize the given dense rational \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p objectiveBound is not exceeded, then the call must be forwarded to the
     * next oracle.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * This implementation calls SCIP and reconstructs rational solutions from the returned
     * floating-point ones.
     */

    virtual void maximize(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

  protected:

    void copySCIP(SCIP* originalSCIP, Space& scipSpace);

    void initializeFromMIP(const MixedIntegerProgram& mip);

    /**
     * \brief Computes a scaling factor for the objective vector.
     *
     * Computes a scaling factor for the objective vector.
     *
     * \sa run()
     */

    soplex::Rational computeVectorScaling(const soplex::VectorRational& vector);

  protected:
    SCIP* _scip; // SCIP instance
    soplex::Rational _maxLargestCoefficient; // Maximum largest coefficient we hand over to SCIP.
    soplex::Rational _minLargestCoefficient; // Minimum largest coefficient we hand over to SCIP.
    soplex::Rational _bestLargestCoefficient; // Best largest coefficient we hand over to SCIP.
    std::vector<SCIP_VAR*> _variables; // SCIP variables.
    SCIP_CONS* _faceConstraint; // Special equation constraint for optimizing over a face.
  };

//   /**
//    * \brief An oracle based on ExactSCIP.
//    *
//    * A oracle for the polyhedron \f$ P \f$
//    * that uses SCIP to optimize.
//    * An instance is constructed from a \ref MixedIntegerProgram.
//    */
//
//   class ExactSCIPOptimizationOracle: public FaceOptimizationOracleBase
//   {
//   public:
//
//     /**
//      * \brief Constructor.
//      *
//      * Constructs an oracle for the given \ref MixedIntegerProgram.
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
