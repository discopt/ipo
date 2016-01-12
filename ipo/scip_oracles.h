#ifndef IPO_SCIP_ORACLES_H_
#define IPO_SCIP_ORACLES_H_

#include <set>
#include <map>
#include <limits>

#include <scip/scip.h>

#include "ipo.h"
#include "oracles.h"
#include "mixed_integer_program.h"

namespace ipo {

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
   * \brief A (heuristic) oracle based on SCIP.
   *
   * A (heuristic) oracle for the polyhedron \f$ P \f$
   * that uses SCIP to optimize.
   * The computed floating-point solutions are turned into
   * rational ones using rational-reconstruction techniques.
   * This may in principle lead to infeasible solutions,
   * in particular if continuous solutions are present.
   * To produce correct solutions one may use a
   * \ref MixedIntegerProgramCorrectorOracle
   * that postprocesses the solutions.
   *
   * An instance is either constructed from a \c SCIP instance
   * or from a \ref MixedIntegerProgram.
   */

  class SCIPOptimizationOracle: public FaceOptimizationOracleBase
  {
  public:

    /**
     * \brief Constructor based on a SCIP instance.
     *
     * Constructs the oracle from a SCIP instance
     * which is copied, that is, the instance may be
     * freed afterwards.
     *
     * \param name
     *   Name of the oracle.
     * \param originalSCIP
     *   SCIP instance that defines \f$ P \f$.
     * \param isHeuristic
     *   If set to \c false, then the oracle will pretend that returned solutions are optimal.
     */

    SCIPOptimizationOracle(const std::string& name, SCIP* originalSCIP, bool isHeuristic = true);

    /**
     * \brief Constructor based on an explicit MIP.
     *
     * Constructs the oracle by creating a SCIP
     * instance for the given \ref MixedIntegerProgram.
     *
     * \param name
     *   Name of the oracle.
     * \param mip
     *   Mixed-integer program to optimize over.
     * \param isHeuristic
     *   If set to \c false, then the oracle will pretend that returned solutions are optimal.
     */

    SCIPOptimizationOracle(const std::string& name, const MixedIntegerProgram& mip, bool isHeuristic = true);

    /**
     * \brief Destructor.
     */

    virtual ~SCIPOptimizationOracle();

  protected:

    /**
     * \brief Actual implementation.
     *
     * The actual implementation.
     * It scales the objective to be suitable for a floating-point solver
     * and uses SCIP to optimize.
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    /**
     * \brief Method that is called when a new face is activated.
     *
     * Method that is called when a new \c face is activated.
     * If a (non-trivial) face was active before,
     * it is ensured that \ref faceDisabled() is called before.
     *
     * The implementation adds an equation constraint.
     */

    virtual void faceEnabled(Face* face);

    /**
     * \brief Method that is called when a face is deactivated.
     *
     * Method that is called when a face is deactivated.
     *
     * The implementation removes the constraint that was added in \ref faceEnabled().
     */

    virtual void faceDisabled(Face* face);

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
    bool _isHeuristic; // If \c false, it pretends that the returned solutions are optimal.
    soplex::Rational _maxLargestCoefficient; // Maximum largest coefficient we hand over to SCIP.
    soplex::Rational _minLargestCoefficient; // Minimum largest coefficient we hand over to SCIP.
    soplex::Rational _bestLargestCoefficient; // Best largest coefficient we hand over to SCIP.
    std::vector<SCIP_VAR*> _variables; // SCIP variables.
    SCIP_CONS* _faceConstraint; // Special equation constraint for optimizing over a face.
  };

  /**
   * \brief An oracle based on ExactSCIP.
   *
   * A oracle for the polyhedron \f$ P \f$
   * that uses SCIP to optimize.
   * An instance is constructed from a \ref MixedIntegerProgram.
   */

  class ExactSCIPOptimizationOracle: public FaceOptimizationOracleBase
  {
  public:

    /**
     * \brief Constructor.
     *
     * Constructs an oracle for the given \ref MixedIntegerProgram.
     * It calls ExactSCIP externally and reads back the optimal
     * solution from a file.
     * For this the constructor creates a temporary directory.
     *
     * \param name
     *   Name of the oracle.
     * \param exactBinary
     *   Binary of ExactSCIP.
     * \param mip
     *   Mixed-integer program to optimize over.
     * \param heuristic
     *   If not \c NULL, this other oracle is called if optimality is not required.
     * \param timeLimit
     *   Time limit (in seconds) for the solver.
     */

    ExactSCIPOptimizationOracle(const std::string& name, const std::string& exactBinary, MixedIntegerProgram& mip,
        FaceOptimizationOracleBase* heuristic, double timeLimit = std::numeric_limits<double>::max());

    /**
     * \brief Destructor.
     *
     * Destructor.
     * It deletes the temporary directory.
     */

    virtual ~ExactSCIPOptimizationOracle();

    /**
     * Whether to scale the objective before passing it to ExactSCIP.
     */

    bool scaleObjective;

  protected:

    /**
     * \brief Actual implementation.
     *
     * The actual implementation.
     * It writes the MIP with the current objective vector to
     * a ZIMPL model and calls the ExactSCIP solver as an external program,
     * redirecting the output into the temporary directory.
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    /**
     * \brief Method that is called when a new face is activated.
     *
     * Method that is called when a new \c face is activated.
     * If a (non-trivial) face was active before,
     * it is ensured that \ref faceDisabled() is called before.
     *
     * The implementation adds an equation constraint
     * to the underlying mixed-integer program.
     */

    virtual void faceEnabled(Face* face);

    /**
     * \brief Method that is called when a face is deactivated.
     *
     * Method that is called when a face is deactivated.
     *
     * The implementation removes the constraint that was added in \ref faceEnabled().
     */

    virtual void faceDisabled(Face* face);

    /**
     * \brief Creates a temporary directory.
     *
     * Creates a temporary directory that contains
     * the model file and logs.
     */

    void createTempDirectory();

    /**
     * \brief Removes the temporary directory.
     *
     * Removes the temporary directory created in \ref createTempDirectory().
     */

    void deleteTempDirectory();

    /**
     * \brief Writes a ZIMPL model into the temporary directory.
     *
     * Writes a ZIMPL model into the temporary directory.
     */

    void writeModel(const soplex::VectorRational& objective);

    /**
     * Calls the ExactSCIP solver.
     */

    void callSolver();

    /**
     * \brief Parses ExactSCIP's output.
     *
     * Parses ExactSCIP's output and extracts the solution vector.
     */

    soplex::DSVectorRational* parseOutput();

  protected:

    const std::string _binary; // The ExactSCIP binary.
    MixedIntegerProgram& _mip; // Reference to the underlying MIP.
    FaceOptimizationOracleBase* _heuristic; // Optional heuristic oracle.
    std::string _path; // Path to the temporary directory.
    double _timeLimit; // Time limit for the solver (in seconds).
  };

} /* namespace ipo */

#endif /* IPO_SCIP_ORACLES_H_ */
