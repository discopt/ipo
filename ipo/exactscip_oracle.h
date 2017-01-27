#ifndef IPO_EXACTSCIP_ORACLE_H_
#define IPO_EXACTSCIP_ORACLE_H_

#include <set>
#include <map>
#include <limits>

#include "common.h"
#include "oracles.h"
#include "mip.h"

namespace ipo {

  /**
   * \brief An oracle based on the exact IP solver scip-ex.
   *
   * An oracle for the given \ref MixedIntegerSet that uses scip-ex (via an external call) to optimize.
   */

  class ExactSCIPOracle : public OracleBase
  {
  public:
#ifdef IPO_WITH_EXACT_SCIP
    /**
     * \brief Constructs an exact SCIP oracle with given \p name in given \p space.
     *
     * Constructs an exact SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is
     * equal to that of \p nextOracle and of the space of the given \p mixedIntegerSet. The oracle is implemented by calling
     * scip-ex (as detected by cmake) to solve mixed-integer programs over the \p mixedIntegerSet.
     *
     * \note This constructor throws an exception if cmake did not find scip-ex during the build.
     */

    ExactSCIPOracle(const std::string& name, const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet,
      const std::shared_ptr<OracleBase>& nextOracle = NULL);
#endif

    /**
     * \brief Constructs an exact SCIP oracle with given \p name in given \p space.
     *
     * Constructs an exact SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is
     * equal to that of \p nextOracle and of the space of the given \p mixedIntegerSet. The oracle is implemented by calling
     * the scip-ex \p binary to solve mixed-integer programs over the \p mixedIntegerSet.
     */

    ExactSCIPOracle(const std::string& binary, const std::string& name, const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet,
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

    virtual void setFace(const LinearConstraint& newFace = completeFaceConstraint());

    /**
     * \brief Sets a time limit for each oracle call.
     *
     * Sets a time limit (in seconds) for each oracle call. If \c heuristicLevel is 0, this raises an exception.
     */

    void setTimeLimit(double timeLimit);

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
    double _timeLimit;
  };

} /* namespace ipo */

#endif /* IPO_EXACTSCIP_ORACLE_H_ */
