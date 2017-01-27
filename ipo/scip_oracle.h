#ifndef IPO_SCIP_ORACLE_H_
#define IPO_SCIP_ORACLE_H_

#include "common.h"

#ifdef IPO_WITH_SCIP

#include <set>
#include <map>
#include <limits>

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
     * \brief Constructs a SCIP oracle from a file specified by \p fileName, optionally associated to \p nextOracle.
     *
     * Constructs a SCIP oracle from a file specified by \p fileName, optionally associated to \p nextOracle. The ambient space is
     * defined by the model read from the file (and must be equal to that of \p nextOracle). The oracle is implemented by calling 
     * SCIP in order to solve mixed-integer programs.
     */

    SCIPOracle(const std::string& fileName, const std::shared_ptr<OracleBase>& nextOracle = NULL);

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

    virtual void setFace(const LinearConstraint& newFace = completeFaceConstraint());

    /**
     * \brief Sets a time limit for each oracle call.
     *
     * Sets a time limit (in seconds) for each oracle call. If \c heuristicLevel is 0, this raises an exception.
     */

    double setTimeLimit(double timeLimit);

    /**
     * \brief Returns the time limit set for each oracle call.
     *
     * Returns the time limit set via setTimeLimit().
     */

    double getTimeLimit();

  protected:

    std::shared_ptr<MixedIntegerSet> constructFromSCIP(SCIP* originalSCIP);

    void constructFromMixedIntegerSet(const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet);

    std::shared_ptr<MixedIntegerSet> constructFromFile(const std::string& fileName);

    virtual void solverMaximize(double* objective, double objectiveBound, std::vector<double*>& points,
      std::vector<double*>& rays, bool& hitLimit);

  protected:
    SCIP* _scip; // SCIP instance
    std::vector<SCIP_VAR*> _variables; // SCIP variables.
    SCIP_CONS* _faceConstraint; // Special equation constraint for optimizing over a face.
  };

} /* namespace ipo */

#endif /* IPO_WITH_SCIP */

#endif /* IPO_SCIP_ORACLE_H_ */
