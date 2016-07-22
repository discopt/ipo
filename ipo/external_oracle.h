#ifndef EXTERNAL_ORACLE_H_
#define EXTERNAL_ORACLE_H_

#include "common.h"
#include "oracles.h"

namespace ipo {

  /**
   * \brief An oracle that uses an external program.
   *
   * An oracle that uses an external program. It calls a wrapper with the following arguments,
   * followed by additional parameters which define the instance. All lists are delimited by
   * whitespace or newline.
   *
   * \li \c --init:
   *   Prints, separated by whitespace, "variables", the number of variables
   *   and a list of all variables names.
   * \li \c --maximize OBJ:
   *   Maximizes objective OBJ, which is a single argument representing a
   *   dense integral vector (as a list). It must print either
   *   'infeasible',
   *   'unbounded' followed by a list of coordinates of an unbounded direction, or
   *   'optimal' followed by a list of coordinates of an optimal solution.
   */

  class ExternalOracle: public FaceOracleBase
  {
  public:
    /**
     * \brief Creates an instance.
     *
     * Creates an instance together with a temporary working directory.
     *
     * \param name               Name of the oracle.
     * \param program            Location of the external program.
     * \param instance           Specification of the instance.
     * \param space              Ambient space.
     * \param numBlindIterations Number of times, \f$ M \f$ is increased, before testing whether
     *                           \f$ F = \emptyset \f$ holds.
     * \param initialM           Initial value of \f$ M \f$.
     */

    ExternalOracle(const std::string& name, const std::string& program,
      const std::string& instance, const Space& space, std::size_t numBlindIterations = 2,
      double initialM = 16);

    /**
     * \brief Creates an instance.
     *
     * Creates an instance together with a temporary working directory.
     *
     * \param name               Name of the oracle.
     * \param program            Location of the external program.
     * \param instance           Specification of the instance.
     * \param space              Reference to ambient space that is defined by the oracle.
     * \param numBlindIterations Number of times, \f$ M \f$ is increased, before testing whether
     *                           \f$ F = \emptyset \f$ holds.
     * \param initialM           Initial value of \f$ M \f$.
     */

    ExternalOracle(const std::string& name, const std::string& program,
      const std::string& instance, Space& space, std::size_t numBlindIterations = 2,
      double initialM = 16);

    /**
     * \brief Creates an instance associated to \p nextOracle.
     *
     * Creates an instance associated to \p nextOracle together with a temporary working directory.
     *
     * \param name               Name of the oracle.
     * \param program            Location of the external program.
     * \param instance           Specification of the instance.
     * \param nextOracle         Next oracle to forward calls to.
     * \param numBlindIterations Number of times, \f$ M \f$ is increased, before testing whether
     *                           \f$ F = \emptyset \f$ holds.
     * \param initialM           Initial value of \f$ M \f$.

     */

    ExternalOracle(const std::string& name, const std::string& program,
      const std::string& instance, OracleBase* nextOracle,
      std::size_t numBlindIterations = 2, double initialM = 16);

    /**
     * \brief Destructor.
     *
     * Destructor. It removes the temporary directory created in the constructor.
     */

    virtual ~ExternalOracle();

  protected:

    /**
     * \brief Implementation of the oracle.
     *
     * Implements the optimization oracle to maximize the given dense rational \p objective over
     * the whole polyhedron \f$ P \f$ (in contrast to the maximize() methods) and returns \p result.
     * If the objective value requested by \p improveValue is not exceeded, then the call must be
     * forwarded to the next oracle using \p originalObjective and \p originalImproveValue to allow
     * the next oracle to use (potentially more efficient) handling of face optimization.
     *
     * \param result               After the call, contains the oracle's answer.
     * \param objective            Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized over
     *                             \f$ P \f$.
     * \param objectiveBound       Objective value \f$ \gamma \in \mathbb{Q} \f$ that should be
     *                             exceeded.
     * \param originalObjective    Objective vector \f$ c' \in \mathbb{Q}^n \f$ to be maximized
     *                             over \f$ F \f$.
     * \param originalImproveValue Objective value \f$ \gamma' \in \mathbb{Q} \f$ that should be
     *                             exceeded when optimizing \f$ c' \f$ over \f$ F \f$.
     * \param maxHeuristic         Requested maximum heuristic level.
     * \param minHeuristic         Requested minimum heuristic level.
     *
     * This implementation calls the external program and reads back its output.
     */

    virtual void unrestrictedMaximize(OracleResult& result,
      const soplex::VectorRational& objective, const ObjectiveBound& improveValue,
      const soplex::VectorRational& originalObjective, const ObjectiveBound& orginalObjectiveBound,
      std::size_t maxHeuristic, std::size_t minHeuristic);

  private:

    void initialize(Space& externalSpace);

    /**
     * \brief Parses the output of the external program.
     *
     * Parses the output of the external program.
     */

    soplex::DSVectorRational* parseSolution(std::stringstream& stream);

    /**
     * \brief Calls the external program.
     *
     * Calls the external program with arguments
     * as specified in the class description.
     * The current directory is the temporary one created in the constructor.
     */

    void call(const std::string& parameters, std::stringstream& output);

  protected:
    std::string _program; // Location of the external program.
    std::string _instance; // Specification of the instance.
    std::string _path; // Temporary directory.
  };

} /* namespace ipo */

#endif /* EXTERNAL_ORACLE_H_ */
