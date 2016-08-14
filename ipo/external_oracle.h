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
   *   'unbounded' followed by a list of coordinates of an unbounded rays, or
   *   'optimal' followed by a list of coordinates of an optimal solution.
   */

  class ExternalOracle: public FaceOracleBase
  {
  public:
    /**
     * \brief Creates an instance associated to \p nextOracle.
     *
     * Creates an instance associated to \p nextOracle together with a temporary working directory.
     *
     * \param name                    Name of the oracle.
     * \param program                 Location of the external program.
     * \param instance                Specification of the instance.
     * \param nextOracle              Next oracle to forward calls to.
     * \param maxInfeasibleIterations Maximum number of iterations before (heuristically) checking if the face is empty.
     * \param initialM                Initial value of \f$ M \f$.
     */

    ExternalOracle(const std::string& name, const std::string& program, const std::string& instance,
      const std::shared_ptr<OracleBase>& nextOracle = NULL, std::size_t maxInfeasibleIterations = 4, double initialM = 16);

    /**
     * \brief Destructor.
     *
     * Destructor. It removes the temporary directory created in the constructor.
     */

    virtual ~ExternalOracle();

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
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual std::size_t maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);

  private:

    void initialize(Space& externalSpace);

    /**
     * \brief Parses the output of the external program.
     *
     * Parses the output of the external program.
     */

    Vector parseSolution(std::stringstream& stream);

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
