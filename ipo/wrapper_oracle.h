#ifndef WRAPPER_ORACLE_H_
#define WRAPPER_ORACLE_H_

#include "ipo.h"
#include "oracles.h"

namespace ipo {

  /**
   * \brief An oracle that uses an external program.
   *
   * An oracle that uses an external program.
   * It calls a wrapper with the following arguments, followed by additional parameters which
   * define the instance. All lists are delimited by whitespace or newline.
   *
   * \li \c --init:
   *   Prints, separated by whitespace, "variables", the number of variables
   *   and a list of all variables names.
   * \li \c --optimize OBJ:
   *   Maximizes objective OBJ, which is a single argument representing a
   *   dense integral vector (as a list). It must print either
   *   'infeasible',
   *   'unbounded' followed by a list of coordinates of an unbounded direction, or
   *   'optimal' followed by a list of coordinates of an optimal solution.
   * \li \c --heuristic OBJ:
   *   Maximizes objective OBJ, which is a single argument representing a
   *   dense integral vector (as a list), heuristically. It must print either
   *   'infeasible',
   *   'unbounded' followed by a list of coordinates of an unbounded direction,
   *   'optimal' followed by a list of coordinates of an optimal solution, or
   *   'feasible' followed by a list of coordinates of a feasible solution.
   */

  class WrapperOptimizationOracle: public OptimizationOracleBase
  {
  public:
    /**
     * \brief Constructor.
     *
     * Constructor.
     * It creates a temporary directory to work in.
     *
     * \param name
     *   Name of the oracle.
     * \param wrapper
     *   Location of the wrapper program.
     * \param instance
     *   Specification of the instance.
     */

    WrapperOptimizationOracle(const std::string& name, const std::string& wrapper, const std::string& instance);

    /**
     * \brief Destructor.
     *
     * Destructor.
     * It removes the temporary directory created in the constructor.
     */

    virtual ~WrapperOptimizationOracle();

  protected:

    /**
     * \brief Implementation of the oracle.
     *
     * Implementation of the oracle.
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

  private:
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
    std::string _wrapper; // Location of the wrapper
    std::string _instance; // Specification of the instance
    std::string _path; // Temporary directory.
  };

} /* namespace ipo */

#endif /* WRAPPER_ORACLE_H_ */
