#ifndef WRAPPER_ORACLE_H_
#define WRAPPER_ORACLE_H_

#include "oracles.h"

namespace ipo {

  /**
   *
   * Calls `wrapper` with the following arguments, followed by additional parameters which
   * define the instance. All lists are delimited by whitespace or newline.
   *
   * --init:          Prints, separated by whitespace, "variables", the number of variables
   *                  and a list of all variables names.
   *
   * --optimize OBJ:  Maximizes objective OBJ, which is a single argument representing a
   *                  dense integral vector (as a list). It must print either
   *                  'infeasible',
   *                  'unbounded' followed by a list of coordinates of an unbounded direction, or
   *                  'optimal' followed by a list of coordinates of an optimal solution.
   *
   * --heuristic OBJ: Maximizes objective OBJ, which is a single argument representing a
   *                  dense integral vector (as a list), heuristically. It must print either
   *                  'infeasible',
   *                  'unbounded' followed by a list of coordinates of an unbounded direction,
   *                  'optimal' followed by a list of coordinates of an optimal solution, or
   *                  'feasible' followed by a list of coordinates of a feasible solution.
   */

  class WrapperOptimizationOracle: public OptimizationOracleBase
  {
  public:
    WrapperOptimizationOracle(const std::string& name, const std::string& wrapper, const std::string& instance);
    virtual ~WrapperOptimizationOracle();

  protected:
    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

  private:
    soplex::DSVectorRational* parseSolution(std::stringstream& stream);
    void call(const std::string& parameters, std::stringstream& output);

  protected:
    std::string _wrapper;
    std::string _instance;
    std::string _path;
  };

} /* namespace ipo */

#endif /* WRAPPER_ORACLE_H_ */
