#include "wrapper_oracle.h"

#include <unistd.h>
#include <sys/stat.h>
#include <string>
#include "cpu_timer.h"

using namespace soplex;

namespace ipo {

  WrapperOptimizationOracle::WrapperOptimizationOracle(const std::string& name, const std::string& wrapper,
      const std::string& instance) :
      OptimizationOracleBase(name), _wrapper(wrapper), _instance(instance)
  {
    /// Create temporary directory.

    char buffer[256] = "/tmp/ipo-WrapperOracle-XXXXXX";
    char* dirName = mkdtemp(buffer);
    if (dirName == NULL)
      throw std::runtime_error("Cannot create temporary directory!");
    _path = std::string(dirName);
    /// Initialize

    std::stringstream initStream;
    call("--init", initStream);
    std::vector<std::string> varNames;
    std::size_t n = 0;
    std::string info;
    initStream >> info;
    if (info == "variables")
    {
      initStream >> n;
      varNames.resize(n);
      for (std::size_t v = 0; v < n; ++v)
      {
        initStream >> varNames[v];
      }
    }
    initialize(varNames);
  }

  WrapperOptimizationOracle::~WrapperOptimizationOracle()
  {
    if (_path.find("/tmp/") != 0)
      throw std::runtime_error(
          "Error in ~WrapperOptimizationOracle: Temporary directory \"" + _path + "\" is not in /tmp/.");
    int status = system(("rm -r " + _path).c_str());
  }

  void WrapperOptimizationOracle::run(OptimizationResult& result, const VectorRational& objective,
      const Rational* improveValue, bool forceOptimal)
  {
    result.reset(numVariables());
    result.optimal = true;

    DVectorRational scaledObjective;
    scaleVectorIntegral(objective, scaledObjective);

    std::stringstream param, output;
    param << (forceOptimal ? "--optimize \"" : "--heuristic \"");
    for (std::size_t v = 0; v < numVariables(); ++v)
    {
      if (v > 0)
        param << ' ';
      param << scaledObjective[v];
    }
    param << "\"";

    call(param.str(), output);

    std::string status;
    output >> status;
    if (status == "infeasible")
      result.setInfeasible();
    else if (status == "unbounded")
    {
      result.rays.push_back(parseSolution(output));
      result.setUnbounded();
    }
    else if (status == "optimal")
    {
      result.points.push_back(parseSolution(output));
      result.setFeasible(objective);
    }
    else if (status == "feasible")
    {
      result.points.push_back(parseSolution(output));
      result.setFeasible(objective);
      result.optimal = false;
    }
    else
    {
      throw std::runtime_error("WrapperOracle: Invalid status \"" + status + "\".");
    }
  }

  DSVectorRational* WrapperOptimizationOracle::parseSolution(std::stringstream& stream)
  {
    DSVectorRational* solution = new DSVectorRational;
    for (std::size_t v = 0; v < numVariables(); ++v)
    {
      Rational x;
      std::string value;
      stream >> value;
      if (!x.readString(value.c_str()))
      {
        delete solution;
        throw std::runtime_error("WrapperOracle: Error while parsing solution vector.");
      }
      if (x != 0)
        solution->add(v, x);
    }
    return solution;
  }

  void WrapperOptimizationOracle::call(const std::string& parameters, std::stringstream& output)
  {
    const std::string cmd = "/usr/bin/time -o '" + _path + "/time.log' -f '%U' " + _wrapper + " " + parameters + " "
        + _path + " " + _instance;
//    std::cerr << "`" << cmd << "`" << std::endl;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe)
      throw std::runtime_error("WrapperOptimizationOracle: Failed to create pipe for wrapper.");
    char buffer[128];
    while (!feof(pipe))
    {
      if (fgets(buffer, 128, pipe) != NULL)
        output << buffer;
    }
    pclose(pipe);

    std::ifstream timing((_path + "/time.log").c_str());
    double time;
    timing >> time;
    addTimeToActiveTimers(time);
  }

} /* namespace ipo */
