#include "external_oracle.h"

#include <unistd.h>
#include <sys/stat.h>
#include <string>

#include "cpu_timer.h"

using namespace soplex;

namespace ipo {

  ExternalOracle::ExternalOracle(const std::string& name, const std::string& program, const std::string& instance, 
    const std::shared_ptr<OracleBase>& nextOracle, std::size_t maxInfeasibleIterations, double initialM)
    : FaceOracleBase(name, nextOracle, maxInfeasibleIterations, initialM)
  {
    Space externalSpace;
    initialize(externalSpace);
//     if (externalSpace != _space)
//       throw std::runtime_error("Spaces differ while constructing ExternalOracle.");

    FaceOracleBase::initializeSpace(externalSpace);
  }

  void ExternalOracle::initialize(Space& externalSpace)
  {
    /// Create temporary directory.

    char buffer[256] = "/tmp/ipo-external-oracle-XXXXXX";
    char* dirName = mkdtemp(buffer);
    if (dirName == NULL)
      throw std::runtime_error("Cannot create temporary directory!");
    _path = std::string(dirName);

    std::stringstream initStream;
    call("--init", initStream);
    std::vector<std::string> varNames;
    std::size_t n = 0;
    std::string info;
    initStream >> info;
    if (info == "variables")
    {
      initStream >> n;
      for (std::size_t v = 0; v < n; ++v)
      {
        std::string varName;
        initStream >> varName;
//         externalSpace.addVariable(varName); // TODO: building of space after oracle construction is not allowed yet!
      }
    }
  }

  ExternalOracle::~ExternalOracle()
  {
    if (_path.find("/tmp/") != 0)
      throw std::runtime_error(
          "Error in ~ExternalOracle: Temporary directory \"" + _path + "\" is not in /tmp/.");
    int status = system(("rm -r " + _path).c_str());
  }
  
  std::size_t ExternalOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    DVectorRational scaledObjective;
    scaleIntegral(objective, scaledObjective);

    std::stringstream param, output;
    param << "--maximize \"";
    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      if (v > 0)
        param << ' ';
      param << scaledObjective[v];
    }
    param << "\"";

    call(param.str(), output);

    std::string status;
    output >> status;

    if (status == "unbounded")
    {
      Vector solution = parseSolution(output);
      result.rays.push_back(OracleResult::Ray(solution));  
    }
    else if (status == "optimal")
    {
      Vector solution = parseSolution(output);
      result.points.push_back(OracleResult::Point(solution));
      sort = true;
      checkDups = true;
    }
    else if (status != "infeasible")
      throw std::runtime_error("ExternalOracle: Invalid status \"" + status + "\".");
    
    return heuristicLevel();
  }

//   void ExternalOracle::unrestrictedMaximize(OracleResult& result, const VectorRational& objective,
//     const ObjectiveBound& improveValue, const VectorRational& originalObjective,
//     const ObjectiveBound& orginalObjectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic)
//   {
//     assert((heuristicLevel() == 0 && _nextOracle == NULL)
//       || heuristicLevel() > 0 && _nextOracle != NULL);
// 
//     // Forward call if requested.
// 
//     if (heuristicLevel() > maxHeuristic)
//     {
//       return _nextOracle->maximize(result, originalObjective, orginalObjectiveBound, maxHeuristic,
//         minHeuristic);
//     }
// 
//     DVectorRational scaledObjective;
//     scaleVectorIntegral(objective, scaledObjective);
// 
//     std::stringstream param, output;
//     param << "--maximize \"";
//     for (std::size_t v = 0; v < space().dimension(); ++v)
//     {
//       if (v > 0)
//         param << ' ';
//       param << scaledObjective[v];
//     }
//     param << "\"";
// 
//     call(param.str(), output);
// 
//     std::string status;
//     output >> status;
//     result.buildStart(objective);
//     if (status == "infeasible")
//     {
//       return result.buildFinish(heuristicLevel(), false, false, false);
//     }
//     else if (status == "unbounded")
//     {
//       result.buildAddDirection(parseSolution(output));
//       return result.buildFinish(heuristicLevel(), false, false, false);
//     }
//     else if (status == "optimal")
//     {
//       result.buildAddPoint(parseSolution(output));
//       return result.buildFinish(heuristicLevel(), true, true, true);
//     }
//     else
//     {
//       throw std::runtime_error("ExternalOracle: Invalid status \"" + status + "\".");
//     }
//   }

  Vector ExternalOracle::parseSolution(std::stringstream& stream)
  {
    VectorData* data = new VectorData(space().dimension());
    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      Rational x;
      std::string value;
      stream >> value;
      if (!x.readString(value.c_str()))
        throw std::runtime_error("ExternalOracle: Error while parsing solution vector.");
      if (x != 0)
        data->add(v, x);
    }
    return Vector(data);
  }

  void ExternalOracle::call(const std::string& parameters, std::stringstream& output)
  {
    const std::string cmd = "/usr/bin/time -o '" + _path + "/time.log' -f '%U' "
      + _program + " " + parameters + " " + _path + " " + _instance;
//    std::cerr << "`" << cmd << "`" << std::endl;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe)
      throw std::runtime_error("ExternalOracle: Failed to create pipe for wrapper.");
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
