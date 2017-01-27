#include "exactscip_oracle.h"

#include <limits>
#include <unistd.h>
#include <sys/stat.h>
#include "timer.h"

using namespace soplex;

namespace ipo {

#ifdef IPO_WITH_EXACT_SCIP

  ExactSCIPOracle::ExactSCIPOracle(const std::string& name, const std::shared_ptr< MixedIntegerSet >& mixedIntegerSet,
    const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase(name, nextOracle), _mixedIntegerSet(mixedIntegerSet), _timeLimit(0)
  {
    _binary = IPO_EXACT_SCIP_PATH;
    createWorkingDirectory();

    initializeSpace(_mixedIntegerSet->space());
  }

#endif

   ExactSCIPOracle::ExactSCIPOracle(const std::string& binary, const std::string& name,
     const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet, const std::shared_ptr<OracleBase>& nextOracle)
     : OracleBase(name, nextOracle), _binary(binary), _timeLimit(0)
   {
     createWorkingDirectory();

     initializeSpace(_mixedIntegerSet->space());
   }


  ExactSCIPOracle::~ExactSCIPOracle()
  {
    deleteWorkingDirectory();
  }

  void ExactSCIPOracle::setFace(const LinearConstraint& newFace)
  {
    OracleBase::setFace(newFace);
  }

  double ExactSCIPOracle::setTimeLimit(double timeLimit)
  {
    assert(timeLimit >= 0);
    if (timeLimit != _timeLimit)
    {
      _timeLimit = timeLimit;
      deleteWorkingDirectory();
      createWorkingDirectory();
    }
    return _timeLimit;
  }

  double ExactSCIPOracle::getTimeLimit()
  {
    return _timeLimit;
  }

  void ExactSCIPOracle::createWorkingDirectory()
  {
    if (_binary.empty())
      throw std::runtime_error("ExactSCIPOracle failed to initialze: Path of binary not specified.");

    char buffer[256] = "/tmp/ipo-scipex-XXXXXX";
    char* name = mkdtemp(buffer);
    if (name == NULL)
      throw std::runtime_error("Cannot create temporary directory!");
    _workingDirectory = std::string(name);

    std::string binary = "/usr/bin/time -o '" + _workingDirectory + "/timing.log' -f '%U' " + _binary + " ";
    std::string parameters = "";
    if (_timeLimit < std::numeric_limits<double>::max() && _timeLimit > 0)
    {
      std::stringstream ss;
      ss << " -c \"set limits time ";
      ss << _timeLimit;
      ss << "\" ";
      parameters += ss.str();
    }
    std::string failureParameters = " -c \"set misc usefprelax FALSE\" -c \"set presolving maxrounds 0\" ";
    std::string commands = "-c \"read model.zpl\" -c optimize -c \"display solution\" -c quit ";

    std::ofstream file((_workingDirectory + "/script.sh").c_str());
    file << "#!/bin/bash\n\n";
    file << "cd " << _workingDirectory << "\n";
    file << binary << " " << parameters << commands << " > solve.log 2>&1\n";
    file << "retcode=$?\n";
    file << "if [[ $retcode != 0 ]]; then\n";
    file << "  " << binary << parameters << failureParameters << commands << " > solve.log 2>&1\n";
    file << "fi\n";
    file.close();
    chmod((_workingDirectory + "/script.sh").c_str(), 00700);
  }

  void ExactSCIPOracle::deleteWorkingDirectory()
  {
    unlink((_workingDirectory + "/solve.log").c_str());
    unlink((_workingDirectory + "/script.sh").c_str());
    unlink((_workingDirectory + "/model.zpl").c_str());
    unlink((_workingDirectory + "/timing.log").c_str());
    rmdir(_workingDirectory.c_str());
  }

  HeuristicLevel ExactSCIPOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    if (_binary.empty())
    {
      throw std::runtime_error("ExactSCIPOracle failed to initialize: Path of binary not specified.");
    }

    sort = !result.points.empty();
    checkDups = !result.points.empty();

    writeModel(objective);
    solveModel();
    VectorData* pointData = parseOutput();
    if (pointData != NULL)
    {
      Vector point(pointData);
      result.points.push_back(OracleResult::Point(point));
      result.computeMissingObjectiveValues();
    }
    return heuristicLevel();
  }

  void ExactSCIPOracle::writeModel(const VectorRational& objective)
  {
    std::ofstream file((_workingDirectory + "/model.zpl").c_str());

    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      const MixedIntegerSet::Variable& var = _mixedIntegerSet->variable(v);
      file << "var x" << v;
      if (var.integral)
        file << " integer";
      else
        file << " real";
      if (var.lowerBound > -infinity)
        file << " >= " << var.lowerBound;
      if (var.upperBound < infinity)
        file << " <= " << var.upperBound;
      file << ";\n";
    }
    file << "\nmaximize cost:";
    bool first = true;
    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      if (objective[v] == 0)
        continue;
      file << "\n";
      if (first)
        first = false;
      else
        file << " + ";
      file << objective[v] << "*x" << v;
    }
    if (first)
      file << "0*x0";
    file << ";\n\n";
    for (std::size_t r = 0; r <= _mixedIntegerSet->numRows(); ++r)
    {
      const LinearConstraint& row = r < _mixedIntegerSet->numRows() ? _mixedIntegerSet->rowConstraint(r) : currentFace();
      file << "\nsubto row" << r << ":";
      first = true;
      for (std::size_t p = 0; p < row.normal().size(); ++p)
      {
        file << "\n";
        if (first)
          first = false;
        else
          file << " + ";
        file << row.normal().value(p) << "*x" << row.normal().index(p);
      }
      if (first)
        file << "x0 - x0";
      file << ' ' << row.type() << "= " << row.rhs() << ";\n\n";
    }
    file.close();
  }

  void ExactSCIPOracle::solveModel()
  {
    if (system((_workingDirectory + "/script.sh").c_str()) != 0)
    {
      throw std::runtime_error("ExactSCIPOracle returned with nonzero exit status (see " + _workingDirectory+ "/solve.log).");
    }
  }

  VectorData* ExactSCIPOracle::parseOutput()
  {
    std::ifstream log((_workingDirectory + "/solve.log").c_str());
    std::string line;
    bool startedSolutionSection = false;
    bool timeLimitReached = false;
    VectorData* pointData = NULL;
    while (std::getline(log, line))
    {
      if (line.substr(0, 16) == "objective value:")
      {
        startedSolutionSection = true;
        pointData = new VectorData();
      }

      if (startedSolutionSection && !line.empty() && line[0] == 'x')
      {
        std::size_t var;
        std::string valueStr;
        std::stringstream ss(line.substr(1, std::string::npos));
        ss >> var >> valueStr;
        Rational value;
        if (!value.readString(valueStr.c_str()))
          throw std::runtime_error("ExactSCIPOracle failed to parse output while reading a number (see " + _workingDirectory+ "/solve.log).");
        assert(pointData != NULL);
        pointData->add(var, value);
      }

      if (line == "no solution available")
      {
        return NULL;
      }
      if (line == "SCIP Status        : solving was interrupted [time limit reached]")
      {
        timeLimitReached = true;
        std::stringstream str;
        str << "Oracle \"" << name() << "\" reached its time limit. Log files remain here: " << _workingDirectory;
        throw std::runtime_error(str.str());
      }
    }

    std::ifstream timing((_workingDirectory + "/timing.log").c_str());
    double time;
    timing >> time;
    addTimeToRunningTimers(time);

    if (!startedSolutionSection)
    {
      throw std::runtime_error("ExactSCIPOracle did not return useful results (see " + _workingDirectory+ "/solve.log).");
    }

    return pointData;
  }

} /* namespace ipo */
