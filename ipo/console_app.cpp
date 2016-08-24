#include "console_app.h"

#include <cassert>
#include <random>
#include <vector>
#include <regex>
#include <cmath>

#include "affine_hull.h"
#include "min_norm_2d.h"
#include "facets.h"
#include "smallest_face.h"
#include "spx_gmp.h"
#include "parser.h"
#include "common.h"

using namespace soplex;

namespace ipo {

  class ConsoleApplicationObjectiveParser : public LPObjectiveParser
  {
  public:
    ConsoleApplicationObjectiveParser(std::istream& stream, const Space& space,
      std::vector<Vector>& objectives, std::vector<std::string>& objectiveNames)
      : LPObjectiveParser(stream), _objectives(objectives), _objectiveNames(objectiveNames)
    {
      for (std::size_t i = 0; i < space.dimension(); ++i)
        _oracleVariables[space[i]] = i;
    }

    virtual ~ConsoleApplicationObjectiveParser()
    {

    }

    virtual void handleObjective(const std::string& name, const std::map<std::string, Rational>& coefficients)
    {
      VectorData* vectorData = new VectorData(_oracleVariables.size());
      for (std::map<std::string, Rational>::const_iterator iter = coefficients.begin(); iter != coefficients.end(); ++iter)
      {
        std::map<std::string, std::size_t>::const_iterator varIter = _oracleVariables.find(iter->first);
        if (varIter != _oracleVariables.end())
        {
          if (iter->second != 0)
            vectorData->add(varIter->second, iter->second);
        }
        else
        {
          std::cerr << "Skipping objective: Unknown variable <" << iter->first << ">.\n" << std::endl;
          delete vectorData;
          return;
        }
      }
      _objectives.push_back(Vector(vectorData));
      _objectiveNames.push_back(name);
    }

  private:
    std::vector<Vector>& _objectives;
    std::vector<std::string>& _objectiveNames;
    std::map<std::string, std::size_t> _oracleVariables;
  };

  class ConsoleApplicationInequalityParser : public LPInequalityParser
  {
  public:
    ConsoleApplicationInequalityParser(std::istream& stream, const Space& space,
      std::vector<LinearConstraint>& faces, std::vector<std::string>& faceNames)
      : LPInequalityParser(stream), _faces(faces), _faceNames(faceNames), _space(space)
    {
      for (std::size_t i = 0; i < space.dimension(); ++i)
        _oracleVariables[space[i]] = i;
    }

    virtual ~ConsoleApplicationInequalityParser()
    {

    }

    virtual void handleInequality(const std::string& name, const Rational& lhs,
      const std::map<std::string, Rational >& coefficients, const Rational& rhs)
    {
      if ((lhs <= -infinity && rhs >= infinity) || lhs == rhs)
        return;

      VectorData* vectorData = new VectorData();
      for (std::map<std::string, Rational>::const_iterator iter = coefficients.begin(); iter != coefficients.end(); ++iter)
      {
        std::map<std::string, std::size_t>::const_iterator varIter = _oracleVariables.find(iter->first);
        if (varIter != _oracleVariables.end())
        {
          vectorData->add(varIter->second, iter->second);
        }
        else
        {
          std::cerr << "Skipping inequality: Unknown variable <" << iter->first << ">.\n" << std::endl;
          return;
        }
      }

      Vector vector(vectorData);
      bool hasBoth = lhs > -infinity && rhs < infinity;
      if (lhs > -infinity)
      {
        _faces.push_back(LinearConstraint(',', -vector, -lhs));
        _faceNames.push_back(hasBoth ? (name + "-lhs") : name);
      }
      if (rhs < infinity)
      {
        _faces.push_back(LinearConstraint('<', vector, rhs));
        _faceNames.push_back(hasBoth ? (name + "-rhs") : name);
      }
    }

  private:
    std::map<std::string, std::size_t> _oracleVariables;
    std::vector<LinearConstraint>& _faces;
    std::vector<std::string>& _faceNames;
    const Space& _space;
  };

  class ConsoleApplicationPointParser : public PointParser
  {
  public:
    ConsoleApplicationPointParser(std::istream& stream, const Space& space,
      std::vector<Vector>& points, std::vector<std::string>& pointNames)
      : PointParser(stream), _points(points), _pointNames(pointNames)
    {
      for (std::size_t i = 0; i < space.dimension(); ++i)
        _oracleVariables[space[i]] = i;
    }

    virtual ~ConsoleApplicationPointParser()
    {

    }

    virtual void handlePoint(const std::string& name, const std::map< std::string, Rational >& values)
    {
      VectorData* data = new VectorData(values.size());
      for (std::map<std::string, Rational>::const_iterator iter = values.begin(); iter != values.end(); ++iter)
      {
        std::map<std::string, std::size_t>::const_iterator varIter = _oracleVariables.find(iter->first);
        if (varIter != _oracleVariables.end())
        {
          data->add(varIter->second, iter->second);
        }
        else
        {
          std::cerr << "Skipping point: Unknown variable <" << iter->first << ">.\n" << std::endl;
          return;
        }
      }

      _points.push_back(Vector(data));
      _pointNames.push_back(name);
    }

  private:
    std::map<std::string, std::size_t> _oracleVariables;
    std::vector<Vector>& _points;
    std::vector<std::string>& _pointNames;
  };


  ConsoleApplicationBase::ConsoleApplicationBase(int numArguments, char** arguments)
  {
    _program = arguments[0];
    for (int i = 1; i < numArguments; ++i)
      _arguments.push_back(arguments[i]);

    _numRandomObjectives = 0;

    _taskPrintAmbientDimension = false;
    _taskPrintVariables = false;
    _taskMaximize = false;
    _taskMinimize = false;
    _taskDimension = false;
    _taskEquations = false;
    _taskSeparateFacet = false;
    _taskSmallestFace = false;
    _taskGenerateFacets = false;
    _taskPrintCached = false;
    _optionReadable = true;
    _optionCertificates = false;
    _optionReuseFacets = true;
    _optionPrintRandom = 0;
    _optionCache = true;

    _oracle = NULL;
    _cacheOracle = NULL;
    _projection = NULL;
    _projectionOracle = NULL;
  }

  ConsoleApplicationBase::~ConsoleApplicationBase()
  {
    if (_projection)
      delete _projection;
  }

  void ConsoleApplicationBase::setBasicOracle(std::shared_ptr<OracleBase> oracle)
  {
    if (_oracle)
      throw std::runtime_error("Error in ConsoleApplicationBase::setBasicOracle: Oracle already set.");
    if (_faceRestrictionArgument != "")
      throw std::runtime_error("Restricting to faces is not implemented, yet.");

    if (_projectionArgument != "")
    {
      const Space& originalSpace = oracle->space();
      std::regex re = std::regex(_projectionArgument);
      std::vector<std::size_t> projectionVariables;
      projectionVariables.reserve(originalSpace.dimension());
      for (std::size_t v = 0; v < originalSpace.dimension(); ++v)
      {
        if (std::regex_match(originalSpace[v], re))
          projectionVariables.push_back(v);
      }

      ProjectionData* projectionData = new ProjectionData(originalSpace, projectionVariables);
      _projection = new Projection(projectionData);
      _projectionOracle =  std::make_shared<ProjectionOracle>(*_projection, oracle);
      oracle = _projectionOracle;
    }

    if (_optionCache)
    {
      _cacheOracle = std::make_shared<CacheOracle>(oracle);
      _oracle = _cacheOracle;
    }
    else
    {
      _oracle = oracle;
    }
  }

  bool ConsoleApplicationBase::parseArguments()
  {
    try
    {
      while (numArguments() > 0)
      {
        std::size_t eat = parseArgument(_arguments.front());
        if (eat == 0)
        {
          std::stringstream ss;
          ss << "Argument \"" << _arguments.front() << "\", followed by " << (_arguments.size() - 1)
              << " further arguments could not be parsed.";
          throw std::runtime_error(ss.str());
        }
        std::min(eat, numArguments());
        for (std::size_t i = 0; i < eat; ++i)
          _arguments.pop_front();
      }
      return true;
    }
    catch (std::runtime_error& exc)
    {
      std::cerr << exc.what() << "\n\nCall again with --help to print usage.\n" << std::flush;
      return false;
    }
  }

  bool ConsoleApplicationBase::runDefaultTasks()
  {
    if (taskPrintAmbientDimension())
    {
      if (!printAmbientDimension())
        return false;
    }

    if (taskPrintVariables())
    {
      if (!printVariables())
        return false;
    }

    if (taskDimension() || taskEquations() || taskFacet() || taskFacets() || taskSmallestFace())
    {
      if (!computeAffineHull(completeFace(), ""))
        return false;
    }

    if (taskDimension() || taskEquations())
    {
      for (std::size_t i = 0; i < numFaces(); ++i)
      {
        std::string name;
        if (_faceNames[i] == "" || _faceNames[i] == "-lhs" || _faceNames[i] == "-rhs")
        {
          std::stringstream stream;
          stream << "#" << (i+1) << _faceNames[i];
          name = stream.str();
        }
        else
          name = _faceNames[i];
        if (!computeAffineHull(_faces[i], name))
          return false;
      }
    }

    if (taskMaximize() || taskMinimize() || taskFacets())
    {
      for (std::size_t i = 0; i < numObjectives(); ++i)
      {
        if (_objectiveNames[i] == "")
          std::cout << "Objective";
        else if (_objectiveNames[i] == ":RANDOM:")
        {
          if (_optionPrintRandom > 0)
            std::cout << "Random objective";
        }
        else
          std::cout << "Objective <" << _objectiveNames[i] << ">";
        int print = (_objectiveNames[i] != ":RANDOM:" ) ? 2 : _optionPrintRandom;
        if (print >= 2)
        {
          std::cout << " ";
          oracle()->space().printLinearForm(std::cout, objective(i));
        }
        if (print >= 1)
          std::cout << ":\n" << std::flush;

        if (taskMaximize())
        {
          if (!optimizeObjective(objective(i), true))
            return false;
        }
        if (taskMinimize())
        {
          if (!optimizeObjective(objective(i), false))
            return false;
        }
        if (taskFacets())
        {
          if (!generateFacets(objective(i), print))
            return false;
        }
      }
    }

    if (taskFacet())
    {
      for (std::size_t i = 0; i < numRays(); ++i)
      {
        std::cout << "Point";
        if (_rayNames[i] != "")
          std::cout << " " << _rayNames[i];
        std::cout << std::flush;

        bool isFeasible = true;
        if (!separateRayFacet(ray(i), isFeasible))
          return false;
      }
    }

    if (taskFacet() || taskSmallestFace())
    {
      for (std::size_t i = 0; i < numPoints(); ++i)
      {
        std::cout << "Point";
        if (_pointNames[i] != "")
          std::cout << " " << _pointNames[i];
        std::cout << std::flush;

        bool isFeasible = true;
        if (taskFacet())
        {
          if (!separatePointFacet(point(i), isFeasible))
            return false;
        }
        else
        {
          std::cout << " is assumed to be feasible.\n" << std::flush;
        }
        if (taskSmallestFace() && isFeasible)
        {
          if (!computeSmallestFace(point(i)))
            return false;
        }
      }
    }

    if (taskPrintCached())
    {
      if (!printCached())
        return false;
    }

    return true;
  }

  bool ConsoleApplicationBase::run()
  {
    if (!parseArguments())
      return false;
    if (!processArguments())
      return false;
    if (!runDefaultTasks())
      return false;

    return true;
  }

  void ConsoleApplicationBase::printUsage()
  {
    std::cerr << "Usage: " << _program << " [OPTIONS] INSTANCE...\n";
    std::cerr << "\n";
    std::cerr << "Options for defining the polyhedron P:\n";
    printAdditionalOptionsPolyhedron(std::cerr);
    std::cerr << "  --projection VARS\n";
    std::cerr << "      Consider the orthogonal projection onto a subset of variables. VARS is a regular expression\n";
    std::cerr << "      that defines the projection onto all variables matching it. Note that this way you can specify\n";
    std::cerr << "      a list of variables using '|' as a separator.\n";
    std::cerr << "  --restrict-face INEQ\n";
    std::cerr << "      Restrict the oracle to the face specified by INEQ. See --face below.\n";
    std::cerr << "\n";
    std::cerr << "Options for tasks (execution order is independent of argument order):\n";
    std::cerr << "  --print-ambient-dimension\n";
    std::cerr << "      Print the ambient dimension.\n";
    std::cerr << "  --print-variables\n";
    std::cerr << "      Print the names of all variables in a comma-separated list.\n";
    std::cerr << "  --maximize\n";
    std::cerr << "      Maximize objectives specified by --objective(s) arguments.\n";
    std::cerr << "  --minimize\n";
    std::cerr << "      Minimize objectives specified by --objective(s) arguments.\n";
    std::cerr << "  --dimension\n";
    std::cerr << "      Compute the dimension.\n";
    std::cerr << "  --equations\n";
    std::cerr << "      Compute valid equations.\n";
    std::cerr << "  --facet\n";
    std::cerr << "      Attempt to separate every given point or ray using a facet-defining inequality.\n";
    std::cerr << "  --smallest-face\n";
    std::cerr << "      Compute, for every given point in P, an objective vector that induces the smallest face\n";
    std::cerr << "      containing the point, and the dimension of that face.\n";
    std::cerr << "  --facets\n";
    std::cerr << "      Generate facets that are useful when maximizing the objectives specified by --objective(s) or\n";
    std::cerr << "      generated by --random.\n";
    std::cerr << "  --print-cached\n";
    std::cerr << "      Print the cached points and rays at the end.\n";
    printAdditionalOptionsTasks(std::cerr);
    std::cerr << "\n";
    std::cerr << "Options for input data (ordering w.r.t. tasks does not matter):\n";
    std::cerr << "  --face INEQ\n";
    std::cerr << "      Compute dimension/equations also for face induced by INEQ. INEQ is an inequality, potentially\n";
    std::cerr << "      preceeded by a name with a colon, e.g., \"myface: x#1 + 3y#5 <= 7\". Usable multiple times.\n";
    std::cerr << "  --faces FILE\n";
    std::cerr << "      The content of FILE is considered as a --face argument. If they do not match the pattern,\n";
    std::cerr << "      tokens are silently ignored, allowing FILE to be in LP format.\n";
    std::cerr << "  --objective OBJ\n";
    std::cerr << "      Use the given objective vector for optimization/facet generation. OBJ is an LP file objective,\n";
    std::cerr << "      e.g., \"max x#1 + 3y#5\". Usable multiple times.\n";
    std::cerr << "  --objectives FILE\n";
    std::cerr << "      The content of FILE is considered as a --objective argument. If they do not match the pattern,\n";
    std::cerr << "      tokens are silently ignored, allowing FILE to be in LP format.\n";
    std::cerr << "  --random N\n";
    std::cerr << "      Sample N objective vectors (uniformly at random from the sphere) to be used for facet\n";
    std::cerr << "      generation.\n";
    std::cerr << "  --point POINT\n";
    std::cerr << "      Use the point for smallest-face computation and facet-separation. POINT is a sparse vector,\n";
    std::cerr << "      e.g., \"(x#1, 3y#5)\". Usable multiple times.\n";
    std::cerr << "  --points FILE\n";
    std::cerr << "      The content of FILE is considered as a --point argument. If they do not match the pattern,\n";
    std::cerr << "      tokens are silently ignored.\n";
    std::cerr << "  --ray RAY\n";
    std::cerr << "      Use the given unbounded ray for facet-separation. RAY is a sparse vector, e.g.,\n";
    std::cerr << "      \"(x#1, 3y#5)\". Usable multiple times.\n";
    std::cerr << "  --rays FILE\n";
    std::cerr << "      The content of FILE is considered as a --ray argument. If they do not match the pattern,\n";
    std::cerr << "      tokens are silently ignored.\n";
    printAdditionalOptionsInput(std::cerr);
    std::cerr << "\n";
    std::cerr << "Further options (ordering does not matter):\n";
    std::cerr << "  --cache on|off\n";
    std::cerr << "      If on (default), cache points and rays produced by IPO to speed up further computations.\n";
    std::cerr << "  --readable on|off\n";
    std::cerr << "      If on (default), improve readability of inequalities/equations produced by IPO.\n";
    std::cerr << "  --certificates on|off\n";
    std::cerr << "      If on, output, for every computed facet, points and rays that span it. Off by default.\n";
    std::cerr << "  --trust-points on|off\n";
    std::cerr << "      If off (default), it is checked that points are in P when computing the smallest containing\n";
    std::cerr << "      faces.\n";
    std::cerr << "  --reuse-facets on|off\n";
    std::cerr << "      If on (default), use computed facets for subsequent facet-generation.\n";
    std::cerr << "  --print-random long|short|none\n";
    std::cerr << "      If none (default), do not print anything, if short, just print \"Random objective\", and if long,\n";
    std::cerr << "      also print the andom objective vector1.\n";
    std::cerr << "  --progress stdout|stderr|off\n";
    std::cerr << "      If on, print progress output to stdout or stderr. Off by default.\n";
    std::cerr << "  --debug stdout|stderr|off\n";
    std::cerr << "      If on, print debug output to stdout or stderr. Off by default.\n";
    printAdditionalOptionsFurther(std::cerr);
    printAdditionalOptionsSpecific(std::cerr);
    std::cerr << std::flush;
  }

  void ConsoleApplicationBase::setRelaxationBounds(const soplex::VectorRational& lowerBounds, const soplex::VectorRational& 
upperBounds)
  {
    if (_projectionOracle)
    {
      // TODO: Some of the bounds may be used in case of projection!
    }
    else
    {
      assert(lowerBounds.dim() == _relaxationColumns.num());
      assert(upperBounds.dim() == _relaxationColumns.num());

      _relaxationColumns.lower_w() = lowerBounds;
      _relaxationColumns.upper_w() = upperBounds;
    }
  }

  void ConsoleApplicationBase::addRelaxationConstraints(const std::vector<LinearConstraint>& constraints)
  {
    if (_projectionOracle)
    {
      // TODO: Some of the rows may be used in case of projection!
    }
    else
    {
      std::copy(constraints.begin(), constraints.end(), std::back_inserter(_relaxationConstraints));
    }
  }

  void ConsoleApplicationBase::printAdditionalOptionsPolyhedron(std::ostream& stream)
  {

  }

  void ConsoleApplicationBase::printAdditionalOptionsTasks(std::ostream& stream)
  {

  }

  void ConsoleApplicationBase::printAdditionalOptionsInput(std::ostream& stream)
  {

  }

  void ConsoleApplicationBase::printAdditionalOptionsFurther(std::ostream& stream)
  {

  }

  void ConsoleApplicationBase::printAdditionalOptionsSpecific(std::ostream& stream)
  {

  }

  int ConsoleApplicationBase::parseArgument(const std::string& firstArgument)
  {
    if (firstArgument == "--help")
    {
      printUsage();
      exit(0);
    }
    if (firstArgument == "--projection")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error(
            "Invalid option: --projection must be followed by a list of variables or a regular expression.");
      }
      _projectionArgument = argument(1);
      return 2;
    }
    if (firstArgument == "--restrict-face")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --restrict-face must be followed by an inequality.");
      }
      _faceRestrictionArgument = argument(1);
      return 2;
    }
    if (firstArgument == "--print-ambient-dimension")
    {
      _taskPrintAmbientDimension = true;
      return 1;
    }
    if (firstArgument == "--print-variables")
    {
      _taskPrintVariables = true;
      return 1;
    }
    if (firstArgument == "--maximize")
    {
      _taskMaximize = true;
      return 1;
    }
    if (firstArgument == "--minimize")
    {
      _taskMinimize = true;
      return 1;
    }
    if (firstArgument == "--dimension")
    {
      _taskDimension = true;
      return 1;
    }
    if (firstArgument == "--equations")
    {
      _taskEquations = true;
      return 1;
    }
    if (firstArgument == "--facet")
    {
      _taskSeparateFacet = true;
      return 1;
    }
    if (firstArgument == "--facets")
    {
      _taskGenerateFacets = true;
      return 1;
    }
    if (firstArgument == "--smallest-face")
    {
      _taskSmallestFace = true;
      return 1;
    }
    if (firstArgument == "--print-cached")
    {
      _taskPrintCached = true;
      return 1;
    }
    if (firstArgument == "--face")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --face must be followed by an inequality.");
      }
      _faceArguments.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--faces")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --faces must be followed by a file name.");
      }
      _faceFiles.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--objective")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --objective must be followed by a vector.");
      }
      _objectiveArguments.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--objectives")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --objectives must be followed by a file name.");
      }
      _objectiveFiles.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--ray")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --ray must be followed by a vector.");
      }
      _rayArguments.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--rays")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --rays must be followed by a file name.");
      }
      _rayFiles.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--point")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --point must be followed by a vector.");
      }
      _pointArguments.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--points")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --points must be followed by a file name.");
      }
      _pointFiles.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--random")
    {
      if (numArguments() == 1)
        throw std::runtime_error("Invalid option: --random must be followed by a positive number.");
      std::stringstream ss(argument(1));
      int n;
      ss >> n;
      if (n <= 0 || !bool(ss) || ss.good())
        throw std::runtime_error("Invalid option: --random must be followed by a positive number.");
      _numRandomObjectives += n;
      return 2;
    }
    if (firstArgument == "--cache")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _optionCache = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _optionCache = false;
          return 2;
        }
      }
      throw std::runtime_error(
        "Invalid option: --cache must be followed by either \"on\" or \"off\".");
    }
    if (firstArgument == "--readable")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _optionReadable = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _optionReadable = false;
          return 2;
        }
      }
      throw std::runtime_error(
        "Invalid option: --readable must be followed by either \"on\" or \"off\".");
    }
    if (firstArgument == "--certificates")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _optionCertificates = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _optionCertificates = false;
          return 2;
        }
      }
      throw std::runtime_error("Invalid option: --certificates must be followed by either \"on\" or \"off\".");
    }
    if (firstArgument == "--reuse-facets")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _optionReuseFacets = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _optionReuseFacets = false;
          return 2;
        }
      }
      throw std::runtime_error("Invalid option: --reuse-facets must be followed by either \"on\" or \"off\".");
    }
    if (firstArgument == "--print-random")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "yes" || argument(1) == "long" || argument(1) == "full")
        {
          _optionPrintRandom = 2;
          return 2;
        }
        else if (argument(1) == "short" || argument(1) == "info")
        {
          _optionPrintRandom = 1;
          return 2;
        }
        else if (argument(1) == "no" || argument(1) == "none" || argument(1) == "quiet")
        {
          _optionPrintRandom = 0;
          return 2;
        }
      }
      throw std::runtime_error("Invalid option: --print-random must be followed by either \"on\" or \"off\".");
    }

    return 0;
  }

  bool ConsoleApplicationBase::processArguments()
  {
    if (_oracle == NULL)
    {
      std::cerr << "Invalid invocation: No instance specified.\n\n";
      std::cerr << "Call again with --help to print usage.\n" << std::flush;
      return false;
    }

    std::size_t n = oracle()->space().dimension();

    _relaxationColumns.clear();
    DSVectorRational zeroVector;
    for (std::size_t c = 0; c < n; ++c)
      _relaxationColumns.add(0, -infinity, zeroVector, infinity);

    /// Add specified objectives.

    for (std::size_t i = 0; i < _objectiveArguments.size(); ++i)
    {
      std::istringstream stream(_objectiveArguments[i]);
      ConsoleApplicationObjectiveParser parser(stream, oracle()->space(), _objectives, _objectiveNames);
      parser.run();
    }

    /// Parse objectives from specified files.

    for (std::size_t i = 0; i < _objectiveFiles.size(); ++i)
    {
      std::ifstream stream(_objectiveFiles[i]);
      ConsoleApplicationObjectiveParser parser(stream, oracle()->space(), _objectives, _objectiveNames);
      parser.run();
    }

    /// Add random objectives.

    if (_numRandomObjectives > 0)
    {
//       std::random_device randomDevice; // TODO: valgrind does not like random_device.
      std::default_random_engine generator(0);
      std::normal_distribution<double> distribution;
      DVectorReal randomVector(n);
      for (std::size_t i = 0; i < _numRandomObjectives; ++i)
      {
        double norm = 0;
        for (std::size_t c = 0; c < n; ++c)
        {
          double x = distribution(generator);
          randomVector[c] = x;
          norm += x*x;
        }
        if (norm > 0)
        {
          norm = std::sqrt(norm);
          VectorData* objectiveVectorData = new VectorData(n);
          for (std::size_t v = 0; v < n; ++v)
            objectiveVectorData->add(v, Rational(randomVector[v] / norm));
          _objectives.push_back(Vector(objectiveVectorData));
          _objectiveNames.push_back(":RANDOM:");
        }
      }
    }

    /// Parse specified faces.

    for (std::size_t i = 0; i < _faceArguments.size(); ++i)
    {
      std::istringstream stream(_faceArguments[i]);
      ConsoleApplicationInequalityParser parser(stream, oracle()->space(), _faces, _faceNames);
      parser.run();
    }

    /// Parse faces from specified files.

    for (std::size_t i = 0; i < _faceFiles.size(); ++i)
    {
      std::ifstream stream(_faceFiles[i]);
      ConsoleApplicationInequalityParser parser(stream, oracle()->space(), _faces, _faceNames);
      parser.run();
    }

    /// Add specified rays.

    for (std::size_t i = 0; i < _rayArguments.size(); ++i)
    {
      std::istringstream stream(_rayArguments[i]);
      ConsoleApplicationPointParser parser(stream, oracle()->space(), _rays, _rayNames);
      parser.run();
    }

    /// Parse rays from specified files.

    for (std::size_t i = 0; i < _rayFiles.size(); ++i)
    {
      std::ifstream stream(_rayFiles[i]);
      ConsoleApplicationPointParser parser(stream, oracle()->space(), _rays, _rayNames);
      parser.run();
    }

     /// Add specified points.

    for (std::size_t i = 0; i < _pointArguments.size(); ++i)
    {
      std::istringstream stream(_pointArguments[i]);
      ConsoleApplicationPointParser parser(stream, oracle()->space(), _points, _pointNames);
      parser.run();
    }

    /// Parse points from specified files.

    for (std::size_t i = 0; i < _pointFiles.size(); ++i)
    {
      std::ifstream stream(_pointFiles[i]);
      ConsoleApplicationPointParser parser(stream, oracle()->space(), _points, _pointNames);
      parser.run();
    }

    return true;
  }

  bool ConsoleApplicationBase::printAmbientDimension()
  {
    std::cout << "Ambient dimension:\n " << oracle()->space().dimension() << "\n" << std::flush;
    return true;
  }

  bool ConsoleApplicationBase::printVariables()
  {
    std::cout << "Variables:\n";
    for (std::size_t i = 0; i < oracle()->space().dimension(); ++i)
    {
      const std::string& name = oracle()->space()[i];
      std::cout << " " << name << (i + 1 < oracle()->space().dimension() ? "," : "\n");
    }
    return true;
  }

  bool ConsoleApplicationBase::computeAffineHull(const LinearConstraint& face, const std::string& faceName)
  {
    if (_taskDimension || _taskEquations)
    {
      std::cout << "Computing the affine hull";
      if (!face.definesCompleteFace())
      {
        std::cout << " of face " << faceName << ":\n" << std::flush;
      }
      else
      {
        std::cout << ":\n" << std::flush;
      }
    }

    std::shared_ptr<OracleBase> orac = oracle();
    std::size_t n = oracle()->space().dimension();

    AffineHull::QuietOutput hullOutput;
    AffineHull::Result hull;

    std::vector<LinearConstraint> faceEquations;
    std::vector<LinearConstraint>* equations;
    if (face.definesCompleteFace())
    {
      hull.run(_equations, orac, hullOutput, 1, true);
      
      _basicColumns = hull.basicColumns();
      _spanningPoints.reserve(hull.numSpanningPoints());
      for (std::size_t i = 0; i < hull.numSpanningPoints(); ++i)
        _spanningPoints.push_back(hull.spanningPoint(i));
      _spanningRays.reserve(hull.numSpanningRays());
      for (std::size_t i = 0; i < hull.numSpanningRays(); ++i)
        _spanningRays.push_back(hull.spanningRay(i));
      equations = &_equations;
    }
    else
    {
      orac->setFace(face);
      hull.run(faceEquations, orac, hullOutput, 1, true);
      equations = &faceEquations;
      orac->setFace(completeFace());
    }

    if (_taskDimension)
    {
      std::cout << " Dimension: " << hull.dimension() << std::endl;
    }

    if (_optionReadable)
      manhattanNormImproveEquations(n, *equations);

    if (_taskEquations)
    {
      for (std::size_t i = 0; i < equations->size(); ++i)
      {
        std::cout << " Equation: ";
        oracle()->space().printLinearConstraint(std::cout, (*equations)[i]);
        std::cout << "\n";
      }
      std::cout << std::flush;
    }

    return true;
  }

  bool ConsoleApplicationBase::optimizeObjective(const Vector& objective, bool maximize)
  {
    std::cout << (maximize ? " Maximum: " : " Minimum: ") << std::flush;

    OracleResult result;
    if (maximize)
    {
      oracle()->maximize(result, objective, ObjectiveBound(), 0);
    }
    else
    {
      oracle()->maximize(result, -objective, ObjectiveBound(), 0);
    }

    if (result.isInfeasible())
    {
      std::cout << "infeasible.\n" << std::flush;
    }
    else if (result.isUnbounded())
    {
      std::cout << "unbounded ray ";
      oracle()->space().printVector(std::cout, result.rays.front().vector);
      std::cout << "\n" << std::flush;
    }
    else
    {
      oracle()->space().printVector(std::cout, result.points.front().vector);
      std::cout << " (value: " << result.points.front().objectiveValue << ")\n" << std::flush;
    }
    return true;
  }

  bool ConsoleApplicationBase::generateFacets(const Vector& objective, bool print)
  {
    std::size_t n = oracle()->space().dimension();

    /// Setup the LP.

//     std::cerr << "Creating new SoPlex instance for generateFacets." << std::endl;
    SoPlex* spx = new SoPlex;
    spx->setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    spx->setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    spx->setRealParam(SoPlex::FEASTOL, 0.0);
    spx->setBoolParam(SoPlex::RATREC, true);
    spx->setBoolParam(SoPlex::RATFAC, true);
    spx->setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    spx->setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    LPColSetRational cols;
    cols = _relaxationColumns;
    for (std::size_t v = 0; v < cols.num(); ++v)
      cols.maxObj_w(v) = 0;
    for (std::size_t p = 0; p < objective.size(); ++p)
      cols.maxObj_w(objective.index(p)) = objective.value(p);
    spx->addColsRational(cols);
    addToLP(*spx, _equations);
    addToLP(*spx, _relaxationConstraints);

    Separation::QuietOutput separateOutput;
    Separation::Result separate(_spanningPoints, _spanningRays, _basicColumns, oracle());

    DVectorRational denseSolution;
    DSVectorRational sparseSolution;
    denseSolution.reDim(n);
    while (true)
    {
      SPxSolver::Status status = spx->solve();
      if (status == SPxSolver::UNBOUNDED)
      {
        spx->getPrimalRayRational(denseSolution);
        Vector ray = denseToVector(denseSolution);

//        std::cout << "Relaxation LP is unbounded with ray ";
//        oracle()->printVector(std::cout, &sparseSolution);
//        std::cout << std::endl;

        separate.separateRay(ray, separateOutput);
        if (separate.violation() <= 0)
          break;
      }
      else if (status == SPxSolver::OPTIMAL)
      {
        spx->getPrimalRational(denseSolution);
        Vector point = denseToVector(denseSolution);

//        std::cout << "Relaxation LP is bounded with optimum ";
//        oracle()->printVector(std::cout, &sparseSolution);
//        std::cout << " of value " << spx->objValueRational() << "." << std::endl;

        separate.separatePoint(point, separateOutput);
        if (separate.violation() <= 0)
          break;
      }
      else
      {
        std::stringstream ss;
        ss << "Cut loop LP could not be solved to optimality. Status is " << status << ".";
        throw std::runtime_error(ss.str());
      }

      /// Obtain inequality and certificate.

      Separation::Certificate certificate;
      LinearConstraint inequality = separate.inequality();
      separate.certificate(certificate);
      addToLP(*spx, inequality);

      /// If it should be reused, record the facet.
      if (_optionReuseFacets)
      {
        _relaxationConstraints.push_back(inequality);
      }

      if (separate.separatedFacet())
        std::cout << " Facet: ";
      else if (separate.separatedEquation())
        std::cout << " Equation: ";
      else
      {
        throw std::runtime_error("A bug in IPO occured: Separated neither a facet nor an equation! Please report.");
      }

      manhattanNormImproveInequality(n, inequality, _equations);

      oracle()->space().printLinearConstraint(std::cout, inequality);

      if (_optionCertificates)
      {
        std::cout << ", certified by " << certificate.points.size() << " points and "
            << certificate.rays.size() << " rays.\n" << std::endl;
        for (std::size_t i = 0; i < certificate.points.size(); ++i)
        {
          std::cout << "  Certifying point: ";
          oracle()->space().printVector(std::cout, certificate.points[i]);
          std::cout << "\n";
        }
        for (std::size_t i = 0; i < certificate.rays.size(); ++i)
        {
          std::cout << "  Certifying ray: ";
          oracle()->space().printVector(std::cout, certificate.rays[i]);
          std::cout << "\n";
        }
        std::cout << "\n" << std::flush;
      }
      else
      {
        std::cout << "\n" << std::flush;
      }
    }

    delete spx;

    return true;
  }

  bool ConsoleApplicationBase::separateRayFacet(const Vector& ray, bool& isFeasible)
  {
    std::size_t n = oracle()->space().dimension();

    Separation::QuietOutput separateOutput;
    Separation::Result separate(_spanningPoints, _spanningRays, _basicColumns, oracle());
    separate.separateRay(ray, separateOutput);
    if (separate.violation() <= 0)
    {
      isFeasible = true;
      std::cout << " is feasible.\n" << std::flush;
      return true;
    }
    else
      isFeasible = false;


    Separation::Certificate certificate;
    LinearConstraint inequality = separate.inequality();
    separate.certificate(certificate);

    if (separate.separatedFacet())
      std::cout << "\n Separated by facet: ";
    else if (separate.separatedEquation())
      std::cout << "\n Separated by equation: ";
    else
    {
      throw std::runtime_error("A bug in IPO occured: Separated neither a facet nor an equation! Please report.");
    }

    manhattanNormImproveInequality(n, inequality, _equations);

    oracle()->space().printLinearConstraint(std::cout, inequality);

    if (_optionCertificates)
    {
      std::cout << ", certified by " << certificate.points.size() << " points and "
          << certificate.rays.size() << " rays.\n" << std::flush;
      for (std::size_t i = 0; i < certificate.points.size(); ++i)
      {
        std::cout << "  Certifying point: ";
        oracle()->space().printVector(std::cout, certificate.points[i]);
        std::cout << "\n";
      }
      for (std::size_t i = 0; i < certificate.rays.size(); ++i)
      {
        std::cout << "  Certifying ray: ";
        oracle()->space().printVector(std::cout, certificate.rays[i]);
        std::cout << "\n";
      }
      std::cout << std::flush;
    }
    else
    {
      std::cout << "\n" << std::flush;
    }

    return true;
  }

  bool ConsoleApplicationBase::separatePointFacet(const Vector& point, bool& isFeasible)
  {
    std::size_t n = oracle()->space().dimension();

    Separation::QuietOutput separateOutput;
    Separation::Result separate(_spanningPoints, _spanningRays, _basicColumns, oracle());
    separate.separatePoint(point, separateOutput);
    if (separate.violation() <= 0)
    {
      isFeasible = true;
      std::cout << " is feasible.\n" << std::flush;
      return true;
    }
    else
      isFeasible = false;

    Separation::Certificate certificate;
    LinearConstraint inequality = separate.inequality();
    separate.certificate(certificate);

    if (separate.separatedFacet())
    {
      std::cout << "\n Separated by facet: ";
      manhattanNormImproveInequality(n, inequality, _equations);
    }
    else if (separate.separatedEquation())
      std::cout << "\n Separated by equation: ";
    else
    {
      throw std::runtime_error("A bug in IPO occured: Separated neither a facet nor an equation! Please report.");
    }

    oracle()->space().printLinearConstraint(std::cout, inequality);

    if (_optionCertificates)
    {
      std::cout << ", certified by " << certificate.points.size() << " points and "
          << certificate.rays.size() << " rays.\n" << std::flush;
      for (std::size_t i = 0; i < certificate.points.size(); ++i)
      {
        std::cout << "  Certifying point: ";
        oracle()->space().printVector(std::cout, certificate.points[i]);
        std::cout << "\n";
      }
      for (std::size_t i = 0; i < certificate.rays.size(); ++i)
      {
        std::cout << "  Certifying ray: ";
        oracle()->space().printVector(std::cout, certificate.rays[i]);
        std::cout << "\n";
      }
      std::cout << std::flush;
    }
    else
    {
      std::cout << "\n" << std::flush;
    }

    return true;
  }

  bool ConsoleApplicationBase::computeSmallestFace(const Vector& point)
  {
    SmallestFace::QuietOutput smallestFaceOutput;
    SmallestFace::Result smallestFace(oracle());
    smallestFace.run(point, smallestFaceOutput);

    std::cout << " Dimension: " << smallestFace.dimension() << std::endl;
    Vector maximizingObjective = smallestFace.getMaximizingObjective();
    std::cout << " Objective : ";
    oracle()->space().printVector(std::cout, maximizingObjective);
    std::cout << "\n" << std::flush;

    return true;
  }

  bool ConsoleApplicationBase::printCached()
  {
    std::cout << "Cached points: " << _cacheOracle->numPoints() << "\n";
//     for (std::size_t i = _cachedPoints->first(); i < _cachedPoints->size(); i = _cachedPoints->next(i))
//     {
//       std::cout << ' ';
//       space().printVector(std::cout, _cachedPoints->get(i));
//       std::cout << '\n';
//     }
    std::cout << "Cached rays: " << _cacheOracle->numRays() << "\n";
//     for (std::size_t i = _cachedDirections->first(); i < _cachedDirections->size(); i = _cachedDirections->next(i))
//     {
//       std::cout << ' ';
//       space().printVector(std::cout, _cachedDirections->get(i));
//       std::cout << '\n';
//     }
    std::cout << std::flush;
    return true;
  }

} /* namespace ipo */
