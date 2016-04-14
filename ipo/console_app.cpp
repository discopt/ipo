#include "console_app.h"

#include <cassert>
#include <random>

#include "affine_hull.h"
#include "min_norm_2d.h"
#include "facets.h"
#include "smallest_face.h"
#include "spx_gmp.h"
#include "parser.h"
#include "ipo.h"

using namespace soplex;

namespace ipo {
  
  class ConsoleApplicationObjectiveParser : public LPObjectiveParser
  {
  public:
    ConsoleApplicationObjectiveParser(std::istream& stream, OptimizationOracleBase* oracle, std::vector<DSVectorRational*>& objectives, std::vector<std::string>& objectiveNames)
      : LPObjectiveParser(stream), _objectives(objectives), _objectiveNames(objectiveNames)
    {
      for (std::size_t i = 0; i < oracle->space().dimension(); ++i)
        _oracleVariables[oracle->space()[i]] = i;
    }

    virtual ~ConsoleApplicationObjectiveParser()
    {
      
    }

    virtual void handleObjective(const std::string& name, const std::map<std::string, Rational>& coefficients)
    {
      DSVectorRational* vector = new DSVectorRational;
      for (std::map<std::string, Rational>::const_iterator iter = coefficients.begin(); iter != coefficients.end(); ++iter)
      {
        std::map<std::string, std::size_t>::const_iterator varIter = _oracleVariables.find(iter->first);
        if (varIter != _oracleVariables.end())
        {
          vector->add(varIter->second, iter->second);
        }
        else
        {
          std::cerr << "Skipping objective: Unknown variable <" << iter->first << ">.\n" << std::endl;
          delete vector;
          return;
        }
      }
      _objectives.push_back(vector);
      _objectiveNames.push_back(name);
    }

  private:
    std::vector<DSVectorRational*>& _objectives;
    std::vector<std::string>& _objectiveNames;
    std::map<std::string, std::size_t> _oracleVariables;
  };
  
  class ConsoleApplicationInequalityParser : public LPInequalityParser
  {
  public:
    ConsoleApplicationInequalityParser(std::istream& stream, OptimizationOracleBase* oracle, std::vector<Face*>& faces, std::vector<std::string>& faceNames)
      : LPInequalityParser(stream), _faces(faces), _faceNames(faceNames)
    {
      for (std::size_t i = 0; i < oracle->space().dimension(); ++i)
        _oracleVariables[oracle->space()[i]] = i;
    }

    virtual ~ConsoleApplicationInequalityParser()
    {
      
    }
    
    virtual void handleInequality(const std::string& name, const Rational& lhs, const std::map< std::string, Rational >& coefficients, const Rational& rhs)
    {
      if ((lhs <= -infinity && rhs >= infinity) || lhs == rhs)
        return;

      DSVectorRational vector;
      for (std::map<std::string, Rational>::const_iterator iter = coefficients.begin(); iter != coefficients.end(); ++iter)
      {
        std::map<std::string, std::size_t>::const_iterator varIter = _oracleVariables.find(iter->first);
        if (varIter != _oracleVariables.end())
        {
          vector.add(varIter->second, iter->second);
        }
        else
        {
          std::cerr << "Skipping inequality: Unknown variable <" << iter->first << ">.\n" << std::endl;
          return;
        }
      }
      
      bool hasBoth = lhs > -infinity && rhs < infinity;
      if (lhs > -infinity)
      {
        LPRowRational row(lhs, vector, infinity);
        _faces.push_back(new Face(_oracleVariables.size(), row));
        _faceNames.push_back(hasBoth ? (name + "-lhs") : name);
      }
      if (rhs < infinity)
      {
        LPRowRational row(-infinity, vector, rhs);
        _faces.push_back(new Face(_oracleVariables.size(), row));
        _faceNames.push_back(hasBoth ? (name + "-rhs") : name);
      }
    }
    
  private:
    std::map<std::string, std::size_t> _oracleVariables;
    std::vector<Face*>& _faces;
    std::vector<std::string>& _faceNames;
  };
  
  class ConsoleApplicationPointParser : public PointParser
  {
  public:
    ConsoleApplicationPointParser(std::istream& stream, OptimizationOracleBase* oracle, std::vector<Point*>& points, std::vector<std::string>& pointNames)
      : PointParser(stream), _points(points), _pointNames(pointNames)
    {
      for (std::size_t i = 0; i < oracle->space().dimension(); ++i)
        _oracleVariables[oracle->space()[i]] = i;
    }

    virtual ~ConsoleApplicationPointParser()
    {
      
    }
    
    virtual void handlePoint(const std::string& name, const std::map< std::string, Rational >& values)
    {
      DSVectorRational* point = new DSVectorRational;
      for (std::map<std::string, Rational>::const_iterator iter = values.begin(); iter != values.end(); ++iter)
      {
        std::map<std::string, std::size_t>::const_iterator varIter = _oracleVariables.find(iter->first);
        if (varIter != _oracleVariables.end())
        {
          point->add(varIter->second, iter->second);
        }
        else
        {
          std::cerr << "Skipping point: Unknown variable <" << iter->first << ">.\n" << std::endl;
          delete point;
          return;
        }
      }
      
      _points.push_back(point);
      _pointNames.push_back(name);
    }
    
  private:
    std::map<std::string, std::size_t> _oracleVariables;
    std::vector<Point*>& _points;
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
    _optionPrintRandom = false;

    _oracle = NULL;
    _cachedDirections = NULL;
    _cachedPoints = NULL;
    _equations = NULL;
  }

  ConsoleApplicationBase::~ConsoleApplicationBase()
  {
    if (_cachedDirections)
      delete _cachedDirections;
    if (_cachedPoints)
      delete _cachedPoints;
    if (_equations)
      delete _equations;

    if (_oracle)
      delete _oracle;
    for (std::size_t i = 0; i < numObjectives(); ++i)
      delete _objectives[i];
    for (std::size_t i = 0; i < numFaces(); ++i)
      delete _faces[i];
    for (std::size_t i = 0; i < numDirections(); ++i)
      delete _directions[i];
    for (std::size_t i = 0; i < numPoints(); ++i)
      delete _points[i];
  }

  void ConsoleApplicationBase::setBasicOracle(FaceOptimizationOracleBase* oracle)
  {
    if (_oracle)
      throw std::runtime_error("Error in ConsoleApplicationBase::setBasicOracle: Oracle already set.");
    if (_faceRestrictionArgument != "" || _projectionArgument != "")
      throw std::runtime_error("Restricting to faces or projecting is not implemented, yet.");
    _oracle = oracle;
    _space = oracle->space();
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
      if (!computeAffineHull(NULL, ""))
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
          std::cout << "Random objective";
        else
          std::cout << "Objective <" << _objectiveNames[i] << ">";
        bool print = (_objectiveNames[i] != ":RANDOM:") || _optionPrintRandom;
        if (print)
        {
          std::cout << " ";
          space().printLinearForm(std::cout, objective(i));
        }
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
      for (std::size_t i = 0; i < numDirections(); ++i)
      {
        std::cout << "Point";
        if (_directionNames[i] != "")
          std::cout << " " << _directionNames[i];
        std::cout << std::flush;
        
        bool isFeasible = true;
        if (!separateDirectionFacet(directions(i), isFeasible))
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
          if (!separatePointFacet(_points[i], isFeasible))
            return false;
        }
        else
        {
          std::cout << " is assumed to be feasible.\n" << std::flush;
        }
        if (taskSmallestFace() && isFeasible)
        {
          if (!computeSmallestFace(points(i)))
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
    std::cerr << "      Consider the orthogonal projection onto a subset of variables. VARS is either a comma-\n";
    std::cerr << "      separated list of variables or a regular expression enclosed with ^ and $. In this case it\n";
    std::cerr << "      defines the projection onto all variables matching it.\n";
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
    std::cerr << "      Attempt to separate every given point or direction using a facet-defining inequality.\n";
    std::cerr << "  --smallest-face\n";
    std::cerr << "      Compute, for every given point in P, an objective vector that induces the smallest face\n";
    std::cerr << "      containing the point, and the dimension of that face.\n";
    std::cerr << "  --facets\n";
    std::cerr << "      Generate facets that are useful when maximizing the objectives specified by --objective(s) or\n";
    std::cerr << "      generated by --random.\n";
    std::cerr << "  --print-cached\n";
    std::cerr << "      Print the cached points and directions at the end.\n";
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
    std::cerr << "  --direction DIR\n";
    std::cerr << "      Use the given unbounded direction for facet-separation. DIR is a sparse vector, e.g.,\n";
    std::cerr << "      \"(x#1, 3y#5)\". Usable multiple times.\n";
    std::cerr << "  --directions FILE\n";
    std::cerr << "      The content of FILE is considered as a --direction argument. If they do not match the pattern,\n";
    std::cerr << "      tokens are silently ignored.\n";
    printAdditionalOptionsInput(std::cerr);
    std::cerr << "\n";
    std::cerr << "Further options (ordering does not matter):\n";
    std::cerr << "  --cache on|off\n";
    std::cerr << "      If on (default), cache points and directions produced by IPO to speed up further computations.\n";
    std::cerr << "  --readable on|off\n";
    std::cerr << "      If on (default), improve readability of inequalities/equations produced by IPO.\n";
    std::cerr << "  --certificates on|off\n";
    std::cerr << "      If on, output, for every computed facet, points and directions that span it. Off by default.\n";
    std::cerr << "  --trust-points on|off\n";
    std::cerr << "      If off (default), it is checked that points are in P when computing the smallest containing\n";
    std::cerr << "      faces.\n";
    std::cerr << "  --reuse-facets on|off\n";
    std::cerr << "      If on (default), use computed facets for subsequent facet-generation.\n";
    std::cerr << "  --print-random on|off\n";
    std::cerr << "      If off (default), do not print the random objectives.\n";
    std::cerr << "  --progress stdout|stderr|off\n";
    std::cerr << "      If on, print progress output to stdout or stderr. Off by default.\n";
    std::cerr << "  --debug stdout|stderr|off\n";
    std::cerr << "      If on, print debug output to stdout or stderr. Off by default.\n";
    printAdditionalOptionsFurther(std::cerr);
    printAdditionalOptionsSpecific(std::cerr);
    std::cerr << std::flush;
  }
  
  void ConsoleApplicationBase::setRelaxationBounds(const soplex::VectorRational& lowerBounds, const soplex::VectorRational& upperBounds)
  {
    assert(lowerBounds.dim() == _relaxationColumns.num());
    assert(upperBounds.dim() == _relaxationColumns.num());
    
    _relaxationColumns.lower_w() = lowerBounds;
    _relaxationColumns.upper_w() = upperBounds;
  }
  
  void ConsoleApplicationBase::addRelaxationRows(const soplex::LPRowSetRational& rows)
  {
    _relaxationRows.add(rows);
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
    if (firstArgument == "--direction")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --direction must be followed by a vector.");
      }
      _directionArguments.push_back(argument(1));
      return 2;
    }
    if (firstArgument == "--directions")
    {
      if (numArguments() == 1)
      {
        throw std::runtime_error("Invalid option: --directions must be followed by a file name.");
      }
      _directionFiles.push_back(argument(1));
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
      throw std::runtime_error("Invalid option: --readable must be followed by either \"on\" or \"off\".");
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
        if (argument(1) == "on")
        {
          _optionPrintRandom = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _optionPrintRandom = false;
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

    std::size_t n = space().dimension();

    _cachedDirections = new UniqueRationalVectors(n);
    _cachedPoints = new UniqueRationalVectors(n);
    _equations = new LPRowSetRational;
    _relaxationColumns.clear();
    DSVectorRational zeroVector;
    for (std::size_t c = 0; c < n; ++c)
      _relaxationColumns.add(0, -infinity, zeroVector, infinity);
    
    /// Add specified objectives.
    
    for (std::size_t i = 0; i < _objectiveArguments.size(); ++i)
    {
      std::istringstream stream(_objectiveArguments[i]);
      ConsoleApplicationObjectiveParser parser(stream, oracle(), _objectives, _objectiveNames);
      parser.run();
    }
    
    /// Parse objectives from specified files.
    
    for (std::size_t i = 0; i < _objectiveFiles.size(); ++i)
    {
      std::ifstream stream(_objectiveFiles[i]);
      ConsoleApplicationObjectiveParser parser(stream, oracle(), _objectives, _objectiveNames);
      parser.run();
    }

    /// Add random objectives.

    if (_numRandomObjectives > 0)
    {
      std::random_device randomDevice;
      std::default_random_engine generator(randomDevice());
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
          DSVectorRational* objectiveVector = new DSVectorRational;
          for (std::size_t c = 0; c < n; ++c)
          {
            if (randomVector[c] != 0.0)
              objectiveVector->add(c, Rational(randomVector[c] / norm));
          }
          _objectives.push_back(objectiveVector);
          _objectiveNames.push_back(":RANDOM:");
        }
      }
    }
    
    /// Parse specified faces.

    for (std::size_t i = 0; i < _faceArguments.size(); ++i)
    {
      std::istringstream stream(_faceArguments[i]);
      ConsoleApplicationInequalityParser parser(stream, oracle(), _faces, _faceNames);
      parser.run();
    }
    
    /// Parse faces from specified files.
    
    for (std::size_t i = 0; i < _faceFiles.size(); ++i)
    {
      std::ifstream stream(_faceFiles[i]);
      ConsoleApplicationInequalityParser parser(stream, oracle(), _faces, _faceNames);
      parser.run();
    }
    
    /// Add specified directions.

    for (std::size_t i = 0; i < _directionArguments.size(); ++i)
    {
      std::istringstream stream(_directionArguments[i]);
      ConsoleApplicationPointParser parser(stream, oracle(), _directions, _directionNames);
      parser.run();
    }
    
    /// Parse directions from specified files.
    
    for (std::size_t i = 0; i < _directionFiles.size(); ++i)
    {
      std::ifstream stream(_directionFiles[i]);
      ConsoleApplicationPointParser parser(stream, oracle(), _directions, _directionNames);
      parser.run();
    }
    
     /// Add specified points.

    for (std::size_t i = 0; i < _pointArguments.size(); ++i)
    {
      std::istringstream stream(_pointArguments[i]);
      ConsoleApplicationPointParser parser(stream, oracle(), _points, _pointNames);
      parser.run();
    }

    /// Parse points from specified files.
    
    for (std::size_t i = 0; i < _pointFiles.size(); ++i)
    {
      std::ifstream stream(_pointFiles[i]);
      ConsoleApplicationPointParser parser(stream, oracle(), _points, _pointNames);
      parser.run();
    }

    return true;
  }

  bool ConsoleApplicationBase::printAmbientDimension()
  {
    std::cout << "Ambient dimension:\n " << space().dimension() << "\n" << std::flush;
    return true;
  }

  bool ConsoleApplicationBase::printVariables()
  {
    std::cout << "Variables:\n";
    for (std::size_t i = 0; i < space().dimension(); ++i)
    {
      const std::string& name = space()[i];
      std::cout << " " << name << (i + 1 < space().dimension() ? "," : "\n");
    }
    return true;
  }

  bool ConsoleApplicationBase::computeAffineHull(Face* face, const std::string& faceName)
  {
    if (_taskDimension || _taskEquations)
    {
      std::cout << "Computing the affine hull";
      if (face)
      {
        std::cout << " of face " << faceName << ":\n" << std::flush;
      }
      else
      {
        std::cout << ":\n" << std::flush;
      }
    }

    FaceOptimizationOracleBase* orac = oracle();
    std::size_t n = space().dimension();

    AffineHull::QuietOutput hullOutput;
    AffineHull::Result hull;
    
    LPRowSetRational faceEquations;
    
    LPRowSetRational* equations;
    if (face == NULL)
    {
      hull.run(*_cachedPoints, *_cachedDirections, *_equations, orac, hullOutput, AffineHull::REDUNDANT_EQUATIONS_REMOVE);
      _basicColumns = hull.basicColumns();
      _spanningCachedDirections = hull.spanningDirections();
      _spanningCachedPoints = hull.spanningPoints();
      equations = _equations;
    }
    else
    {
      orac->setFace(face);
      FilteredUniqueRationalVectors filteredPoints(*_cachedPoints);
      FilteredUniqueRationalVectors filteredDirections(*_cachedDirections);

      /// Filter points and rays.

      soplex::DVectorRational faceNormal(n);
      faceNormal.clear();
      faceNormal.assign(face->normal());

      for (std::size_t p = _cachedPoints->first(); p < _cachedPoints->size(); p = _cachedPoints->next(p))
      {
        soplex::Rational activity = *(*_cachedPoints)[p] * faceNormal;
        filteredPoints.set(p, activity == face->rhs());
      }
      for (std::size_t r = _cachedDirections->first(); r < _cachedDirections->size(); r = _cachedDirections->next(r))
      {
        soplex::Rational activity = *(*_cachedDirections)[r] * faceNormal;
        filteredDirections.set(r, activity == 0);
      }
      
      faceEquations = LPRowSetRational(*_equations);
      faceEquations.add(face->rhs(), face->normal(), face->rhs());
      
      hull.run(filteredPoints, filteredDirections, faceEquations, orac, hullOutput, AffineHull::REDUNDANT_EQUATIONS_REMOVE);
      
      equations = &faceEquations;
      orac->setFace();
    }

    if (_taskDimension)
    {
      std::cout << " Dimension: " << hull.dimension() << std::endl;
    }

    if (_optionReadable)
      manhattanNormImproveEquations(n, *_equations);

    if (_taskEquations)
    {
      for (int i = 0; i < equations->num(); ++i)
      {
        std::cout << " Equation: ";
        space().printRow(std::cout, *equations, i);
        std::cout << "\n";
      }
      std::cout << std::flush;
    }

    return true;
  }

  bool ConsoleApplicationBase::optimizeObjective(const SVectorRational* objective, bool maximize)
  {
    std::cout << (maximize ? " Maximum: " : " Minimum: ") << std::flush;
    
    OptimizationResult result;
    if (maximize)
      oracle()->maximize(result, *objective);
    else
    {
      DSVectorRational negated;
      negated = *objective;
      for (int p = negated.size() - 1; p >= 0; --p)
        negated.value(p) *= -1;
      oracle()->maximize(result, negated);
    }
    
    if (result.isInfeasible())
    {
      std::cout << "infeasible.\n" << std::flush;
    }
    else if (result.isUnbounded())
    {
      std::cout << "unbounded direction ";
      space().printVector(std::cout, result.directions.front());
      std::cout << "\n" << std::flush;
    }
    else
    {
      space().printVector(std::cout, result.points.front());
      std::cout << "\n" << std::flush;
    }
    return true;
  }

  bool ConsoleApplicationBase::generateFacets(const SVectorRational* objective, bool print)
  {
    std::size_t n = space().dimension();

    /// Setup the LP.

    SoPlex spx;
    spx.setIntParam(SoPlex::SOLVEMODE, SoPlex::SOLVEMODE_RATIONAL);
    spx.setIntParam(SoPlex::SYNCMODE, SoPlex::SYNCMODE_AUTO);
    spx.setRealParam(SoPlex::FEASTOL, 0.0);
    spx.setBoolParam(SoPlex::RATREC, true);
    spx.setBoolParam(SoPlex::RATFAC, true);
    spx.setIntParam(SoPlex::OBJSENSE, SoPlex::OBJSENSE_MAXIMIZE);
    spx.setIntParam(SoPlex::VERBOSITY, SoPlex::VERBOSITY_ERROR);
    LPColSetRational cols;
    cols = _relaxationColumns;
    cols.maxObj_w().assign(*objective);
    spx.addColsRational(cols);
    spx.addRowsRational(*_equations);
    spx.addRowsRational(_relaxationRows);
    
    Separation::QuietOutput separateOutput;
    Separation::Result separate(*_cachedPoints, *_cachedDirections, _spanningCachedPoints, _spanningCachedDirections,
        _basicColumns, oracle());

    DVectorRational denseSolution;
    DSVectorRational sparseSolution;
    denseSolution.reDim(n);
    while (true)
    {
      SPxSolver::Status status = spx.solve();
      if (status == SPxSolver::UNBOUNDED)
      {
        spx.getPrimalRayRational(denseSolution);
        sparseSolution.clear();
        sparseSolution = denseSolution;

//        std::cout << "Relaxation LP is unbounded with ray ";
//        oracle()->printVector(std::cout, &sparseSolution);
//        std::cout << std::endl;

        separate.separateRay(&sparseSolution, separateOutput);
        if (separate.violation() <= 0)
          break;
      }
      else if (status == SPxSolver::OPTIMAL)
      {
        spx.getPrimalRational(denseSolution);
        sparseSolution.clear();
        sparseSolution = denseSolution;

//        std::cout << "Relaxation LP is bounded with optimum ";
//        oracle()->printVector(std::cout, &sparseSolution);
//        std::cout << " of value " << spx.objValueRational() << "." << std::endl;

        separate.separatePoint(&sparseSolution, separateOutput);
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

      LPRowRational inequality;
      Separation::Certificate certificate;
      separate.inequality(inequality);
      separate.certificate(certificate);
      spx.addRowRational(inequality);
      
      /// If it should be reused, record the facet.
      if (_optionReuseFacets)
      {
        _relaxationRows.add(inequality);
      }

      if (separate.separatedFacet())
        std::cout << " Facet: ";
      else if (separate.separatedEquation())
        std::cout << " Equation: ";
      else
      {
        throw std::runtime_error("A bug in IPO occured: Separated neither a facet nor an equation! Please report.");
      }

      manhattanNormImproveInequality(n, inequality, *_equations);

      space().printRow(std::cout, inequality);

      if (_optionCertificates)
      {
        std::cout << ", certified by " << certificate.pointIndices.size() << " points and "
            << certificate.directionIndices.size() << " rays.\n" << std::endl;
        for (std::size_t i = 0; i < certificate.pointIndices.size(); ++i)
        {
          std::cout << "  Certifying point: ";
          space().printVector(std::cout, (*_cachedPoints)[certificate.pointIndices[i]]);
          std::cout << "\n";
        }
        for (std::size_t i = 0; i < certificate.directionIndices.size(); ++i)
        {
          std::cout << "  Certifying ray: ";
          space().printVector(std::cout, (*_cachedDirections)[certificate.directionIndices[i]]);
          std::cout << "\n";
        }
        std::cout << "\n" << std::flush;
      }
      else
      {
        std::cout << "\n" << std::flush;
      }
    }

    return true;
  }

  bool ConsoleApplicationBase::separateDirectionFacet(const Direction* direction, bool& isFeasible)
  {
    std::size_t n = space().dimension();
    
    Separation::QuietOutput separateOutput;
    Separation::Result separate(*_cachedPoints, *_cachedDirections, _spanningCachedPoints, _spanningCachedDirections,
        _basicColumns, oracle());
    separate.separateRay(direction, separateOutput);
    if (separate.violation() <= 0)
    {
      isFeasible = true;
      std::cout << " is feasible.\n" << std::flush;
      return true;
    }
    else
      isFeasible = false;
 
 
    LPRowRational inequality;
    Separation::Certificate certificate;
    separate.inequality(inequality);
    separate.certificate(certificate);
    
    if (separate.separatedFacet())
      std::cout << "\n Separated by facet: ";
    else if (separate.separatedEquation())
      std::cout << "\n Separated by equation: ";
    else
    {
      throw std::runtime_error("A bug in IPO occured: Separated neither a facet nor an equation! Please report.");
    }
    
    manhattanNormImproveInequality(n, inequality, *_equations);

    space().printRow(std::cout, inequality);

    if (_optionCertificates)
    {
      std::cout << ", certified by " << certificate.pointIndices.size() << " points and "
          << certificate.directionIndices.size() << " rays.\n" << std::flush;
      for (std::size_t i = 0; i < certificate.pointIndices.size(); ++i)
      {
        std::cout << "  Certifying point: ";
        space().printVector(std::cout, (*_cachedPoints)[certificate.pointIndices[i]]);
        std::cout << "\n";
      }
      for (std::size_t i = 0; i < certificate.directionIndices.size(); ++i)
      {
        std::cout << "  Certifying direction: ";
        space().printVector(std::cout, (*_cachedDirections)[certificate.directionIndices[i]]);
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

  bool ConsoleApplicationBase::separatePointFacet(const Point* point, bool& isFeasible)
  {
    std::size_t n = space().dimension();
    
    Separation::QuietOutput separateOutput;
    Separation::Result separate(*_cachedPoints, *_cachedDirections, _spanningCachedPoints, _spanningCachedDirections,
        _basicColumns, oracle());
    separate.separatePoint(point, separateOutput);
    if (separate.violation() <= 0)
    {
      isFeasible = true;
      std::cout << " is feasible.\n" << std::flush;
      return true;
    }
    else
      isFeasible = false;

    LPRowRational inequality;
    Separation::Certificate certificate;
    separate.inequality(inequality);
    separate.certificate(certificate);
    
    if (separate.separatedFacet())
      std::cout << "\n Separated by facet: ";
    else if (separate.separatedEquation())
      std::cout << "\n Separated by equation: ";
    else
    {
      throw std::runtime_error("A bug in IPO occured: Separated neither a facet nor an equation! Please report.");
    }

    manhattanNormImproveInequality(n, inequality, *_equations);

    space().printRow(std::cout, inequality);

    if (_optionCertificates)
    {
      std::cout << ", certified by " << certificate.pointIndices.size() << " points and "
          << certificate.directionIndices.size() << " rays.\n" << std::flush;
      for (std::size_t i = 0; i < certificate.pointIndices.size(); ++i)
      {
        std::cout << "  Certifying point: ";
        space().printVector(std::cout, (*_cachedPoints)[certificate.pointIndices[i]]);
        std::cout << "\n";
      }
      for (std::size_t i = 0; i < certificate.directionIndices.size(); ++i)
      {
        std::cout << "  Certifying direction: ";
        space().printVector(std::cout, (*_cachedDirections)[certificate.directionIndices[i]]);
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

  bool ConsoleApplicationBase::computeSmallestFace(const Point* point)
  {
    SmallestFace::QuietOutput smallestFaceOutput;
    SmallestFace::Result smallestFace(*_cachedPoints, *_cachedDirections, oracle());
    smallestFace.run(point, smallestFaceOutput);

    std::cout << " Dimension: " << smallestFace.dimension() << std::endl;
    Point maximizingObjective;
    smallestFace.getMaximizingObjective(maximizingObjective);
    std::cout << " Objective : ";
    space().printVector(std::cout, &maximizingObjective);
    std::cout << "\n" << std::flush;
    
    return true;
  }

  bool ConsoleApplicationBase::printCached()
  {
    std::cout << "Cached points: " << _cachedPoints->size() << "\n";
    for (std::size_t i = _cachedPoints->first(); i < _cachedPoints->size(); i = _cachedPoints->next(i))
    {
      std::cout << ' ';
      space().printVector(std::cout, _cachedPoints->get(i));
      std::cout << '\n';
    }
    std::cout << "Cached directions: " << _cachedDirections->size() << "\n";
    for (std::size_t i = _cachedDirections->first(); i < _cachedDirections->size(); i = _cachedDirections->next(i))
    {
      std::cout << ' ';
      space().printVector(std::cout, _cachedDirections->get(i));
      std::cout << '\n';
    }
    std::cout << std::flush;
    return true;
  }

} /* namespace ipo */
