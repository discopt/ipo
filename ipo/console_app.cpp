#include "console_app.h"

#include <cassert>
#include <random>

#include "affine_hull.h"
#include "min_norm_2d.h"
#include "facets.h"
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
      for (std::size_t i = 0; i < oracle->numVariables(); ++i)
        _oracleVariables[oracle->variableName(i)] = i;
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
      for (std::size_t i = 0; i < oracle->numVariables(); ++i)
        _oracleVariables[oracle->variableName(i)] = i;
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

    for (std::size_t i = 0; i < numObjectives(); ++i)
    {
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
        bool print = (_objectiveNames[i] != ":RANDOM:") || _optionPrintRandom;
        if (!generateFacets(objective(i), _objectiveNames[i], print))
          return false;
      }
    }

    for (std::size_t i = 0; i < numDirections(); ++i)
    {
      if (taskFacet())
      {
        if (!separateDirectionFacet(directions(i)))
          return false;
      }
    }

    for (std::size_t i = 0; i < numPoints(); ++i)
    {
      if (taskFacet())
      {
        if (!separatePointFacet(points(i)))
          return false;
      }
      if (taskSmallestFace())
      {
        if (!computeSmallestFace(points(i)))
          return false;
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
    std::cerr << " --projection VARS            Consider the orthogonal projection onto a subset of variables.\n";
    std::cerr << "                              VARS is either a comma-separated list of variables or a regular\n";
    std::cerr << "                              expression enclosed with ^ and $. In this case it defines the\n";
    std::cerr << "                              projection onto all variables matching it.\n";
    std::cerr << " --restrict-face INEQ         Restrict the oracle to the specified face. See --face below.\n";
    std::cerr << "\n";
    std::cerr << "Options for tasks (execution order is independent of argument order):\n";
    std::cerr << " --print-ambient-dimension    Print the ambient dimension.\n";
    std::cerr << " --print-variables            Print the names of all variables in a comma-separated list.\n";
    std::cerr << " --maximize                   Maximize objectives specified by --objective(s) arguments.\n";
    std::cerr << " --minimize                   Minimize objectives specified by --objective(s) arguments.\n";
    std::cerr << " --dimension                  Compute the dimension.\n";
    std::cerr << " --equations                  Compute valid equations.\n";
    std::cerr << " --facet                      Attempt to separate every given point or direction using a\n";
    std::cerr << "                              facet-defining inequality.\n";
    std::cerr
        << " --smallest-face              Computes, for every given point in P, an objective vector that\n";
    std::cerr << "                              induces the smallest face containing the point.\n";
    std::cerr << " --facets                     Generate facets that are useful when maximizing the objectives\n";
    std::cerr << "                              specified by --objective(s) or generated by --random.\n";
    std::cerr << " --print-cached               Print the cached points and directions at the end.\n";
    printAdditionalOptionsTasks(std::cerr);
    std::cerr << "\n";
    std::cerr << "Options for input data (ordering w.r.t. tasks does not matter):\n";
    std::cerr << " --face INEQ                  Compute dimension/equations also for face induced by INEQ.\n";
    std::cerr << "                              INEQ is an inequality, potentially preceeded by a name with a colon,\n";
    std::cerr << "                              e.g., \"myface: x#1 + 3y#5 <= 7\". Usable multiple times.\n";
    std::cerr
        << " --faces FILE                 Each line of the file is considered as a single --face argument. If\n";
    std::cerr
        << "                              it does not match the pattern, the line is silently ignored, allowing\n";
    std::cerr << "                              it to be in LP format.\n";
    std::cerr << " --objective OBJ              Use the given objective vector for optimization/facet generation.\n";
    std::cerr << "                              OBJ is an LP file objective, e.g., \"max x#1 + 3y#5\".\n";
    std::cerr << "                              Usable multiple times.\n";
    std::cerr
        << " --objectives FILE            The file is considered as multiple --objective arguments.\n";
    std::cerr << "                              Invalid lines are silently ignored. See LP format.\n";
    std::cerr
        << " --random n                   Sample n objective vectors (uniformly at random from the sphere) to\n";
    std::cerr << "                              be used for facet generation.\n";
    std::cerr << " --point POINT                Use the point for smallest-face computation and facet-separation.\n";
    std::cerr
        << "                              POINT is a sparse vector, e.g., \"(x#1, 3y#5)\". Usable multiple times.\n";
    std::cerr << " --points FILE                Each line of the file is considered as a single --point argument.\n";

    std::cerr << " --direction DIR              Use the given unbounded direction for facet-separation.\n";
    std::cerr
        << "                              DIR is a sparse vector, e.g., \"(x#1, 3y#5)\". Usable multiple times.\n";
    std::cerr
        << " --directions FILE            Each line of the file is considered as a single --direction argument.\n";
    printAdditionalOptionsInput(std::cerr);
    std::cerr << "\n";
    std::cerr << "Further options (ordering does not matter):\n";
    std::cerr << " --cache on|off               Cache points and directions produced by IPO to speed up further\n";
    std::cerr << "                              computations. On by default.\n";
    std::cerr << " --readable on|off            Improve readability of inequalities/equations produced by IPO.\n";
    std::cerr << "                              On by default.\n";
    std::cerr
        << " --certificates on|off        Output, for every computed facet, points and directions that span it.\n";
    std::cerr << "                              Off by default.\n";
    std::cerr << " --trust-points on|off        If off (default), it is checked that points are in P when computing\n";
    std::cerr << "                              the smallest containing faces.\n";
    std::cerr << " --reuse-facets on|off        If on (default), it uses computed facets for subsequent facet-\n";
    std::cerr << "                              generation.\n";
    std::cerr << " --print-random on|off        If off (default), it does not print the random objectives.\n";
    std::cerr << " --progress stdout|stderr|off Prints progress output to stdout or stderr. Off by default.\n";
    std::cerr << " --debug stdout|stderr|off    Prints debug output to stdout or stderr. Off by default.\n";
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

    std::size_t n = oracle()->numVariables();

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

    return true;
  }

  bool ConsoleApplicationBase::printAmbientDimension()
  {
    std::cout << "Ambient dimension:\n " << _oracle->numVariables() << "\n" << std::flush;
    return true;
  }

  bool ConsoleApplicationBase::printVariables()
  {
    std::cout << "Variables:\n";
    for (std::size_t i = 0; i < _oracle->numVariables(); ++i)
    {
      const std::string& name = _oracle->variableName(i);
      std::cout << " " << name << (i + 1 < _oracle->numVariables() ? "," : "\n");
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
    std::size_t n = orac->numVariables();

    AffineHull::QuietOutput hullOutput;
    AffineHull::Result hull;
    
    LPRowSetRational properFaceEquations;
    
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
      properFaceEquations = LPRowSetRational(*_equations);
      
      // TODO: init from filtered!
      
      UniqueRationalVectors points(n);
      UniqueRationalVectors directions(n);
      
      hull.run(points, directions, properFaceEquations, orac, hullOutput, AffineHull::REDUNDANT_EQUATIONS_REMOVE);
      
      equations = &properFaceEquations;
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
        orac->printRow(std::cout, *equations, i);
        std::cout << "\n";
      }
      std::cout << std::flush;
    }

    return true;
  }

  bool ConsoleApplicationBase::optimizeObjective(const SVectorRational* objective, bool maximize)
  {
    return true;
  }

  bool ConsoleApplicationBase::generateFacets(const SVectorRational* objective, const std::string& objectiveName, bool print)
  {
    if (objectiveName == "")
      std::cout << "Generating facets for objective";
    else if (objectiveName == ":RANDOM:")
      std::cout << "Generating facets for random objective";
    else
      std::cout << "Generating facets for objective <" << objectiveName << ">";
    if (print)
    {
      std::cout << " ";
      oracle()->printLinearForm(std::cout, objective);
    }
    std::cout << ".\n" << std::flush;
    std::size_t n = oracle()->numVariables();

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
//     Separation::ProgressOutput separateOutput;
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

      oracle()->printRow(std::cout, inequality);

      if (_optionCertificates)
      {
        std::cout << ", certified by " << certificate.pointIndices.size() << " points and "
            << certificate.directionIndices.size() << " rays.\n" << std::endl;
        for (std::size_t i = 0; i < certificate.pointIndices.size(); ++i)
        {
          std::cout << "  Certifying point: ";
          oracle()->printVector(std::cout, (*_cachedPoints)[certificate.pointIndices[i]]);
          std::cout << "\n";
        }
        for (std::size_t i = 0; i < certificate.directionIndices.size(); ++i)
        {
          std::cout << "  Certifying ray: ";
          oracle()->printVector(std::cout, (*_cachedDirections)[certificate.directionIndices[i]]);
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

  bool ConsoleApplicationBase::separateDirectionFacet(const Direction* direction)
  {
    return true;
  }

  bool ConsoleApplicationBase::separatePointFacet(const Point* point)
  {
    return true;
  }

  bool ConsoleApplicationBase::computeSmallestFace(const Point* point)
  {
    return true;
  }

  bool ConsoleApplicationBase::printCached()
  {
    return true;
  }

} /* namespace ipo */
