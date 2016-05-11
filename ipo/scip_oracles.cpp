#include "scip_oracles.h"

#include <limits>
#include <unistd.h>
#include <sys/stat.h>

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include "scip_exception.hpp"
  #define NDEBUG
#else
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include "scip_exception.hpp"
#endif

#include "reconstruct.h"
#include "cpu_timer.h"

using namespace soplex;

/// TODO: Handle unbounded instances.

namespace ipo {

  void getSCIPObjective(SCIP* scip, soplex::DVectorRational& objective, bool makeMaximization)
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    objective.reDim(n);
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    for (std::size_t v = 0; v < n; ++v)
      reconstruct(SCIPvarGetObj(vars[v]), objective[v], SCIPfeastol(scip));

    if (makeMaximization && SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE)
    {
      for (std::size_t v = 0; v < n; ++v)
        objective[v] *= -1;
    }
  }

  void getSCIPvarToIndexMap(SCIP* scip, SCIPvarToIndexMap& map)
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    for (std::size_t v = 0; v < n; ++v)
      map[vars[v]] = v;
    if (SCIPisTransformed(scip))
    {
      SCIP_VAR* transVar = NULL;
      for (std::size_t v = 0; v < n; ++v)
      {
        SCIP_CALL_EXC(SCIPgetTransformedVar(scip, vars[v], &transVar));
        map[transVar] = v;
      }
    }
  }

  SCIPOracle::SCIPOracle(const std::string& name, const Space& space, SCIP* originalSCIP)
    : OracleBase(name, space), _maxLargestCoefficient(1000), _minLargestCoefficient(1),
    _bestLargestCoefficient(1000)
  {
    Space scipSpace;
    copySCIP(originalSCIP, scipSpace);
    if (space != scipSpace)
      throw std::runtime_error("Spaces differ while constructing SCIPOracle.");

    initializedSpace();
  }

  SCIPOracle::SCIPOracle(const std::string& name, Space& space, SCIP* originalSCIP)
    : OracleBase(name, space), _maxLargestCoefficient(1000), _minLargestCoefficient(1),
    _bestLargestCoefficient(1000)
  {
    Space scipSpace;
    copySCIP(originalSCIP, scipSpace);
    if (space.dimension() == 0)
      space = scipSpace;
    else if (space != scipSpace)
      throw std::runtime_error("Spaces differ while constructing SCIPOracle.");

    initializedSpace();
  }

  SCIPOracle::SCIPOracle(const std::string& name, OracleBase* nextOracle, SCIP* originalSCIP)
    : OracleBase(name, nextOracle), _maxLargestCoefficient(1000), _minLargestCoefficient(1),
    _bestLargestCoefficient(1000)
  {
    Space scipSpace;
    copySCIP(originalSCIP, scipSpace);
    if (_space != scipSpace)
      throw std::runtime_error("Spaces differ while constructing SCIPOracle.");

    initializedSpace();
  }

  SCIPOracle::SCIPOracle(const std::string& name, const MixedIntegerProgram& mip)
    : OracleBase(name, mip.space()), _maxLargestCoefficient(1000), _minLargestCoefficient(1),
    _bestLargestCoefficient(1000)
  {
    initializeFromMIP(mip);

    initializedSpace();
  }

  SCIPOracle::SCIPOracle(const std::string& name, OracleBase* nextOracle,
    const MixedIntegerProgram& mip)
    : OracleBase(name, nextOracle)
  {
    if (mip.space() != nextOracle->space())
      throw std::runtime_error("Spaces differ while constructing SCIPOracle.");

    initializeFromMIP(mip);

    initializedSpace();
  }

  SCIPOracle::~SCIPOracle()
  {
    for (std::size_t c = 0; c < _variables.size(); ++c)
    {
      SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_variables[c]));
    }
    SCIP_CALL_EXC(SCIPfree(&_scip));
  }

  void SCIPOracle::copySCIP(SCIP* originalSCIP, Space& scipSpace)
  {
    if (!SCIPisTransformed(originalSCIP))
    {
      throw std::runtime_error(
        "Initialization of SCIPOracle requires transformed problem!");
    }

    std::size_t n = SCIPgetNOrigVars(originalSCIP);
    SCIP_VAR** origVars = SCIPgetOrigVars(originalSCIP);

    for (std::size_t v = 0; v < n; ++v)
    {
      const std::string varName = SCIPvarGetName(origVars[v]);
      scipSpace.addVariable(varName);
    }

    /// Create SCIP instance via copy.

    unsigned int validSCIP = 0;
    SCIP_HASHMAP* hashMap = NULL;
    SCIP_CALL_EXC(SCIPhashmapCreate(&hashMap, SCIPblkmem(originalSCIP), n));
    SCIP_CALL_EXC(SCIPcreate(&_scip));
    SCIP_CALL_EXC(SCIPcopy(originalSCIP, _scip, hashMap, NULL, "-oracle", TRUE, FALSE, FALSE,
      &validSCIP));
    if (!validSCIP)
      throw std::runtime_error("SCIPcopy failed while constructing oracle!");
    SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
    SCIP_CALL_EXC(SCIPsetBoolParam(_scip, "misc/catchctrlc", 0));

    _variables.resize(n);
    for (std::size_t v = 0; v < n; ++v)
    {
      SCIP_VAR* transVar = NULL;
      SCIP_CALL_EXC(SCIPgetTransformedVar(originalSCIP, origVars[v], &transVar));
      _variables[v] = static_cast<SCIP_VAR*>(SCIPhashmapGetImage(hashMap, transVar));
    }
    SCIPhashmapFree(&hashMap);
  }

  void SCIPOracle::initializeFromMIP(const MixedIntegerProgram& mip)
  {
    std::size_t n = mip.space().dimension();

    /// Initialize SCIP.

    SCIP_CALL_EXC(SCIPcreate(&_scip));
    SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
    SCIP_CALL_EXC(SCIPcreateProbBasic(_scip, _name.c_str()));
    SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
    SCIP_CALL_EXC(SCIPsetBoolParam(_scip, "misc/catchctrlc", 0));
    SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0));

    /// Create variables.

    const LPColSetRational& cols = mip.columns();
    _variables.resize(cols.num());
    for (std::size_t c = 0; c < cols.num(); ++c)
    {
      SCIP_CALL_EXC(SCIPcreateVarBasic(_scip, &_variables[c], space()[c].c_str(),
        double(cols.lower(c)), double(cols.upper(c)), double(cols.maxObj(c)),
        mip.isIntegral(c) ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_CONTINUOUS));
      SCIP_CALL_EXC(SCIPaddVar(_scip, _variables[c]));
    }

    /// Create constraints.

    const LPRowSetRational& rows = mip.rows();
    for (std::size_t r = 0; r < rows.num(); ++r)
    {
      SCIP_CONS* cons = NULL;
      SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &cons, mip.rowName(r).c_str(), 0, 0, 0,
        double(rows.lhs(r)), double(rows.rhs(r))));
      const SVectorRational& row = rows.rowVector(r);
      for (int p = row.size() - 1; p >= 0; --p)
      {
        SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _variables[row.index(p)],
double(row.value(p))));
      }
      SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
      SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
    }
  }

  Rational SCIPOracle::computeVectorScaling(const VectorRational& vector)
  {
    assert(vector.dim() == space().dimension());
    Rational largest = 0;
    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      if (vector[v] > 0 && vector[v] > largest)
        largest = vector[v];
      else if (vector[v] < 0 && vector[v] < -largest)
        largest = -vector[v];
    }
    if (largest > 0 && (largest > _maxLargestCoefficient
      || largest < _minLargestCoefficient))
    {
      return _bestLargestCoefficient / largest;
    }
    else
      return 1;
  }

  void SCIPOracle::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;

    if (currentFace() != NULL)
    {
      SCIP_CALL_EXC(SCIPdelCons(_scip, _faceConstraint));
    }

    OracleBase::setFace(newFace);

    if (currentFace() != NULL)
    {
      DSVectorRational normalCopy;
      const Rational& largest = currentFace()->maxNorm();
      Rational scaling = 1;
      if (largest != 0 && (largest < _minLargestCoefficient || largest > _maxLargestCoefficient))
      {
        scaling = _bestLargestCoefficient / largest;
        normalCopy = currentFace()->sparseNormal() * scaling;
      }
      const SVectorRational& normal = (scaling == 1) ? currentFace()->sparseNormal() : normalCopy;
      Rational rhs = currentFace()->rhs() * scaling;

      /// Add constraint to SCIP.

      SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &_faceConstraint, "face", 0, 0, 0,
        double(rhs),   double(rhs)));
      for (int p = normal.size() - 1; p >= 0; --p)
      {
        SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, _faceConstraint, _variables[normal.index(p)],
          double(normal.value(p))));
      }

      SCIP_CALL_EXC(SCIPaddCons(_scip, _faceConstraint));
    }
  }

  void SCIPOracle::maximize(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic)
  {
    std::size_t n = space().dimension();
    if (objective.dim() != n)
      throw std::runtime_error("Oracle called with objective vector of wrong dimension!");

    /// Set SCIP objective to scaled objective vector.

    Rational scalingFactor = computeVectorScaling(objective);
    for (std::size_t v = 0; v < n; ++v)
    {
      SCIP_CALL_EXC(SCIPchgVarObj(_scip, _variables[v], double(scalingFactor * objective[v])));
    }

    result.buildStart(objective);
    int oldMaxRounds;
    SCIP_CALL_EXC(SCIPgetIntParam(_scip, "presolving/maxrounds", &oldMaxRounds));
    for (int attempt = 1; attempt <= 2; ++attempt)
    {
      try
      {
        SCIP_CALL_EXC(SCIPsolve(_scip));
      }
      catch (SCIPException& ex)
      {
        std::cerr << "SCIPOracle failed. Writing problem to oracle-solve-failed.lp!" << std::endl;
        SCIP_CALL_EXC(SCIPwriteOrigProblem(_scip, "oracle-solve-failed.lp", NULL, false));
        throw ex;
      }

      bool hasRay = SCIPhasPrimalRay(_scip);
      if (hasRay)
      {
        DSVectorRational* direction = new DSVectorRational;
        for (std::size_t v = 0; v < n; ++v)
        {
          double apxCoeff = SCIPgetPrimalRayVal(_scip, _variables[v]);
          if (SCIPisZero(_scip, apxCoeff))
            continue;
          Rational exCoeff;
          reconstruct(apxCoeff, exCoeff, 10 * SCIPfeastol(_scip));
          if (exCoeff != 0)
            direction->add(v, exCoeff);
        }
        result.buildAddDirection(direction);
        break;
      }

      SCIP_STATUS status = SCIPgetStatus(_scip);
      std::size_t numSolutions = SCIPgetNSols(_scip);
      if (status != SCIP_STATUS_INFEASIBLE && status != SCIP_STATUS_INFORUNBD && numSolutions > 0)
      {
        SCIP_SOL** solutions = SCIPgetSols(_scip);
        for (std::size_t solIndex = 0; solIndex < numSolutions; ++solIndex)
        {
          SCIP_SOL* sol = solutions[solIndex];
          DSVectorRational* point = new DSVectorRational;
          for (std::size_t v = 0; v < n; ++v)
          {
            double apxValue = SCIPgetSolVal(_scip, sol, _variables[v]);
            Rational exValue;
            if (SCIPvarIsIntegral(_variables[v]))
              exValue = int(apxValue + 0.5);
            else
              reconstruct(apxValue, exValue, SCIPfeastol(_scip));
            if (exValue != 0)
              point->add(v, exValue);
          }
          result.buildAddPoint(point);
        }

        /// TODO: Use exact primal SCIP functionality instead of reconstruction.

        break;
      }

      // Disable presolving for the second round.

      SCIP_CALL_EXC(SCIPsetIntParam(_scip, "presolving/maxrounds", 0));
      SCIP_CALL_EXC(SCIPfreeSolve(_scip, true));
      SCIP_CALL_EXC(SCIPfreeTransform(_scip));
    }

    SCIP_CALL_EXC(SCIPsetIntParam(_scip, "presolving/maxrounds", oldMaxRounds));
    SCIP_CALL_EXC(SCIPfreeSolve(_scip, true));

    // TODO: Currently needed for SCIPchgVarObj, but isn't there a warm-start for that?
    SCIP_CALL_EXC(SCIPfreeTransform(_scip));

    result.buildFinish(thisHeuristic(), true, true, true);
  }

//   ExactSCIPOptimizationOracle::ExactSCIPOptimizationOracle(const std::string& name, const std::string& exactBinary,
//       MixedIntegerProgram& mip, FaceOptimizationOracleBase* heuristic, double timeLimit) :
//       FaceOptimizationOracleBase(name, mip.space()), scaleObjective(true), _binary(exactBinary),
// _mip(mip), _heuristic(heuristic), _timeLimit(
//           timeLimit)
//   {
//     if (_heuristic != NULL)
//     {
//       if (_heuristic->space() != space())
//         throw std::runtime_error("Spaces of MixedIntegerProgram and heuristic differ.");
//     }
//
//     createTempDirectory();
//   }
//
//   ExactSCIPOptimizationOracle::~ExactSCIPOptimizationOracle()
//   {
//     deleteTempDirectory();
//   }
//
//   void ExactSCIPOptimizationOracle::run(OptimizationResult& result, const VectorRational& objective,
//       const Rational* improveValue, bool forceOptimal)
//   {
//     /// Run approximate oracle first.
//
//     Rational primalObjective = improveValue != NULL ? *improveValue : Rational(-infinity);
//     if (_heuristic && !forceOptimal)
//     {
// //       TODO:
// //       if (improveValue != NULL)
// //         _heuristic->improve(result, objective, *improveValue, forceOptimal);
// //       else
//         _heuristic->maximize(result, objective, forceOptimal);
//
//       if (result.isUnbounded() || result.isInfeasible())
//         return;
//
//       if (result.isFeasible())
//       {
//         if (!forceOptimal || result.bestValue > *improveValue)
//           return;
//         primalObjective = result.bestValue;
//       }
//     }
//
//     /// TODO: Change model such that only solutions better than primalObjective are searched for.
//
//     /// Solve model.
//
//     writeModel(objective);
//     callSolver();
//     DSVectorRational* optimum = parseOutput();
//     result.optimal = true;
//     if (optimum != NULL)
//     {
//       Rational optimalValue = *optimum * objective;
//       if (optimalValue > result.bestValue)
//       {
//         result.points.push_back(optimum);
//         result.objectives.push_back(optimalValue);
//         result.bestIndex = 0;
//         result.bestValue = optimalValue;
//       }
//     }
//     else if (optimum == NULL && result.isInfeasible())
//     {
//       result.optimal = true;
//     }
//     else
//     {
//       throw std::runtime_error("Feasibility stati of oracles differ!");
//     }
//   }
//
//   void ExactSCIPOptimizationOracle::faceEnabled(Face* face)
//   {
//     _mip.faceEnabled(face);
//     if (_heuristic != NULL)
//       _heuristic->setFace(face);
//   }
//
//   void ExactSCIPOptimizationOracle::faceDisabled(Face* face)
//   {
//     if (_heuristic != NULL)
//       _heuristic->setFace(NULL);
//     _mip.faceDisabled(face);
//   }
//
//   void ExactSCIPOptimizationOracle::createTempDirectory()
//   {
//     char buffer[256] = "/tmp/ipo-ExactSCIPOptimizationOracle-XXXXXX";
//     char* name = mkdtemp(buffer);
//     if (name == NULL)
//       throw std::runtime_error("Cannot create temporary directory!");
//     _path = std::string(name);
//
//     std::string binary = "/usr/bin/time -o '" + _path + "/timing.log' -f '%U' " + _binary + " ";
//     std::string parameters = "";
//     if (_timeLimit < std::numeric_limits<double>::max())
//     {
//       std::stringstream ss;
//       ss << " -c \"set limits time ";
//       ss << _timeLimit;
//       ss << "\"";
//       parameters += ss.str();
//     }
//     std::string failureParameters = " -c \"set misc usefprelax FALSE\" -c \"set presolving maxrounds 0\" ";
//     std::string commands = "-c \"read model.zpl\" -c optimize -c \"display solution\" -c quit ";
//
//     std::ofstream file((_path + "/script.sh").c_str());
//     file << "#!/bin/bash\n\n";
//     file << "cd " << _path << "\n";
//     file << binary << " " << parameters << " " << commands << " > solve.log 2>&1\n";
//     file << "retcode=$?\n";
//     file << "if [[ $retcode != 0 ]]; then\n";
//     file << "  " << binary << parameters << failureParameters << commands << " > solve.log 2>&1\n";
//     file << "fi\n";
//     file.close();
//     chmod((_path + "/script.sh").c_str(), 00700);
//   }
//
//   void ExactSCIPOptimizationOracle::deleteTempDirectory()
//   {
//     unlink((_path + "/solve.log").c_str());
//     unlink((_path + "/script.sh").c_str());
//     unlink((_path + "/model.zpl").c_str());
//     rmdir(_path.c_str());
//   }
//
//   void ExactSCIPOptimizationOracle::writeModel(const VectorRational& objective)
//   {
//     std::ofstream file((_path + "/model.zpl").c_str());
//
//     const LPColSetRational& columns = _mip.columns();
//     const LPRowSetRational& rows = _mip.rows();
//     for (std::size_t v = 0; v < space().dimension(); ++v)
//     {
//       file << "var x" << v;
//       if (_mip.isIntegral(v))
//         file << " integer";
//       else
//         file << " real";
//       if (columns.lower(v) > -infinity)
//         file << " >= " << columns.lower(v);
//       if (columns.upper(v) < infinity)
//         file << " <= " << columns.upper(v);
//       file << ";\n";
//     }
//     file << "\nmaximize cost:";
//     bool first = true;
//     for (std::size_t v = 0; v < space().dimension(); ++v)
//     {
//       if (objective[v] == 0)
//         continue;
//       file << "\n";
//       if (first)
//         first = false;
//       else
//         file << " + ";
//       file << objective[v] << "*x" << v;
//     }
//     if (first)
//       file << "0*x0";
//     file << ";\n\n";
//     for (int r = 0; r < rows.num(); ++r)
//     {
//       const SVectorRational& vector = rows.rowVector(r);
//       if (rows.lhs(r) == rows.rhs(r))
//       {
//         file << "\nsubto row" << r << ":";
//         first = true;
//         for (int p = vector.size() - 1; p >= 0; --p)
//         {
//           file << "\n";
//           if (first)
//             first = false;
//           else
//             file << " + ";
//           file << vector.value(p) << "*x" << vector.index(p);
//         }
//         if (first)
//           file << "x0 - x0";
//         file << " == " << rows.lhs(r) << ";\n\n";
//         continue;
//       }
//       if (rows.lhs(r) > -infinity)
//       {
//         file << "\nsubto row" << r << "lhs:";
//         first = true;
//         for (int p = vector.size() - 1; p >= 0; --p)
//         {
//           file << "\n";
//           if (first)
//             first = false;
//           else
//             file << " + ";
//           file << vector.value(p) << "*x" << vector.index(p);
//         }
//         if (first)
//           file << "x0 - x0";
//         file << " >= " << rows.lhs(r) << ";\n\n";
//       }
//       if (rows.rhs(r) < infinity)
//       {
//         file << "\nsubto row" << r << "rhs:";
//         first = true;
//         for (int p = vector.size() - 1; p >= 0; --p)
//         {
//           file << "\n";
//           if (first)
//             first = false;
//           else
//             file << " + ";
//           file << vector.value(p) << "*x" << vector.index(p);
//         }
//
//         /// Handle trivially violated constraints explicitly.
//
//         if (first)
//           file << "x0 - x0";
//         file << " <= " << rows.rhs(r) << ";\n\n";
//       }
//     }
//     file.close();
//   }
//
//   void ExactSCIPOptimizationOracle::callSolver()
//   {
//     if (system((_path + "/script.sh").c_str()) != 0)
//     {
//       throw std::runtime_error("Exact SCIP solver returned with nonzero exit status.");
//     }
//   }
//
//   DSVectorRational* ExactSCIPOptimizationOracle::parseOutput()
//   {
//     std::ifstream log((_path + "/solve.log").c_str());
//     std::string line;
//     bool startedSolutionSection = false;
//     bool timeLimitReached = false;
//     DSVectorRational* result = NULL;
//     while (std::getline(log, line))
//     {
//       if (line.substr(0, 16) == "objective value:")
//       {
//         startedSolutionSection = true;
//         result = new DSVectorRational;
//       }
//
//       if (startedSolutionSection && !line.empty() && line[0] == 'x')
//       {
//         std::size_t var;
//         std::string valueStr;
//         std::stringstream ss(line.substr(1, std::string::npos));
//         ss >> var >> valueStr;
//         Rational value;
//         if (!value.readString(valueStr.c_str()))
//           throw std::runtime_error("parseOutput failed when reading a number.");
//         assert(result != NULL);
//         result->add(var, value);
//       }
//
//       if (line == "no solution available")
//       {
//         return NULL;
//       }
//       if (line == "SCIP Status        : solving was interrupted [time limit reached]")
//       {
//         timeLimitReached = true;
//         throw std::runtime_error("Time limit for ExactSCIPOptimizationOracle reached.");
//       }
//     }
//
//     std::ifstream timing((_path + "/timing.log").c_str());
//     double time;
//     timing >> time;
//     addTimeToActiveTimers(time);
//
//     if (!startedSolutionSection)
//       throw std::runtime_error(
//           "ExactSCIPOptimizationOracle did not return useful results (see " + _path + "/solve.log)");
//
//     return result;
//   }

} /* namespace ipo */
