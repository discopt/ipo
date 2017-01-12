#include "scip_oracle.h"

#ifdef IPO_WITH_SCIP

#include <limits>
#include <unistd.h>
#include <sys/stat.h>

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #define NDEBUG
#else
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
#endif

#include "scip_exception.h"

#include "reconstruct.h"
#include "timer.h"

using namespace soplex;

namespace ipo {

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

  Vector getSCIPObjective(SCIP* scip, bool makeMaximization)
  {
    std::size_t n = SCIPgetNOrigVars(scip);
    VectorData* data = new VectorData(n);
    SCIP_VAR** vars = SCIPgetOrigVars(scip);
    for (std::size_t v = 0; v < n; ++v)
    {
      double value = SCIPvarGetObj(vars[v]);
      if (!SCIPisZero(scip, value))
      {
        Rational x;
        reconstruct(value, x, SCIPfeastol(scip));
        if (x != 0) // Nonzero-check above does not guarantee that reconstruct yields a nonzero.
          data->add(v, x);
      }
    }

    if (makeMaximization && SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE)
    {
      for (std::size_t p = 0; p < data->size(); ++p)
        data->value(p) *= -1;
    }

    return Vector(data);
  }

  SCIPOracle::SCIPOracle(const std::string& name, SCIP* originalSCIP, const std::shared_ptr<OracleBase>& nextOracle)
    : MIPOracleBase(name, nextOracle)
  {
    std::shared_ptr<MixedIntegerSet> mixedIntegerSet = constructFromSCIP(originalSCIP);

    initialize(mixedIntegerSet);
  }

  SCIPOracle::SCIPOracle(const std::string& name, const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet,
    const std::shared_ptr<OracleBase>& nextOracle)
    : MIPOracleBase(name, nextOracle)
  {
    constructFromMixedIntegerSet(mixedIntegerSet);

    initialize(mixedIntegerSet);
  }

  SCIPOracle::SCIPOracle(const std::string& fileName, const std::shared_ptr<OracleBase>& nextOracle)
    : MIPOracleBase(fileName, nextOracle)
  {
    std::shared_ptr<MixedIntegerSet> mixedIntegerSet = constructFromFile(fileName);

    initialize(mixedIntegerSet);
  }

  SCIPOracle::~SCIPOracle()
  {
    for (std::size_t v = 0; v < space().dimension(); ++v)
      SCIP_CALL_EXC(SCIPreleaseVar(_scip, &_variables[v]));
    SCIP_CALL_EXC(SCIPfree(&_scip));
  }

  std::shared_ptr<MixedIntegerSet> SCIPOracle::constructFromSCIP(SCIP* originalSCIP)
  {
    if (!SCIPisTransformed(originalSCIP))
    {
      throw std::runtime_error(
        "Initialization of SCIPOracle requires transformed problem!");
    }

    std::size_t n = SCIPgetNOrigVars(originalSCIP);
    SCIP_VAR** origVars = SCIPgetOrigVars(originalSCIP);

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
    SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0));

    _variables.resize(n);
    for (std::size_t v = 0; v < n; ++v)
    {
      SCIP_VAR* transVar = NULL;
      SCIP_CALL_EXC(SCIPgetTransformedVar(originalSCIP, origVars[v], &transVar));
      _variables[v] = static_cast<SCIP_VAR*>(SCIPhashmapGetImage(hashMap, transVar));
      SCIP_CALL_EXC(SCIPcaptureVar(_scip, _variables[v]));
    }
    SCIPhashmapFree(&hashMap);

    return std::make_shared<MixedIntegerSet>(originalSCIP);
  }

  void SCIPOracle::constructFromMixedIntegerSet(const std::shared_ptr<MixedIntegerSet>& mixedIntegerSet)
  {
    std::size_t n = mixedIntegerSet->space().dimension();

    // Initialize SCIP.

    SCIP_CALL_EXC(SCIPcreate(&_scip));
    SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
    SCIP_CALL_EXC(SCIPcreateProbBasic(_scip, _name.c_str()));
    SCIP_CALL_EXC(SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE));
    SCIP_CALL_EXC(SCIPsetBoolParam(_scip, "misc/catchctrlc", 0));
    SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0));

    // Create variables.

    _variables.resize(mixedIntegerSet->numVariables());
    for (std::size_t v = 0; v < mixedIntegerSet->numVariables(); ++v)
    {
      const MixedIntegerSet::Variable& variable = mixedIntegerSet->variable(v);
      SCIP_CALL_EXC(SCIPcreateVarBasic(_scip, &_variables[v], mixedIntegerSet->space()[v].c_str(), double(variable.lowerBound),
        double(variable.upperBound), 0.0, mixedIntegerSet->isIntegral(v) ? SCIP_VARTYPE_INTEGER : SCIP_VARTYPE_CONTINUOUS));
      SCIP_CALL_EXC(SCIPaddVar(_scip, _variables[v]));
    }

    // Create row constraints.

    const std::vector<LinearConstraint>& rows = mixedIntegerSet->rowConstraints();
    for (std::size_t r = 0; r < rows.size(); ++r)
    {
      const LinearConstraint& constraint = rows[r];
      SCIP_CONS* cons = NULL;
      double lhs, rhs;
      if (constraint.type() == '<')
      {
        lhs = -infinity;
        rhs = double(constraint.rhs());
      }
      else if (constraint.type() == '>')
      {
        rhs = infinity;
        lhs = double(constraint.rhs());
      }
      else
      {
        assert(constraint.type() == '=');
        lhs = rhs = constraint.rhs();
      }
      SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &cons, mixedIntegerSet->rowName(r).c_str(), 0, 0, 0, lhs, rhs));
      const Vector& normal = constraint.normal();
      for (std::size_t p = 0; p < normal.size(); ++p)
      {
        SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, cons, _variables[normal.index(p)], normal.approximation(p)));
      }
      SCIP_CALL_EXC(SCIPaddCons(_scip, cons));
      SCIP_CALL_EXC(SCIPreleaseCons(_scip, &cons));
    }
  }
  
  std::shared_ptr<MixedIntegerSet> SCIPOracle::constructFromFile(const std::string& fileName)
  {
    SCIP_CALL_EXC(SCIPcreate(&_scip));
    SCIP_CALL_EXC(SCIPincludeDefaultPlugins(_scip));
    SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0));
    SCIP_CALL_EXC(SCIPreadProb(_scip, fileName.c_str(), NULL));
    SCIP_CALL_EXC(SCIPtransformProb(_scip));

    std::size_t n = SCIPgetNOrigVars(_scip);
    SCIP_VAR** origVars = SCIPgetOrigVars(_scip);
    for (std::size_t v = 0; v < n; ++v)
      _variables[v] = origVars[v];

    return std::make_shared<MixedIntegerSet>(_scip);
  }

  void SCIPOracle::setFace(const LinearConstraint& newFace)
  {
    if (newFace == currentFace())
      return;

    if (!currentFace().definesCompleteFace())
    {
      SCIP_CALL_EXC(SCIPdelCons(_scip, _faceConstraint));
      SCIP_CALL_EXC(SCIPreleaseCons(_scip, &_faceConstraint));
    }

    OracleBase::setFace(newFace);

    if (!currentFace().definesCompleteFace())
    {
      const Vector& normal = currentFace().normal();
      double rhs = double(currentFace().rhs());

      // Add constraint to SCIP.

      assert(currentFace().isEquation());
      SCIP_CALL_EXC(SCIPcreateConsBasicLinear(_scip, &_faceConstraint, "face", 0, 0, 0, rhs, rhs));
      for (std::size_t p = 0; p < normal.size(); ++p)
      {
        SCIP_CALL_EXC(SCIPaddCoefLinear(_scip, _faceConstraint, _variables[normal.index(p)],
          normal.approximation(p)));
      }

      SCIP_CALL_EXC(SCIPaddCons(_scip, _faceConstraint));
    }
  }

  void SCIPOracle::solverMaximize(double* objective, double objectiveBound, std::vector<double*>& points,
    std::vector<double*>& rays)
  {
    std::cout << "DEBUG: Query to " << this->name() << std::endl;

    std::size_t n = space().dimension();

    for (std::size_t v = 0; v < n; ++v)
    {
      SCIP_CALL_EXC(SCIPchgVarObj(_scip, _variables[v], objective[v]));
    }

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
        rays.push_back(NULL);
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
          double* point = new double[n];
          for (std::size_t v = 0; v < n; ++v)
            point[v] = SCIPgetSolVal(_scip, sol, _variables[v]);
          points.push_back(point);
        }

        // TODO: Use exact primal SCIP functionality instead of reconstruction.

        break;
      }

      // Disable presolving for the second round.
      SCIP_CALL_EXC(SCIPsetIntParam(_scip, "presolving/maxrounds", 0));
      SCIP_CALL_EXC(SCIPfreeSolve(_scip, true));
      SCIP_CALL_EXC(SCIPfreeTransform(_scip));
    }

    SCIP_CALL_EXC(SCIPsetIntParam(_scip, "presolving/maxrounds", oldMaxRounds));
    SCIP_CALL_EXC(SCIPfreeSolve(_scip, true));
    SCIP_CALL_EXC(SCIPfreeTransform(_scip));
  }

} /* namespace ipo */

#endif /* IPO_WITH_SCIP */
