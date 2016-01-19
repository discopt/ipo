#include "common.h"

#ifdef WITH_SCIP
#include <scip/scipdefplugins.h>
#include "ipo/scip_exception.hpp"
#include "ipo/scip_oracles.h"
#endif

#include "ipo/wrapper_oracle.h"

using namespace ipo;

MIPOracleInfo::MIPOracleInfo() :
    mip(NULL), corrector(NULL), inexact(NULL), deleteCorrector(false)
{

}

MIPOracleInfo::~MIPOracleInfo()
{
  if (deleteCorrector)
    delete corrector;
  if (inexact != NULL)
    delete inexact;
  if (mip)
    delete mip;
}

void printOracleArguments(std::ostream& stream)
{
  stream << "The ORACLE arguments specify the type and instance of the oracle to be used:\n";
  stream << " --soplex MIP             Use SoPlex for the LP relaxation of MIP.\n";
  stream << " --scip MIP               Use SCIP for the MIP.\n";
  stream << " --scip+ MIP              Use SCIP for the MIP followed by a corrector for continuous variables.\n";
  stream << " --scipex=BINARY MIP      Use exact SCIP (binary specified) for the MIP.\n";
  stream
      << " --scip+scipex=BINARY MIP Use exact SCIP (binary specified) and regular SCIP as a heuristic for the MIP.\n";
  stream << " WRAPPER INSTANCE...      Call the oracle WRAPPER and pass INSTANCE parameters.\n";
}

#ifdef WITH_SCIP

bool readSCIP(const std::string& name, SCIP*& scip, bool relaxation = false)
{
  scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  SCIP_CALL_EXC(SCIPreadProb(scip, name.c_str(), NULL));

  if (relaxation)
  {
    int nvars = SCIPgetNOrigVars(scip);
    bool changed = true;
    while (changed)
    {
      changed = false;
      SCIP_VAR** vars = SCIPgetOrigVars(scip);
      unsigned int infeasible;
      for (int i = 0; i < nvars; ++i)
      {
        if (SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS)
        {
          SCIP_CALL_EXC(SCIPchgVarType(scip, vars[i], SCIP_VARTYPE_CONTINUOUS, &infeasible));
          changed = true;
        }
      }
    }
  }

  SCIP_CALL_EXC(SCIPtransformProb(scip));
  return true;
}

bool readMIP(const std::string& name, MixedIntegerProgram*& mip, bool relaxation = false)
{
  SCIP* scip;
  readSCIP(name, scip, relaxation);
  mip = new MixedIntegerProgram(scip);
  SCIP_CALL_EXC(SCIPfree(&scip));

  return true;
}
#endif

bool parseOracleArguments(int numArgs, char** args, OptimizationOracleBase*& oracle, MIPOracleInfo& info)
{
#ifdef WITH_SCIP

  if (numArgs < 1)
    return false;
  if (args[0] == std::string("--soplex"))
  {
    if (!readMIP(args[1], info.mip, true))
      return false;
    info.inexact = new SCIPOptimizationOracle(args[1], *info.mip);
    info.corrector = new MixedIntegerProgramCorrectorOracle("corrected-" + info.inexact->name(), *info.mip,
        info.inexact);
    oracle = info.corrector;
    return true;
  }
  if (args[0] == std::string("--scip"))
  {
    SCIP* scip;
    readSCIP(args[1], scip);
    oracle = new SCIPOptimizationOracle(args[1], scip);
    SCIP_CALL_EXC(SCIPfree(&scip));

    return true;
  }
  else if (args[0] == std::string("--scip+"))
  {
    if (!readMIP(args[1], info.mip))
      return false;
    info.inexact = new SCIPOptimizationOracle(args[1], *info.mip);
    info.corrector = new MixedIntegerProgramCorrectorOracle("corrected-" + info.inexact->name(), *info.mip,
        info.inexact);
    oracle = info.corrector;
    return true;
  }
  else if (std::string(args[0]).find("--scipex=") == 0)
  {
    if (!readMIP(args[1], info.mip))
      return false;
    std::string binary = std::string(args[0]).substr(strlen("--scipex="));
    oracle = new ExactSCIPOptimizationOracle("exact-" + std::string(args[1]), binary, *info.mip, NULL, 3600.0);
    return true;
  }
  else if (std::string(args[0]).find("--scip+scipex=") == 0)
  {
    if (!readMIP(args[1], info.mip))
      return false;
    info.inexact = new SCIPOptimizationOracle(args[1], *info.mip);
    info.corrector = new MixedIntegerProgramCorrectorOracle("corrected-" + info.inexact->name(), *info.mip,
        info.inexact);
    info.deleteCorrector = true;
    std::string binary = std::string(args[0]).substr(strlen("--scip+scipex="));
    oracle = new ExactSCIPOptimizationOracle("exact-" + info.inexact->name(), binary, *info.mip, info.corrector,
        3600.0);
    return true;
  }
  else
  {
    std::stringstream stream;
    for (int a = 1; a < numArgs; ++a)
      stream << args[a];
    oracle = new WrapperOptimizationOracle(args[0], args[0], stream.str());
    return true;
  }
#endif
}

