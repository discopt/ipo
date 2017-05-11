#include <iostream>

#include "ipo/scip_exception.h"
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

int printUsage(const std::string& program)
{
  std::cerr << "Usage: " << program << " INPUT OUTPUT\n" << std::flush;
  return EXIT_FAILURE;
}

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    return printUsage(argv[0]);
  }

  std::cout << "Writing LP relaxation of " << argv[1] << " to " << argv[2] << std::endl;

  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPreadProb(scip, argv[1], NULL));


  SCIP_VAR** vars = SCIPgetOrigVars(scip);
  int nvars = SCIPgetNOrigVars(scip);
  for (int i = nvars - 1; i >= 0; --i)
  {
    SCIP_Bool infeasible;
    SCIP_CALL_EXC(SCIPchgVarType(scip, vars[i], SCIP_VARTYPE_CONTINUOUS, &infeasible));
  }
  
  SCIP_CALL_EXC(SCIPwriteOrigProblem(scip, argv[2], NULL, false));

  SCIP_CALL_EXC(SCIPfree(&scip));

  return 0;
}
