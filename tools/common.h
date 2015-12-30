#ifndef COMMON_H_
#define COMMON_H_

#include "ipo/oracles.h"
#include "ipo/mixed_integer_program.h"

struct MIPOracleInfo
{
  ipo::MixedIntegerProgram* mip;
  ipo::MixedIntegerProgramCorrectorOracle* corrector;
  ipo::FaceOptimizationOracleBase* inexact;
  bool deleteCorrector;

  MIPOracleInfo();
  virtual ~MIPOracleInfo();
};

void printOracleArguments(std::ostream& stream);
bool parseOracleArguments(int numArgs, char** args, ipo::OptimizationOracleBase*& oracle,
    MIPOracleInfo& info);

#endif /* COMMON_H_ */
