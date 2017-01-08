#include "python_wrapper.h"

namespace ipo {
//Base class

OracleControllerBase::OracleControllerBase(){
}

OracleControllerBase::~OracleControllerBase(){
}


const std::shared_ptr<MixedIntegerSet> getMixedIntegerSet(std::string filename){
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<MixedIntegerSet> mixedIntegerSet = std::make_shared<MixedIntegerSet>(scip);

  SCIP_CALL_EXC(SCIPfree(&scip));
return mixedIntegerSet;
}

ScipOracleController::ScipOracleController(std::string filename){
  //Create oracle triple from ScipOracle, StatisticsOracle, CacheOracle
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));

  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  //create SCIPOracle
  scipOracleImpl = std::make_shared<SCIPOracle>("SCIPOracle("+filename+")", scip);
  SCIP_CALL_EXC(SCIPfree(&scip));
  //create corresponding StatisticsOracle
  scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  //create CacheOracle with attached StatisticsOracle
  cacheOracleImpl = std::make_shared<CacheOracle>(scipOracle);
  cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);
}

ScipOracleController::ScipOracleController(std::string filename, OracleControllerBase* next){
  //Create oracle triple from ScipOracle, StatisticsOracle, CacheOracle
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));

  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  //Create SCIPOracle with previously existing Oracle Bundle
  scipOracleImpl = std::make_shared<SCIPOracle>("SCIPOracle("+filename+")", scip, next->getConnectionOracle());
  SCIP_CALL_EXC(SCIPfree(&scip));
  scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  //Create Cache Oracle
  cacheOracleImpl = std::make_shared<CacheOracle>(scipOracle);
  cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);
}

ScipOracleController::~ScipOracleController(){
}


std::shared_ptr<StatisticsOracle> ScipOracleController::getConnectionOracle(){
  return this->scipOracle;
}

std::shared_ptr<OracleBase> ScipOracleController::getCacheOracle(){
  std::shared_ptr<OracleBase> oracle = cacheOracle;
  return oracle;
}

int ScipOracleController::heuristicLevel_ScipOracle(){
  return scipOracle->heuristicLevel();
}

int ScipOracleController::heuristicLevel_CacheOracle(){
  return cacheOracle->heuristicLevel();
}

std::string ScipOracleController::name(){
  return scipOracle->name();
}

Space ScipOracleController::space(){
  return cacheOracle->space();
}

InnerDescription ScipOracleController::affineHullInner(int outputMode){
  InnerDescription inner;
  AffineOuterDescription outer;
  if(cacheInner == nullptr || cacheOuter == nullptr){

  //Create Handlers
  std::vector<AffineHullHandler*> handlers;
  DebugAffineHullHandler debugHandler(std::cout);

  StatisticsAffineHullHandler statsHandler;
  handlers.push_back(&debugHandler);
  handlers.push_back(&statsHandler);

  affineHull(cacheOracle, inner, outer, handlers, 2, 1);
  cacheInner = std::make_shared<InnerDescription>(inner);
  cacheOuter = std::make_shared<AffineOuterDescription>(outer);
  }
  else{
    inner = *cacheInner.get();
  }

  return inner;
}

AffineOuterDescription ScipOracleController::affineHullOuter(int outputMode){
  InnerDescription inner;
  AffineOuterDescription outer;
  if(cacheInner == nullptr || cacheOuter == nullptr){
  //Create Handlers
  std::vector<AffineHullHandler*> handlers;
  DebugAffineHullHandler debugHandler(std::cout);
  StatisticsAffineHullHandler statsHandler;
  handlers.push_back(&debugHandler);
  handlers.push_back(&statsHandler);

  affineHull(cacheOracle, inner, outer, handlers, 2, 1);
  cacheInner = std::make_shared<InnerDescription>(inner);
  cacheOuter = std::make_shared<AffineOuterDescription>(outer);

  }
  else{
    outer = *cacheOuter.get();
  }

  return outer;
}
}

#ifdef IPO_WITH_EXACT_SCIP
namespace ipo{
ExactScipOracleController::ExactScipOracleController(std::string filename, const std::shared_ptr<MixedIntegerSet> mixedIntegerSet){
  //Create oracle triple from ExactScipOracleController, StatisticsOracle, CacheOracle
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));

  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));

  std::shared_ptr<MixedIntegerSet> mixedIntegerSetlocal = std::make_shared<MixedIntegerSet>(scip);

  //create ExactScipOracleController
  const std::string name = "ExactSCIPOracle("+filename+")";
  scipOracleImpl = std::make_shared<ExactSCIPOracle>(name, mixedIntegerSetlocal);
  SCIP_CALL_EXC(SCIPfree(&scip));
  //create corresponding StatisticsOracle
  scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  //create CacheOracle with attached StatisticsOracle
  cacheOracleImpl = std::make_shared<CacheOracle>(scipOracle);
  cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);
}

ExactScipOracleController::ExactScipOracleController(std::string filename, const std::shared_ptr<MixedIntegerSet> mixedIntegerSet, OracleControllerBase* next){
  //Create oracle triple from ExactScipOracleController, StatisticsOracle, CacheOracle
  /*SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));

  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));*/

  //Create SCIPOracle with previously existing Oracle Bundle
  const std::string name = "ExactSCIPOracle("+filename+")";
  scipOracleImpl = std::make_shared<ExactSCIPOracle>(name, mixedIntegerSet, next->getConnectionOracle());
  //SCIP_CALL_EXC(SCIPfree(&scip));
  scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  //Create Cache Oracle
  cacheOracleImpl = std::make_shared<CacheOracle>(scipOracle);
  cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);
}

ExactScipOracleController::~ExactScipOracleController(){
}


std::shared_ptr<StatisticsOracle> ExactScipOracleController::getConnectionOracle(){
  return this->scipOracle;
}

std::shared_ptr<OracleBase> ExactScipOracleController::getCacheOracle(){
  std::shared_ptr<OracleBase> oracle = cacheOracle;
  return oracle;
}

int ExactScipOracleController::heuristicLevel_ExactScipOracle(){
  return scipOracle->heuristicLevel();
}

int ExactScipOracleController::heuristicLevel_CacheOracle(){
  return cacheOracle->heuristicLevel();
}

std::string ExactScipOracleController::name(){
  return scipOracle->name();
}

Space ExactScipOracleController::space(){
  return cacheOracle->space();
}

}
#endif
