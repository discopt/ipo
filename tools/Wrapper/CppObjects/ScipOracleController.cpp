#include "ScipOracleController.h"

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <ipo/scip_exception.hpp>
  #include <ipo/scip_oracles.h>
  #include <scip/debug.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <ipo/scip_exception.hpp>
  #include <ipo/scip_oracles.h>
  #include <scip/debug.h>
#endif

ScipOracleController::ScipOracleController(std::string filename){
  //Create oracle triple from ScipOracle, StatisticsOracle, CacheOracle
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  //const char* filename
  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));
  scipOracleImpl = std::make_shared<SCIPOracle>("SCIPOracle("+filename+")", scip);
  SCIP_CALL_EXC(SCIPfree(&scip));
  scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  cacheOracleImpl = std::make_shared<CacheOracle>(scipOracle);
  cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);
}

ScipOracleController::ScipOracleController(std::string filename, ScipOracleController prev){
  //Create oracle triple from ScipOracle, StatisticsOracle, CacheOracle
  SCIP* scip = NULL;
  SCIP_CALL_EXC(SCIPcreate(&scip));
  SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
  SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
  //const char* filename
  const char *c_filename = filename.c_str();
  SCIP_CALL_EXC(SCIPreadProb(scip, c_filename, NULL));
  SCIP_CALL_EXC(SCIPtransformProb(scip));
  scipOracleImpl = std::make_shared<SCIPOracle>("SCIPOracle("+filename+")", scip, prev.getConnectionOracle());
  SCIP_CALL_EXC(SCIPfree(&scip));
  scipOracle = std::make_shared<StatisticsOracle>(scipOracleImpl);

  cacheOracleImpl = std::make_shared<CacheOracle>(scipOracle);
  cacheOracle = std::make_shared<StatisticsOracle>(cacheOracleImpl);
}

std::shared_ptr<StatisticsOracle> ScipOracleController::getConnectionOracle(){
  return this->scipOracle;
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

  std::vector<AffineHullHandler*> handlers; // Soll pro Aufruf lokal erzeugt werden.
  DebugAffineHullHandler debugHandler(std::cout);

  StatisticsAffineHullHandler statsHandler;
  handlers.push_back(&debugHandler);
  handlers.push_back(&statsHandler);

  AffineOuterDescription outer;

  affineHull(cacheOracle, inner, outer, handlers, 2, 1);

  return inner;
}

AffineOuterDescription ScipOracleController::affineHullOuter(int outputMode){
  InnerDescription inner;

  std::vector<AffineHullHandler*> handlers; // Soll pro Aufruf lokal erzeugt werden.
    
  DebugAffineHullHandler debugHandler(std::cout);
  StatisticsAffineHullHandler statsHandler;
  handlers.push_back(&debugHandler);
  handlers.push_back(&statsHandler);

  AffineOuterDescription outer;

  affineHull(cacheOracle, inner, outer, handlers, 2, 1);

  return outer;
}
