#ifndef IPO_PYTHON_WRAPPER_H_
#define IPO_PYTHON_WRAPPER_H_

#include "common.h"

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
#endif

#include "scip_oracle.h"
#include "exactscip_oracle.h"
#include "scip_exception.h"
#include "affine_hull.h"
#include "cache_oracle.h"
#include "statistics_oracle.h"
#include "vectors.h"
#include "mip.h"
#include "polyhedron.h"

namespace ipo {

const std::shared_ptr<MixedIntegerSet> getMixedIntegerSet(std::string filename);

class OracleControllerBase
  {
    public:
    OracleControllerBase();
    //Destructor
    virtual ~OracleControllerBase();

    /*
     * Methods which can be used by the oracles
     * name(): returns the name of the SCIPOracle in this bundle
     *
     * heuristicLevel_ScipOracle(): 
     * \brief returns the heuristiclevel of the SCIPOracle in this bundle
     *
     */
    virtual std::string name() = 0;

    /* heuristicLevel_CacheOracle(): 
     * \brief returns the heuristiclevel of the CacheOracle in this bundle
     */
    virtual int heuristicLevel_ScipOracle() = 0;


    /* heuristicLevel_CacheOracle(): 
     * \brief returns the heuristiclevel of the CacheOracle in this bundle
     */
    virtual int heuristicLevel_CacheOracle() = 0;

    /* space(): 
     * \brief returns the Space that comes from the CacheOracle->space() call
     */
    virtual Space space() = 0;


    /* getConnectionOracle():
     * \brief returns the Oracle which can be used in the creation of another to create a link (Constructor parameter)
     */
    virtual std::shared_ptr<StatisticsOracle> getConnectionOracle() = 0;
  
      /* getCacheOracle():
     * \brief returns the CacheOracle
     */
    virtual std::shared_ptr<OracleBase> getCacheOracle() = 0;

    private:

  };

  class ScipOracleController: public OracleControllerBase
  {
    public:
    //Constructors
    //creates ScipOracle with Statistics Wrapper and Cache
    ScipOracleController(std::string name);
    ScipOracleController(std::string name, const std::shared_ptr<MixedIntegerSet> mixedIntegerSet);

    ScipOracleController(std::string name, OracleControllerBase* next);
    ScipOracleController(std::string name, const std::shared_ptr<MixedIntegerSet> mixedIntegerSet, OracleControllerBase* next);

    virtual ~ScipOracleController();

    /*
     * Methods which can be used by the oracles
     * name(): returns the name of the SCIPOracle in this bundle
     *
     * heuristicLevel_ScipOracle(): 
     * \brief returns the heuristiclevel of the SCIPOracle in this bundle
     *
     */
    std::string name();

    /* heuristicLevel_CacheOracle(): 
     * \brief returns the heuristiclevel of the CacheOracle in this bundle
     */
    int heuristicLevel_ScipOracle();


    /* heuristicLevel_CacheOracle(): 
     * \brief returns the heuristiclevel of the CacheOracle in this bundle
     */
    int heuristicLevel_CacheOracle();

    /* space(): 
     * \brief returns the Space that comes from the CacheOracle->space() call
     */
    Space space();

    /* affineHullInner(): 
     * \brief creates handlers and returns the InnerDescription
     */
    InnerDescription affineHullInner(int outputMode);

    /* affineHullOuter(): 
     * \brief creates handlers and returns the AffineOuterDescription
     */
    AffineOuterDescription affineHullOuter(int outputMode);

    /* getConnectionOracle():
     * \brief returns the Oracle which can be used in the creation of another to create a link (Constructor parameter)
     */
    std::shared_ptr<StatisticsOracle> getConnectionOracle();
  
      /* getCacheOracle():
     * \brief returns the CacheOracle
     */
    std::shared_ptr<OracleBase> getCacheOracle();

    private:
    std::shared_ptr<StatisticsOracle> scipOracle;
    std::shared_ptr<SCIPOracle> scipOracleImpl;
    std::shared_ptr<CacheOracle> cacheOracleImpl;
    std::shared_ptr<StatisticsOracle> cacheOracle;
    std::shared_ptr<InnerDescription> cacheInner;
    std::shared_ptr<AffineOuterDescription> cacheOuter;
  };


#ifdef IPO_WITH_EXACT_SCIP
  class ExactScipOracleController: public OracleControllerBase
  {
    public:
    //Constructors
    //creates ScipOracle with Statistics Wrapper and Cache
    ExactScipOracleController(std::string name, const std::shared_ptr<MixedIntegerSet> mixedIntegerSet);

    ExactScipOracleController(std::string name, const std::shared_ptr<MixedIntegerSet> mixedIntegerSet, OracleControllerBase* next);

    virtual ~ExactScipOracleController();

    /*
     * Methods which can be used by the oracles
     * name(): returns the name of the SCIPOracle in this bundle
     *
     * heuristicLevel_ScipOracle(): 
     * \brief returns the heuristiclevel of the SCIPOracle in this bundle
     *
     */
    std::string name();

    /* heuristicLevel_CacheOracle(): 
     * \brief returns the heuristiclevel of the CacheOracle in this bundle
     */
    int heuristicLevel_ExactScipOracle();


    /* heuristicLevel_CacheOracle(): 
     * \brief returns the heuristiclevel of the CacheOracle in this bundle
     */
    int heuristicLevel_CacheOracle();

    /* space(): 
     * \brief returns the Space that comes from the CacheOracle->space() call
     */
    Space space();

    /* affineHullInner(): 
     * \brief creates handlers and returns the InnerDescription
     */
    InnerDescription affineHullInner(int outputMode);

    /* affineHullOuter(): 
     * \brief creates handlers and returns the AffineOuterDescription
     */
    AffineOuterDescription affineHullOuter(int outputMode);

    /* getConnectionOracle():
     * \brief returns the Oracle which can be used in the creation of another to create a link (Constructor parameter)
     */
    std::shared_ptr<StatisticsOracle> getConnectionOracle();
  
      /* getCacheOracle():
     * \brief returns the CacheOracle
     */
    std::shared_ptr<OracleBase> getCacheOracle();

    private:
    std::shared_ptr<StatisticsOracle> scipOracle;
    std::shared_ptr<ExactSCIPOracle> scipOracleImpl;
    std::shared_ptr<CacheOracle> cacheOracleImpl;
    std::shared_ptr<StatisticsOracle> cacheOracle;
    std::shared_ptr<InnerDescription> cacheInner;
    std::shared_ptr<AffineOuterDescription> cacheOuter;
  };

#endif
}


#endif /* IPO_PYTHON_WRAPPER_H_ */
