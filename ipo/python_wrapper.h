#ifndef IPO_PYTHON_WRAPPER_H_
#define IPO_PYTHON_WRAPPER_H_

#include <ipo/common.h>

#include "ipo/scip_exception.hpp"
#include "ipo/affine_hull.h"
#include "ipo/scip_oracles.h"
#include "ipo/cache_oracle.h"
#include "ipo/statistics_oracle.h"
#include <ipo/vectors.h>

namespace ipo {

  class ScipOracleController
  {
    public:
    //Constructors
    //creates ScipOracle with Statistics Wrapper and Cache
    ScipOracleController(std::string name);
    /*
     * \brief 
     */
    ScipOracleController(std::string name, ScipOracleController prev);

    ~ScipOracleController();

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
  
  

    private:
    std::shared_ptr<StatisticsOracle> scipOracle;
    std::shared_ptr<SCIPOracle> scipOracleImpl;
    std::shared_ptr<CacheOracle> cacheOracleImpl;
    std::shared_ptr<StatisticsOracle> cacheOracle;
};

}

#endif /* IPO_PYTHON_WRAPPER_H_ */