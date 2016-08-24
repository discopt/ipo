#include "statistics_oracle.h"

namespace ipo {

  StatisticsOracle::ForwardWrapperOracle::ForwardWrapperOracle(StatisticsOracle& host,
    const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase("ForwardWrapper(" + nextOracle->name() + ")", nextOracle), _host(host)
  {
    initializeSpace(nextOracle->space());
  }

  StatisticsOracle::ForwardWrapperOracle::~ForwardWrapperOracle()
  {

  }

  HeuristicLevel StatisticsOracle::ForwardWrapperOracle::maximizeController(OracleResult& result,
    const soplex::VectorRational& objective, const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic,
    HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    assert((heuristicLevel() == 0 && _nextOracle == NULL)
      || heuristicLevel() > 0 && _nextOracle != NULL);

    _host.beforeForward();

    HeuristicLevel resultLevel = _nextOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, 
      sort, checkDups);

    _host.afterForward();
    
    return resultLevel;
  }

  HeuristicLevel StatisticsOracle::ForwardWrapperOracle::maximizeImplementation(OracleResult& result,
    const soplex::VectorRational& objective, const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic,
    HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    assert(false);
  }
  
  StatisticsOracle::StatisticsOracle(const std::shared_ptr<OracleBase>& targetOracle)
    : OracleBase(targetOracle->name(), targetOracle->nextOracle()), _targetOracle(targetOracle)
  {
    if (_targetOracle->nextOracle() != NULL)
    {
      targetOracle->_nextOracle = std::make_shared<ForwardWrapperOracle>(*this, _targetOracle->_nextOracle);
    }

    initializeSpace(targetOracle->space());

    _numCalls = 0;
    _numForwards = 0;
  }

  StatisticsOracle::~StatisticsOracle()
  {

  }
  
  void StatisticsOracle::setFace(const LinearConstraint& newFace)
  {
    _targetOracle->setFace(newFace);

    // We do not call OracleBase::setFace() since the above call will call setFace() for the ForwardWrapperOracle, which does it.
  }

  HeuristicLevel StatisticsOracle::maximizeController(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    assert((heuristicLevel() == 0 && _nextOracle == NULL)
      || heuristicLevel() > 0 && _nextOracle != NULL);
    
    // If requested, forward directly to next oracle, and do not count this as a call.

    if (heuristicLevel() > maxHeuristic)
    {
      return _nextOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, sort, checkDups);
    }

    // The target oracle is supposed to be called.

    _numCalls++;
    _timer.start();

    // We call the target oracle.

    HeuristicLevel resultLevel = _targetOracle->maximizeController(result, objective, objectiveBound, maxHeuristic, minHeuristic, 
      sort, checkDups);

    // The target oracle returned.

    _timer.stop();

    return resultLevel;
  }

  void StatisticsOracle::beforeForward()
  {
    // The target oracle will forward to its next associated oracle.

    _timer.stop();
    _numForwards++;
  }

  void StatisticsOracle::afterForward()
  {
    // The target oracle did forward to its next associated oracle. 

    _timer.start();
  }

  HeuristicLevel StatisticsOracle::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    assert(false);
  }

}