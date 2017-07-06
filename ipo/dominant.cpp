#include "dominant.h"

using namespace soplex;

namespace ipo {

  DominantOracle::DominantOracle(const std::string& name, const std::shared_ptr<OracleBase>& sourceOracle,
    const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase(name, nextOracle), _sourceOracle(sourceOracle), _isFeasible(false)
  {
    OracleBase::initializeSpace(sourceOracle->space());
  }

  DominantOracle::~DominantOracle()
  {

  }

  HeuristicLevel DominantOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    if (!_isFeasible)
    {
      DVectorRational zero;
      zero.reDim(space().dimension(), true);
      OracleResult zeroResult;
      _sourceOracle->maximize(zeroResult, zero);
      if (zeroResult.isFeasible())
        _isFeasible = true;
      else
        return 0;
    }

    for (std::size_t v = 0; v < space().dimension(); ++v)
    {
      if (objective[v] > 0)
      {
        VectorData* data = new VectorData(1);
        data->add(v, 1);
        Vector vector(data);
        result.rays.push_back(OracleResult::Ray(vector));
        return 0;
      }
    }

    _sourceOracle->maximize(result, objective, objectiveBound, minHeuristic, maxHeuristic);
  }

}
