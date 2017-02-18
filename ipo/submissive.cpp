#include "submissive.h"

using namespace soplex;

namespace ipo {

  SubmissiveOracle::SubmissiveOracle(const std::string& name, const std::shared_ptr<OracleBase>& sourceOracle,
    const std::shared_ptr<OracleBase>& nextOracle)
    : OracleBase(name, nextOracle), _sourceOracle(sourceOracle)
  {
    OracleBase::initializeSpace(sourceOracle->space());
    _objective.reDim(space().dimension());
  }

  SubmissiveOracle::~SubmissiveOracle()
  {

  }

  HeuristicLevel SubmissiveOracle::maximizeImplementation(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
  {
    for (std::size_t v = 0; v < space().dimension(); ++v)
      _objective[v] = (objective[v] >= 0) ? objective[v] : Rational(0);

    OracleResult sourceResult;
    _sourceOracle->maximize(sourceResult, _objective, objectiveBound, minHeuristic, maxHeuristic);
    if (sourceResult.isFeasible())
    {
      for (std::size_t i = 0; i < sourceResult.points.size(); ++i)
      {
        const Vector& sourcePoint = sourceResult.points[i].vector;
        VectorData* data = new VectorData(sourcePoint.size());
        for (std::size_t p = 0; p < sourcePoint.size(); ++p)
        {
          std::size_t index = sourcePoint.index(p);
          if (_objective[index] > 0)
            data->add(index, sourcePoint.value(p));
        }
        Vector vector(data);
        result.points.push_back(OracleResult::Point(vector));
      }
      sort = true;
      checkDups = true;
    }
    else if (sourceResult.isUnbounded())
      result.rays = sourceResult.rays;
  }

}
