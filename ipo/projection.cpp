#include "projection.h"

using namespace soplex;

namespace ipo {

  Projection::Projection(const Space& sourceSpace): _sourceSpace(sourceSpace)
  {

  }

  Projection::Projection(const Space& sourceSpace, const std::vector< std::size_t >& variableSubset)
    : _sourceSpace(sourceSpace)
  {
    _variables.reserve(variableSubset.size());
    _map.reserve(variableSubset.size());
    _shift.reserve(variableSubset.size());
    for (std::size_t i = 0; i < variableSubset.size(); ++i)
      addVariable(variableSubset[i]);
  }

  void Projection::addVariable(const std::string& variableName,
    const SVectorRational& coefficients, const Rational& shift)
  {
    _variables.push_back(variableName);
    _map.push_back(DSVectorRational(coefficients));
    _shift.push_back(shift);
  }

  void Projection::addVariable(std::size_t sourceVariable, const soplex::Rational& shift)
  {
    _variables.push_back(_sourceSpace[sourceVariable]);
    _map.push_back(DSVectorRational(UnitVectorRational(sourceVariable)));
    _shift.push_back(shift);
  }

  void Projection::projectPoint(const VectorRational& point, DSVectorRational& projectedPoint) const
  {
    projectedPoint.clear();
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      Rational x = _shift[v] + _map[v] * point;
      if (x != 0)
        projectedPoint.add(v, x);
    }
  }

  void Projection::projectDirection(const VectorRational& direction,
    DSVectorRational& projectedDirection) const
  {
    projectedDirection.clear();
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      Rational x = _map[v] * direction;
      if (x != 0)
        projectedDirection.add(v, x);
    }
  }

  void Projection::liftHyperplane(const VectorRational& normal, const Rational& rhs,
    DVectorRational& liftedNormal, Rational& liftedRhs) const
  {
    liftedNormal.reDim(sourceSpace().dimension());
    liftedNormal.clear();
    liftedRhs = rhs;
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      const SVectorRational& row = _map[v];
      for (int p = row.size() - 1; p >= 0; --p)
        liftedNormal[row.index(p)] += normal[v] * row.value(p);
      liftedRhs -= normal[v] * _shift[v];
    }
  }

  ProjectedOracle::ProjectedOracle(const std::string& name,
    const Projection& projection, OracleBase* oracle)
    : OracleBase(name, projection), _projection(projection), _oracle(oracle),
    _liftedFace(NULL)
  {
    OracleBase::initializedSpace();

    _liftedVector.reDim(projection.sourceSpace().dimension(), false);
    _projectedVector.reDim(projection.dimension(), false);
  }


  ProjectedOracle::~ProjectedOracle()
  {
    if (_liftedFace)
      delete _liftedFace;
  }


  void ProjectedOracle::setFace(Face* newFace)
  {
    if (newFace == currentFace())
      return;

    _currentFace = newFace;
    if (_liftedFace)
      delete _liftedFace;

    if (newFace)
    {
      DVectorRational denseLiftedNormal;
      Rational liftedRhs;
      _projection.liftHyperplane(newFace->denseNormal(), newFace->rhs(), denseLiftedNormal,
        liftedRhs);
      DSVectorRational sparseLiftedNormal;
      sparseLiftedNormal = denseLiftedNormal;
      _liftedFace = new Face(_projection.sourceSpace(),
        LPRowRational(-infinity, sparseLiftedNormal, liftedRhs));
    }
    else
      _liftedFace = NULL;

    _oracle->setFace(_liftedFace);
  }

  void ProjectedOracle::maximize(OracleResult& result, const VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t maxHeuristic, std::size_t minHeuristic)
  {
    _liftedVector.clear();
    ObjectiveBound liftedObjectiveBound;
    liftedObjectiveBound.strict = objectiveBound.strict;
       _projection.liftHyperplane(objective, objectiveBound.value, _liftedVector,
      liftedObjectiveBound.value);

    OracleResult sourceResult;
    _oracle->maximize(sourceResult, _liftedVector, liftedObjectiveBound, maxHeuristic,
      minHeuristic);
    if (sourceResult.isFeasible())
    {
      result.buildStart(objective);
      for (std::size_t i = 0; i < sourceResult.points.size(); ++i)
      {
        _liftedVector.clear();
        _liftedVector.assign(*sourceResult.points[i].point);
        DSVectorRational* point = new DSVectorRational;
        _projection.projectPoint(_liftedVector, *point);
        result.buildAddPoint(point);
      }
      return result.buildFinish(sourceResult.heuristic(), true, false, true);
    }
    else if (sourceResult.isUnbounded())
    {
      result.buildStart(objective);
      for (std::size_t i = 0; i < sourceResult.directions.size(); ++i)
      {
        _liftedVector.clear();
        _liftedVector.assign(*sourceResult.directions[i].direction);
        DSVectorRational* direction = new DSVectorRational;
        _projection.projectDirection(_liftedVector, *direction);
        result.buildAddDirection(direction);
      }
      return result.buildFinish(sourceResult.heuristic(), true, false, true);
    }
    else
    {
      assert(sourceResult.isInfeasible());
      result.buildStart(objective);
      result.buildFinish(sourceResult.heuristic(), false, false, false);
    }
  }

}
