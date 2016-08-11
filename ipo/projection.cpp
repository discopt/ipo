#include "projection.h"

using namespace soplex;

namespace ipo {

  Projection::Projection(const Space& sourceSpace)
    : _sourceSpace(sourceSpace)
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

  Projection::~Projection()
  {

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
  
  Vector Projection::projectPoint(const VectorRational& point) const
  {
    VectorData* data = new VectorData(dimension());
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      Rational x = _shift[v] + _map[v] * point;
      if (x != 0)
        data->add(v, x);
    }
    return Vector(data);
  }

  Vector Projection::projectDirection(const VectorRational& direction) const
  {
    VectorData* data = new VectorData(dimension());
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      Rational x = _map[v] * direction;
      if (x != 0)
        data->add(v, x);
    }
    return Vector(data);
  }
  
  bool Projection::projectHyperplane(const soplex::VectorRational& normal, const Rational& rhs,
    soplex::DVectorRational& projectedNormal,
    Rational& projectedRhs) const
  {
    std::cerr << "Projection of hyperplanes not implemented - claiming not projectible!" << std::endl;

    return false;
  }

  void Projection::liftHyperplane(const soplex::VectorRational& normal, const Rational& rhs,
    soplex::DVectorRational& liftedNormal, Rational& liftedRhs) const
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

    OracleBase::setFace(newFace);

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
  
  std::size_t ProjectedOracle::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    _liftedVector.clear();
    ObjectiveBound liftedObjectiveBound;
    liftedObjectiveBound.strict = objectiveBound.strict;
    _projection.liftHyperplane(objective, objectiveBound.value, _liftedVector, liftedObjectiveBound.value);
    
    OracleResult sourceResult;
    _oracle->maximize(sourceResult, _liftedVector, liftedObjectiveBound, minHeuristic, maxHeuristic);
    if (sourceResult.isFeasible())
    {
      for (std::size_t i = 0; i < sourceResult.points.size(); ++i)
      {
        vectorToDense(sourceResult.points[i].vector, _liftedVector);
        Vector projectedPoint = _projection.projectPoint(_liftedVector);
        result.points.push_back(OracleResult::Point(projectedPoint));
      }
      checkDups = true;
    }
    else if (sourceResult.isUnbounded())
    {
      for (std::size_t i = 0; i < sourceResult.directions.size(); ++i)
      {
        vectorToDense(sourceResult.directions[i].vector, _liftedVector);
        Vector projectedDirection = _projection.projectDirection(_liftedVector);
        result.directions.push_back(OracleResult::Direction(projectedDirection));
      }
    }
  }

}
