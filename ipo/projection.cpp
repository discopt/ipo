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
    const Vector& coefficients, const Rational& shift)
  {
    _variables.push_back(variableName);
    _map.push_back(coefficients);
    _shift.push_back(shift);
  }

  void Projection::addVariable(std::size_t sourceVariable, const soplex::Rational& shift)
  {
    _variables.push_back(_sourceSpace[sourceVariable]);
    _map.push_back(unitVector(sourceVariable));
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

  Vector Projection::projectPoint(const Vector& point) const
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

  Vector Projection::projectRay(const Vector& ray) const
  {
    VectorData* data = new VectorData(dimension());
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      Rational x = _map[v] * ray;
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

//   void Projection::liftHyperplane(const soplex::VectorRational& normal, const Rational& rhs,
//     soplex::DVectorRational& liftedNormal, Rational& liftedRhs) const
//   {
//     liftedNormal.reDim(sourceSpace().dimension());
//     liftedNormal.clear();
//     liftedRhs = rhs;
//     for (std::size_t v = 0; v < dimension(); ++v)
//     {
//       const SVectorRational& row = _map[v];
//       for (int p = row.size() - 1; p >= 0; --p)
//         liftedNormal[row.index(p)] += normal[v] * row.value(p);
//       liftedRhs -= normal[v] * _shift[v];
//     }
//   }

  LinearConstraint Projection::liftLinearConstraint(const LinearConstraint& projectedConstraint) const
  {
    soplex::DVectorRational liftedNormal(sourceSpace().dimension());
    Rational liftedRhs = projectedConstraint.rhs();
    for (std::size_t p = 0; p < projectedConstraint.normal().size(); ++p)
    {
      std::size_t v = projectedConstraint.normal().index(p);
      const Rational& x = projectedConstraint.normal().value(p);
      const Vector& row = _map[v];
      for (int q = row.size() - 1; q >= 0; --q)
        liftedNormal[row.index(q)] += x * row.value(q);
      liftedRhs -= x * _shift[v]; 
    }
    return LinearConstraint(projectedConstraint.type(), denseToVector(liftedNormal, true), liftedRhs);
  }

  ProjectedOracle::ProjectedOracle(const std::string& name,
    const Projection& projection, OracleBase* oracle)
    : OracleBase(name, projection), _projection(projection), _oracle(oracle)
  {
    OracleBase::initializedSpace();

    _liftedVector.reDim(projection.sourceSpace().dimension(), false);
    _projectedVector.reDim(projection.dimension(), false);
  }


  ProjectedOracle::~ProjectedOracle()
  {
    
  }


  void ProjectedOracle::setFace(const LinearConstraint& newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);

    if (newFace.definesCompleteFace())
      _liftedFace = completeFace();
    else
      _liftedFace = _projection.liftLinearConstraint(newFace);

    _oracle->setFace(_liftedFace);
  }
  
  std::size_t ProjectedOracle::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups)
  {
    Vector objectiveVector = denseToVector(objective);
    LinearConstraint objectiveConstraint('<' , objectiveVector, objectiveBound.value);
    LinearConstraint liftedObjectiveConstraint =  _projection.liftLinearConstraint(objectiveConstraint);

    Vector liftedObjective = liftedObjectiveConstraint.normal();
    Rational liftedRhs = liftedObjectiveConstraint.rhs();

    ObjectiveBound liftedObjectiveBound;
    liftedObjectiveBound.strict = objectiveBound.strict;
    
    OracleResult sourceResult;
    _oracle->maximize(sourceResult, liftedObjective, liftedObjectiveBound, minHeuristic, maxHeuristic);
    if (sourceResult.isFeasible())
    {
      for (std::size_t i = 0; i < sourceResult.points.size(); ++i)
      {
        Vector projectedPoint = _projection.projectPoint(sourceResult.points[i].vector);
        result.points.push_back(OracleResult::Point(projectedPoint));
      }
      checkDups = true;
    }
    else if (sourceResult.isUnbounded())
    {
      for (std::size_t i = 0; i < sourceResult.directions.size(); ++i)
      {
        Vector projectedRay = _projection.projectRay(sourceResult.directions[i].vector);
        result.directions.push_back(OracleResult::Direction(projectedRay));
      }
    }
  }

}
