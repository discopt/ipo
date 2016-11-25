#include "projection.h"

using namespace soplex;

namespace ipo {

  ProjectionData::ProjectionData(const Space& sourceSpace)
    : _sourceSpace(sourceSpace)
  {
    _imageSpaceData = new SpaceData;
    _imageSpaceData->markUsed();
  }

  ProjectionData::ProjectionData(const Space& sourceSpace, const std::vector< std::size_t >& variableSubset)
    : _sourceSpace(sourceSpace)
  {
    std::vector<std::string> variableNames;
    for (std::size_t iv = 0; iv < variableSubset.size(); ++iv)
      variableNames.push_back(sourceSpace[variableSubset[iv]]);

    for (std::size_t iv = 0; iv < variableSubset.size(); ++iv)
      _map.push_back(unitVector(variableSubset[iv]));

    _shift.resize(variableNames.size(), Rational(0));

    _imageSpaceData = new SpaceData(variableNames);
    _imageSpaceData->markUsed();
  }

  ProjectionData::~ProjectionData()
  {
    _imageSpaceData->unmarkUsed();
  }

  bool ProjectionData::operator==(const ProjectionData& other) const
  {
    if (_sourceSpace != other._sourceSpace)
      return false;

    if (*_imageSpaceData != *other._imageSpaceData)
      return false;

    if (_shift != other._shift)
      return false;

    for (std::size_t v = 0; v < _map.size(); ++v)
    {
      if (_map[v] != other._map[v])
        return false;
    }
    return true;
  }

  void ProjectionData::addVariable(const std::string& variableName, const Vector& coefficients, const Rational& shift)
  {
    _imageSpaceData->addVariable(variableName);
    _map.push_back(coefficients);
    _shift.push_back(shift);
  }

  void ProjectionData::addVariable(std::size_t sourceVariable, const Rational& shift)
  {
    _imageSpaceData->addVariable(_sourceSpace[sourceVariable]);
    _map.push_back(unitVector(sourceVariable));
    _shift.push_back(shift);
  }

  void ProjectionData::unmarkUsed()
  {
    _usage--;
    if (_usage == 0)
      delete this;
  }

  Vector Projection::projectPoint(const Vector& point) const
  {
    VectorData* data = new VectorData(imageSpace().dimension());
    for (std::size_t v = 0; v < imageSpace().dimension(); ++v)
    {
      Rational x = shift(v) + row(v) * point;
      if (x != 0)
        data->add(v, x);
    }
    return Vector(data);
  }

  Vector Projection::projectRay(const Vector& ray) const
  {
    VectorData* data = new VectorData(imageSpace().dimension());
    for (std::size_t v = 0; v < imageSpace().dimension(); ++v)
    {
      Rational x = row(v) * ray;
      if (x != 0)
        data->add(v, x);
    }
    return Vector(data);
  }

  LinearConstraint Projection::projectLinearConstraint(const LinearConstraint& constraint) const
  {
    std::cerr << "Projection of linear constraints not implemented - claiming not projectible!" << std::endl;

    return completeFaceConstraint();
  }

  LinearConstraint Projection::liftLinearConstraint(const LinearConstraint& constraint) const
  {
    soplex::DVectorRational liftedNormal(sourceSpace().dimension());
    Rational liftedRhs = constraint.rhs();
    for (std::size_t p = 0; p < constraint.normal().size(); ++p)
    {
      std::size_t variable = constraint.normal().index(p);
      const Rational& x = constraint.normal().value(p);
      const Vector& rowVector = row(variable);
      for (int q = rowVector.size() - 1; q >= 0; --q)
        liftedNormal[rowVector.index(q)] += x * rowVector.value(q);
      liftedRhs -= x * shift(variable);
    }
    return LinearConstraint(constraint.type(), denseToVector(liftedNormal, true), liftedRhs);
  }

  ProjectionOracle::ProjectionOracle(const std::string& name, const Projection& projection,
    const std::shared_ptr<OracleBase>& oracle)
    : OracleBase(name), _projection(projection), _oracle(oracle)
  {
    OracleBase::initializeSpace(projection.imageSpace());

    _projectedVector.reDim(space().dimension(), false);
  }

  ProjectionOracle::ProjectionOracle(const Projection& projection, const std::shared_ptr<OracleBase>& oracle)
    : OracleBase("Projection(" + oracle->name() + ")"), _projection(projection), _oracle(oracle)
  {
    OracleBase::initializeSpace(projection.imageSpace());

    _projectedVector.reDim(space().dimension(), false);
  }

  ProjectionOracle::~ProjectionOracle()
  {

  }


  void ProjectionOracle::setFace(const LinearConstraint& newFace)
  {
    if (newFace == currentFace())
      return;

    OracleBase::setFace(newFace);

    if (newFace.definesCompleteFace())
      _liftedFace = completeFaceConstraint();
    else
      _liftedFace = _projection.liftLinearConstraint(newFace);

    _oracle->setFace(_liftedFace);
  }

  HeuristicLevel ProjectionOracle::maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
    const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort, bool& checkDups)
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
      for (std::size_t i = 0; i < sourceResult.rays.size(); ++i)
      {
        Vector projectedRay = _projection.projectRay(sourceResult.rays[i].vector);
        result.rays.push_back(OracleResult::Ray(projectedRay));
      }
    }
  }

}
