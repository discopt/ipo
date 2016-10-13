#include "linear_constraint.h"

namespace ipo {

  LinearConstraint::LinearConstraint()
    : _type('<'), _normal(), _rhs(0)
  {

  }

  LinearConstraint::LinearConstraint(char type, const Vector& normal, const Rational& rhs)
    : _type(type), _normal(normal), _rhs(rhs)
  {
    assert(type == '<' || type == '>' || type == '=');
  }

  LinearConstraint::~LinearConstraint()
  {

  }

  int LinearConstraint::evaluatePoint(const Vector& point) const
  {
    int cmp = compareRational(_rhs, _normal * point);
    if (_type == '<')
      return cmp;
    else if (_type == '>')
      return -cmp;
    else if (cmp == 0)
      return 0;
    else
      return -1;
  }

  int LinearConstraint::evaluateRay(const Vector& ray) const
  {
    int cmp = compareRational(0, _normal * ray);
    if (_type == '<')
      return cmp;
    else if (_type == '>')
      return -cmp;
    else if (cmp == 0)
      return 0;
    else
      return -1;
  }

  LinearConstraint completeFace()
  {
    return LinearConstraint('<', Vector(), Rational(0));
  }

  LinearConstraint emptyFace()
  {
    return LinearConstraint('<', Vector(), Rational(-1));
  }

  LinearConstraint operator+(const LinearConstraint& a, const LinearConstraint& b)
  {
    if (a.type() == '>')
    {
      if (b.type() == '>')
        return addScaled('<', a, -1, b, -1);
      else
        return addScaled('<', a, -1, b, 1);
    }
    else
    {
      assert(a.type() == '=' || a.type() == '<');
      if (b.type() == '>')
        return addScaled('<', a, 1, b, -1);
      else
        return addScaled(b.type(), a, 1, b, 1);
    }
  }

  LinearConstraint operator-(const LinearConstraint& a, const LinearConstraint& b)
  {
    if (a.type() == '>')
    {
      if (b.type() == '>')
        return addScaled('<', a, -1, b, -1);
      else if (b.type() == '=')
        return addScaled('<', a, -1, b, -1);
      else
        return addScaled('<', a, -1, b, 1);
    }
    else if (a.type() == '<')
    {
      if (b.type() == '>')
        return addScaled('<', a, 1, b, -1);
      else if (b.type() == '<')
        return addScaled('<', a, 1, b, 1);
      else
        return addScaled('<', a, 1, b, -1);
    }
    else
    {
      assert(a.type() == '=');
      if (b.type() == '>')
        return addScaled('<', a, 1, b, -1);
      else if (b.type() == '<')
        return addScaled('<', a, 1, b, -1);
      else
        return addScaled(b.type(), a, 1, b, -1);
    }
  }

  LinearConstraint addScaled(char type, const LinearConstraint& a, int scaleA, const LinearConstraint& b, int scaleB)
  {
    return LinearConstraint(b.type(), addScaled(a.normal(), scaleA, b.normal(), scaleB), scaleA * a.rhs() + scaleB * b.rhs());
  }

  void addToLP(soplex::SoPlex& spx, const LinearConstraint& constraint)
  {
    const Rational& rhs = constraint.rhs();
    soplex::DSVectorRational normal(constraint.normal().size());
    vectorToSparse(constraint.normal(), normal);
    if (constraint.type() == '<')
      spx.addRowRational(soplex::LPRowRational(-soplex::infinity, normal, rhs));
    else if (constraint.type() == '>')
      spx.addRowRational(soplex::LPRowRational(rhs, normal, soplex::infinity));
    else
      spx.addRowRational(soplex::LPRowRational(rhs, normal, rhs));
  }

  void addToLP(soplex::SoPlex& spx, const std::vector< LinearConstraint >& constraints)
  {
    soplex::DSVectorRational normal(spx.numColsRational());
    for (std::size_t i = 0; i < constraints.size(); ++i)
    {
      const Rational& rhs = constraints[i].rhs();
      vectorToSparse(constraints[i].normal(), normal);
      if (constraints[i].type() == '<')
        spx.addRowRational(soplex::LPRowRational(-soplex::infinity, normal, rhs));
      else if (constraints[i].type() == '>')
        spx.addRowRational(soplex::LPRowRational(rhs, normal, soplex::infinity));
      else
        spx.addRowRational(soplex::LPRowRational(rhs, normal, rhs));
    }
  }

  LinearConstraint integralScaled(const LinearConstraint& constraint)
  {
    // Compute scaling factor.

    IntegralScaler scaler;
    for (std::size_t p = 0; p < constraint.normal().size(); ++p)
      scaler(constraint.normal().value(p));
    Rational factor = scaler.factor();

    // Scale it.

    VectorData* data = new VectorData(constraint.normal().size());
    for (std::size_t p = 0; p < constraint.normal().size(); ++p)
      data->add(constraint.normal().index(p), factor * constraint.normal().value(p));
    return LinearConstraint(constraint.type(), Vector(data), factor * constraint.rhs());
  }

  void scaleIntegral(LinearConstraint& constraint)
  {
    constraint = integralScaled(constraint);
  }

  void scaleIntegral(std::vector<LinearConstraint>& constraints)
  {
    for (std::size_t i = 0; i < constraints.size(); ++i)
      scaleIntegral(constraints[i]);
  }


} /* namespace ipo */
