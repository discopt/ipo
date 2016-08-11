#include "linear_constraint.h"

namespace ipo {
    
  LinearConstraint::LinearConstraint(char type, Vector& normal, const Rational& rhs)
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

} /* namespace ipo */