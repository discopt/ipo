#include "linear_constraint.h"

namespace ipo {
  
  LinearConstraint::LinearConstraint()
    : _type(','), _normal(), _rhs(0)
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


} /* namespace ipo */