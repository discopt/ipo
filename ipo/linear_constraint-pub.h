#ifndef IPO_LINEAR_CONSTRAINT_H_
#define IPO_LINEAR_CONSTRAINT_H_

#include "common.h"
#include "vectors.h"

namespace ipo {

  /**
   * \brief Linear constraint (inequality or equation).
   *
   * Linear constraint (inequality or equation) whose normal vector is stored in a \ref Vector object.
   */

  class LinearConstraint
  {
  public:
    /**
     * \brief Constructs constraint <0,x> <= 0.
     *
     * Constructs constraint <0,x> <= 0.
     */

    LinearConstraint();

    /**
     * \brief Constructs constraint of given \p type with given \p normal vector, right-hand side \p rhs.
     *
     * Constructs constraint of given \p type with given \p normal vector, right-hand side \p rhs.
     */

    LinearConstraint(char type, const Vector& normal, const Rational& rhs);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    ~LinearConstraint();

    /**
     * \brief Returns true iff the constraints are equal.
     *
     * Returns true iff the constraints are equal. Does not take into account the scaling.
     */

    inline bool operator==(const LinearConstraint& other) const
    {
      if (_type != other._type)
        return false;
      if (_rhs != other._rhs)
        return false;
      return _normal == other._normal;
    }

    /**
     * \brief Fast ordering comparison.
     *
     * Fast ordering comparison. First compares the right-hand sides and types and finally the normal vectors.
     */

    inline bool operator<(const LinearConstraint& other) const
    {
      if (_type < other._type)
        return true;
      if (_type > other._type)
        return false;
      if (_rhs < other._rhs)
        return true;
      if (_rhs > other._rhs)
        return false;
      return _normal < other._normal;
    }

    /**
     * \brief Returns true for equations.
     *
     * Returns true for equations.
     */

    inline bool isEquation() const
    {
      return _type == '=';
    }

    /**
     * \brief Returns the type of the constraint, i.e., '<', '>' or '='.
     */

    inline char type() const
    {
      return _type;
    }

    /**
     * \brief Returns the normal vector of the constraint.
     *
     * Returns the normal vector of the constraint.
     */

    inline const Vector& normal() const
    {
      return _normal;
    }

    /**
     * \brief Returns the right-hande side of the constraint.
     *
     * Returns the right-hande side of the constraint.
     */

    inline const Rational& rhs() const
    {
      return _rhs;
    }

    /**
     * \brief Returns the maximum norm of the normal vector.
     *
     * Returns the maximum norm of the normal vector, i.e., the largest absolute value of an entry.
     */

    inline Rational getMaximumNorm() const
    {
      Rational result = 0;
      for (std::size_t p = 0; p < _normal.size(); ++p)
      {
        const Rational& x = _normal.value(p);
        if (x > 0 && x > result)
          result = x;
        if (x < 0 && x < -result)
          result = -x;
      }
      return result;
    }

    /**
     * \brief Returns true iff this inequality constraint defines the complete face.
     *
     * Returns true iff the constraint is an inequality constraint that defines the complete face.
     */

    inline bool definesCompleteFace() const
    {
      if (_type == '=')
        return false;
      else if (_type == '<')
        return _normal.size() == 0 && _rhs >= 0;
      else
        return _normal.size() == 0 && _rhs <= 0;
    }

    /**
     * \brief Returns true iff this inequality constraint defines the empty face.
     *
     * Returns true iff the constraint is an inequality constraint that defines the empty face.
     */

    inline bool definesEmptyFace() const
    {
      if (_type == '=')
        return false;
      else if (_type == '<')
        return _normal.size() == 0 && _rhs < 0;
      else
        return _normal.size() == 0 && _rhs > 0;
    }

    /**
     * \brief Returns true iff this inequality constraint defines a trivial face.
     *
     * Returns true iff the constraint is an inequality constraint that defines the trivial (complete or empty) face.
     */

    inline bool definesTrivialFace() const
    {
      return _type != '=' && _normal.size() ==0;
    }

    /**
     * \brief Returns an integer indicating where the point lies relative to the halfspace/hyperplane.
     *
     * Returns 0 iff the point lies on the hyperplane. It returns a positive number iff the point lies in the interior of the
     * halfspace (happens never for equations). Otherwise, returns a negative number.
     */

    int evaluatePoint(const Vector& point) const;

    /**
     * \brief Returns an integer indicating where the ray lies relative to the halfspace/hyperplane.
     *
     * Returns 0 iff the ray is parallel to the hyperplane. It returns a positive number iff the ray has positive scalar product
     * with the halfspace's outer normal vector (happens never for equations). Otherwise, returns a negative number.
     */

    int evaluateRay(const Vector& ray) const;

  private:
    char _type; // '<', '>' or '='
    Vector _normal; // Normal vector of constraint.
    Rational _rhs; // Right-hand side of constraint.
  };

  LinearConstraint completeFaceConstraint();
  LinearConstraint emptyFaceConstraint();

  LinearConstraint operator+(const LinearConstraint& a, const LinearConstraint& b);
  LinearConstraint operator-(const LinearConstraint& a, const LinearConstraint& b);

  typedef std::vector<LinearConstraint> AffineOuterDescription;
  

} /* namespace ipo */

#endif /* IPO_LINEAR_CONSTRAINT_H_ */
