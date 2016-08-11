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
     * \brief Constructs constraint of given \p type with given \p normal vector, right-hand side \p rhs.
     * 
     * Constructs constraint of given \p type with given \p normal vector, right-hand side \p rhs.
     */

    LinearConstraint(char type, Vector& normal, const Rational& rhs);

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    ~LinearConstraint();

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
  
} /* namespace ipo */

#endif /* IPO_LINEAR_CONSTRAINT_H_ */