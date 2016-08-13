#ifndef IPO_PROJECTION_H_
#define IPO_PROJECTION_H_

#include "oracles.h"

namespace ipo {

  /**
   * \brief An affine projection map, together with the image space.
   *
   * Defines an affine projection \f$ \pi(x) = Ax + b \f$ together with the image space.
   * Its main purpose is to create a \ref ProjectedOracle from another oracle.
   */

  class Projection : public Space
  {
  public:
    /**
     * \brief Constructs a projection of \p sourceSpace into the 0-dimensional image space.
     *
     * Constructs a projection of \p sourceSpace into the 0-dimensional image space. Use
     * \ref addVariable() to increase its dimension.
     *
     * \param sourceSpace Source space.
     */

    Projection(const Space& sourceSpace);

    /**
     * \brief Constructs an orthogonal projection for a given space.
     *
     * Constructs an orthogonal projection of the \c space onto a subset of the variables.
     *
     * \param sourceSpace    Source space of the projection.
     * \param variableSubset Indices of the subset of variables to project to.
     */

    Projection(const Space& sourceSpace, const std::vector<std::size_t>& variableSubset);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~Projection();

    /**
     * \brief Returns the source space.
     *
     * Returns a const reference to the source space.
     */

    const Space& sourceSpace() const
    {
      return _sourceSpace;
    }

    /**
     * \brief Adds a variable to the image space.
     *
     * Adds a variable to the image space by adding a given row to the projection matrix and to the
     * constant term.
     *
     * \param variableName Name of the new variable.
     * \param coefficients Sparse rational row of the projection matrix \f$ A \f$.
     * \param shift        Rational entry for shift vector \f$ b \f$.
     */

    void addVariable(const std::string& variableName, const Vector& coefficients, const Rational& shift = Rational(0));

    /**
     * \brief Copies a source variable to the image space.
     *
     * Copies a source variable to the image space.
     *
     * \param sourceVariable Index of the variable in the source space.
     * \param shift Rational entry for shift vector \f$ b \f$.
     */

    void addVariable(std::size_t sourceVariable, const Rational& shift = Rational(0));

    /**
     * \brief Projects a \c point.
     *
     * Projects a \c point.
     */

    Vector projectPoint(const soplex::VectorRational& point) const;

    /**
     * \brief Projects a \c point.
     *
     * Projects a \c point.
     */

    Vector projectPoint(const Vector& point) const;

    /**
     * \brief Projects a \c ray.
     *
     * Projects a \c ray.
     */

    Vector projectRay(const soplex::VectorRational& ray) const;

    /**
     * \brief Projects a \c ray.
     *
     * Projects a \c ray.
     */

    Vector projectRay(const Vector& ray) const;

    /**
     * \brief Projects a linear \p constraint if possible.
     *
     * If possible, projects a linear constraint into the projected space. Otherwise, returns the completeFace().
     *
     * Let the projection map be \f$ y = Ax + b \f$ and let \f$ \left<a,x\right> = \beta \f$ be the equation defining the 
     * given (boundary) hyperplane. The hyperplane defined by \f$ \left<c,x\right> = \delta \f$ is its projection if
     * \f$ c^\intercal A  = a^\intercal \f$ and \f$\delta - c^\intercal b = \beta\f$ hold. If this system has no solution, then
     * completeFace() is returned.
     */

    LinearConstraint projectLinearConstraint(const LinearConstraint& constraint) const;

    /**
     * \brief Lifts a linear \p constraint into the source space.
     *
     * Lifts a linear \p constraint into the source space.
     *
     * Let the projection map be \f$ y = Ax + b \f$ and let \f$ \left<a,y\right> = \beta \f$ be the equation defining the 
     * (boundary) hyperplane. Its lifted version is then defined by \f$ \left<a^\intercal A,x\right> = \beta - a^\intercal b \f$.
     */
  
    LinearConstraint liftLinearConstraint(const LinearConstraint& constraint) const;

    /**
     * \brief Returns the dimension of the image space.
     *
     * Returns the dimension of the image space, i.e., the number of variables.
     */

    inline std::size_t dimension() const
    {
      return _map.size();
    }

    /**
     * \brief Returns the name of the image variable indexed by \c var.
     *
     * Returns the name of the image variable indexed by \c var.
     */

    inline const std::string& operator[](std::size_t var) const
    {
      return _variables[var];
    }

    /**
     * \brief Returns the name of the image variable indexed by \c var.
     *
     * Returns the name of the image variable indexed by \c var.
     */

    inline const std::string& variable(std::size_t var) const
    {
      return _variables[var];
    }

    /**
     * \brief Returns a row of the projection matrix.
     *
     * Returns a row of the projection matrix \f$ A \f$ for a given \p variable
     * as a const reference to a sparse rational vector.
     */

    inline const Vector& map(std::size_t variable) const
    {
      return _map[variable];
    }

    /**
     * \brief Returns an entry of the projection vector.
     *
     * Returns an entry of the projection vector \f$ b \f$ for a given \p variable.
     */

    inline const Rational& shift(std::size_t var) const
    {
      return _shift[var];
    }

  private:

    /**
     * \brief The default constructor is disabled.
     *
     * The default constructor is disabled.
     */

    Projection();

  protected:
    const Space& _sourceSpace;
    std::vector<Vector> _map; // Rows of the projection matrix.
    std::vector<Rational> _shift; // Entries of the projection vector.
  };

  /**
   * \brief An oracle for the projection of a polyhedron defined by another oracle.
   *
   * Defines an optimization oracle for the projection \f$ \pi(P) \f$ of a polyhedron \f$ P \f$ defined by another oracle
   * for given affine projection map \f$ \pi \f$.
   *
   * \sa Projection
   */

  class ProjectedOracle: public OracleBase
  {
  public:
    /**
     * \brief Constructor for given projection and source oracle.
     *
     * Constructor for given \c projection and source \c oracle.
     *
     * \param name       Name of the new oracle.
     * \param projection Projection map.
     * \param oracle     Oracle in the source space.
     */

    ProjectedOracle(const std::string& name, const Projection& projection,
      OracleBase* oracle);

    /**
     * \brief Destructor.
     */

    virtual ~ProjectedOracle();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation creates a new face which represents the given face in the preimage
     * space.
     */

    virtual void setFace(const LinearConstraint& newFace = completeFace());
    
    /**
     * \brief Returns the ambient \c space.
     *
     * Returns a reference to the ambient \c space.
     */

    inline const Projection& space() const
    {
      return _projection;
    }

  protected:

    /**
     * \brief Oracle's implementation to maximize the dense rational \p objective.
     *
     * This method is called by maximizeController() and contains the implementation of the oracle.
     * 
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param sort           Set this variable to true if points must be sorted.
     * \param checkDups      Set this variable to true if points or rays must be checked for duplicates.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual std::size_t maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);

  protected:

    const Projection& _projection; // Projection.
    OracleBase* _oracle; // Source oracle.
    soplex::DVectorRational _projectedVector; // Projected point or ray.
    LinearConstraint _liftedFace;
  };

} /* namespace ipo */



#endif /* IPO_PROJECTION_H_ */
