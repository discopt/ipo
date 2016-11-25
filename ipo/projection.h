#ifndef IPO_PROJECTION_H_
#define IPO_PROJECTION_H_

#include "oracles.h"

namespace ipo {

  class ProjectionData
  {
  public:
    ProjectionData(const Space& sourceSpace);
    ProjectionData(const Space& sourceSpace, const std::vector<std::size_t>& variableSubset);
    ~ProjectionData();

    bool operator==(const ProjectionData& other) const;

    inline const Space& sourceSpace() const
    {
      return _sourceSpace;
    }

    inline Space imageSpace()
    {
      return Space(_imageSpaceData);
    }

    inline const Vector& row(std::size_t variable) const
    {
      return _map[variable];
    }

    inline const Rational& shift(std::size_t variable) const
    {
      return _shift[variable];
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
     * \brief Increases the usage counter by 1.
     *
     * Increases the usage counter by 1.
     */

    inline void markUsed()
    {
      _usage++;
    }

    /**
     * \brief Decreases the usage counter by 1.
     *
     * Decreases the usage counter by 1 and frees it if the latter reached 0.
     */

    void unmarkUsed();

  protected:
    SpaceData* _imageSpaceData;
    const Space _sourceSpace;
    std::vector<Vector> _map; // Rows of the projection matrix.
    std::vector<Rational> _shift; // Entries of the projection vector.
    std::size_t _usage;
  };

  class Projection
  {
  public:
    inline Projection(ProjectionData* data)
      : _data(data), _imageSpace(_data->imageSpace())
    {
      _data->markUsed();
    }

    inline Projection(const Projection& other)
      : _data(other._data), _imageSpace(other.imageSpace())
    {
      _data->markUsed();
    }

    ~Projection()
    {
      _data->unmarkUsed();
    }

    inline Projection& operator=(const Projection& other)
    {
      if (_data != other._data)
      {
        _data->unmarkUsed();
        _data = other._data;
        _data->markUsed();
        _imageSpace = other.imageSpace();
      }
      return *this;
    }

    inline bool operator==(const Projection& other)
    {
      if (_data == other._data)
        return true;

      return *_data == *other._data;
    }

    inline const Space& sourceSpace() const
    {
      return _data->sourceSpace();
    }

    inline const Space& imageSpace() const
    {
      return _imageSpace;
    }

    inline const Vector& row(std::size_t variable) const
    {
      return _data->row(variable);
    }

    inline const Rational& shift(std::size_t variable) const
    {
      return _data->shift(variable);
    }

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

  protected:
    ProjectionData* _data;
    Space _imageSpace;
  };

  /**
   * \brief An oracle for the projection of a polyhedron defined by another oracle.
   *
   * Defines an optimization oracle for the projection \f$ \pi(P) \f$ of a polyhedron \f$ P \f$ defined by another oracle
   * for given affine projection map \f$ \pi \f$.
   *
   * \sa Projection
   */

  class ProjectionOracle: public OracleBase
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

    ProjectionOracle(const std::string& name, const Projection& projection, const std::shared_ptr<OracleBase>& oracle);

    /**
     * \brief Constructor for given projection and source oracle.
     *
     * Constructor for given \p projection and source \p oracle. The oracle is named "Projection(\p oracle->name())".
     *
     * \param projection Projection map.
     * \param oracle     Oracle in the source space.
     */

    ProjectionOracle(const Projection& projection, const std::shared_ptr<OracleBase>& oracle);

    /**
     * \brief Destructor.
     */

    virtual ~ProjectionOracle();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation creates a new face which represents the given face in the preimage space.
     */

    virtual void setFace(const LinearConstraint& newFace = completeFaceConstraint());

    /**
     * \brief Returns the ambient \c space.
     *
     * Returns a const-reference to the ambient \c space.
     */

    inline const Space& space() const
    {
      return _projection.imageSpace();
    }

    /**
     * \brief Returns the associated projection.
     *
     * Returns a const-reference to the associated projection.
     */

    inline const Projection& projection() const
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

    virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
      bool& checkDups);

  protected:

    Projection _projection; // Projection map and image space.
    std::shared_ptr<OracleBase> _oracle; // Source oracle.
    soplex::DVectorRational _projectedVector; // Projected point or ray.
    LinearConstraint _liftedFace;
  };

} /* namespace ipo */



#endif /* IPO_PROJECTION_H_ */
