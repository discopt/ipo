#ifndef IPO_MIXED_INTEGER_PROGRAM_H_
#define IPO_MIXED_INTEGER_PROGRAM_H_

#include <string>
#include <vector>

#include <soplex.h>

#ifdef WITH_SCIP
#include <scip/scip.h>
#include "scip_oracles.h"
#endif

#include "ipo.h"
#include "oracles.h"

namespace ipo {

  /**
   * \brief A mixed-integer program.
   *
   * A mixed-integer program.
   * It is typically used to create oracles for optimizing over it,
   * e.g., using \ref SCIPOptimizationOracle and \ref ExactSCIPOptimizationOracle.
   */

  class MixedIntegerProgram
  {
  public:
#ifdef WITH_SCIP
    /**
     * \brief Constructs a \c MixedIntegerProgram from a \c SCIP instance.
     *
     * Constructs a \c MixedIntegerProgram from a \c SCIP instance.
     * Only explicitly stated linear constraints of the instance are considered.
     */

    MixedIntegerProgram(SCIP* scip);
#endif

    /**
     * \brief Destructor.
     */

    virtual ~MixedIntegerProgram();

    /**
     * \brief Returns the number of variables.
     *
     * Returns the number of variables.
     */

    inline std::size_t numVariables() const
    {
      return _variableNames.size();
    }

    /**
     * \brief Returns the number of constraints.
     *
     * Returns the number of constraints.
     */

    inline std::size_t numConstraints() const
    {
      return _constraintNames.size();
    }

    /**
     * \brief Returns the columns.
     *
     * Returns a const reference to the columns.
     */

    inline const soplex::LPColSetRational& columns() const
    {
      return _columns;
    }

    /**
     * \brief Returns the rows.
     *
     * Returns a const reference to the rows.
     */

    inline const soplex::LPRowSetRational& rows() const
    {
      return _rows;
    }

    /**
     * \brief Returns whether the given variable is integral.
     *
     * Returns whether the given \c variable has an integrality constraint.
     */

    inline bool isIntegral(std::size_t variable) const
    {
      return _integrality[variable];
    }

    /**
     * \brief Returns the name of the given variable.
     *
     * Returns the name of the given \c variable.
     */

    inline const std::string& variableName(std::size_t variable) const
    {
      return _variableNames[variable];
    }

    /**
     * \brief Returns the name of the given constraint.
     *
     * Returns the name of the given \c constraint.
     */

    inline const std::string& constraintName(std::size_t constraint) const
    {
      return _constraintNames[constraint];
    }

    /**
     * \brief Returns true if the given \c point satisfies all bound constraints.
     *
     * Returns true if the given \c point satisfies all bound constraints.
     */

    bool checkPointBounds(const soplex::SVectorRational* point) const;

    /**
     * \brief Returns true if the given \c point satisfies all row constraints.
     *
     * Returns true if the given \c point satisfies all row constraints.
     */

    bool checkPointConstraints(const soplex::SVectorRational* point);

    /**
     * \brief Returns true if the given \c point satisfies all integrality constraints.
     *
     * Returns true if the given \c point satisfies all integrality constraints.
     */

    bool checkPointIntegral(const soplex::SVectorRational* point) const;

    /**
     * \brief Returns true if the given \c point is feasible.
     *
     * Returns true if the given \c point is feasible.
     */

    bool checkPoint(const soplex::SVectorRational* point);

    /**
     * \brief Returns true if the given \c ray satisfies all bound constraints.
     *
     * Returns true if the given \c ray satisfies all bound constraints.
     */

    bool checkRayBounds(const soplex::SVectorRational* ray) const;

    /**
     * \brief Returns true if the given \c ray satisfies all row constraints.
     *
     * Returns true if the given \c ray satisfies all row constraints.
     */

    bool checkRayConstraints(const soplex::SVectorRational* ray);

    /**
     * \brief Returns true if the given \c ray is feasible.
     *
     * Returns true if the given \c ray is feasible.
     */

    bool checkRay(const soplex::SVectorRational* ray);

    /**
     * \brief Restricts the MIP to the given face.
     *
     * Restricts the MIP to the given face by adding an equation constraint.
     */

    void faceEnabled(Face* face);

    /**
     * \brief Removes a face restriction.
     *
     * Removes a face restriction added by \ref faceEnabled().
     */

    void faceDisabled(Face* face);

    /**
     * \brief Returns all row constraints.
     *
     * Returns all row constraints or only those that are inequalities or equations,
     * respectively.
     *
     * \param rows
     *   Rows structure to write to.
     * \param inequalities
     *   Whether to extract inequalities.
     * \param equations
     *   Whether to extract equations.
     * \param names
     *   If not \c NULL, writes corresponding row names.
     */

    void getConstraints(soplex::LPRowSetRational& rows, bool inequalities, bool equations,
        std::vector<std::string>* names = NULL);

    /**
     * \brief Extracts equations that correspond to fixed variables.
     *
     * Extracts equations that correspond to fixed variables.
     *
     * \param rows
     *   Rows structure to write to.
     * \param names
     *   If not \c NULL, writes corresponding row names.
     */

    void getFixedVariableEquations(soplex::LPRowSetRational& rows, std::vector<std::string>* names = NULL);

  protected:
    soplex::LPColSetRational _columns; // Columns
    soplex::LPRowSetRational _rows; // Rows
    std::vector<std::string> _variableNames; // Variable names
    std::vector<std::string> _constraintNames; // Row names
    std::vector<bool> _integrality; // Integrality constraints
    soplex::DVectorRational _worker; // Temporary dense rational vector.
    Face* _face; // Currently active face.
  };

  /**
   * \brief An oracle that postprocesses solutions from an inexact oracle.
   *
   * An oracle that can postprocess floating-point solutions
   * of a mixed-integer program to get correct rational ones.
   */

  class MixedIntegerProgramCorrectorOracle: public FaceOptimizationOracleBase
  {
  public:
    /**
     * \brief Constructor.
     *
     * Constructs the oracle for a given \c mip.
     * It calls another \c inexact oracle.
     * If a returned solution (point or direction) is
     * infeasible or if \c correctAlways is \c true,
     * then it rounds the integer variables of the \c mip,
     * and recomputes optimal continuous variables
     * by solving a linear program.
     */

    MixedIntegerProgramCorrectorOracle(const std::string& name, MixedIntegerProgram& mip,
        FaceOptimizationOracleBase* inexact, bool correctAlways = true);

    /**
     * \brief Destructor.
     */

    virtual ~MixedIntegerProgramCorrectorOracle();

  protected:

    /**
     * \brief Corrects a given point.
     *
     * Corrects a given point by rounding integer variables
     * and solving an LP for the continuous ones.
     */

    soplex::DSVectorRational* correctPoint(const soplex::SVectorRational* point,
        const soplex::VectorRational& objective);

    /**
     * \brief Corrects a given ray.
     *
     * Corrects a given ray (unbounded direction)
     * by solving an LP.
     */

    soplex::DSVectorRational* correctRay(const soplex::SVectorRational* ray, const soplex::VectorRational& objective);

    /**
     * \brief Actual implementation.
     *
     * The actual implementation calls the inexact oracle,
     * checks feasibility of the returned points or directions,
     * and finally corrects them, if required.
     */

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    /**
     * \brief Method that is called when a new face is activated.
     *
     * Method that is called when a new \c face is activated.
     * If a (non-trivial) face was active before,
     * it is ensured that \ref faceDisabled() is called before.
     *
     * The implementation calls the corresponding method
     * for the inexact oracle and adds an equation constraint.
     */

    virtual void faceEnabled(Face* face);

    /**
     * \brief Method that is called when a face is deactivated.
     *
     * Method that is called when a face is deactivated.
     *
     * The implementation removes the constraint that was added in \ref faceEnabled()
     * and calls the corresponding method for the inexact oracle.
     */

    virtual void faceDisabled(Face* face);

  protected:

    MixedIntegerProgram& _mip; // Mixed-integer program.
    FaceOptimizationOracleBase* _inexact; // Inexact oracle.
    soplex::SoPlex _spx; // LP solver
    soplex::DVectorRational _worker; // Temporary dense rational vector.
    bool _correctAlways; // Whether to always correct, regardless of feasibility.
  };

} /* namespace ipo */

#endif /* IPO_MIXED_INTEGER_PROGRAM_H_ */
