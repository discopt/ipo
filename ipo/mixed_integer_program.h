#ifndef IPO_MIXED_INTEGER_PROGRAM_H_
#define IPO_MIXED_INTEGER_PROGRAM_H_

#include "ipo.h"

#include <string>
#include <vector>

#ifdef WITH_SCIP
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include "scip_oracles.h"
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include "scip_oracles.h"
#endif
#endif

#include "oracles.h"

namespace ipo {

  /**
   * \brief A mixed-integer program.
   *
   * A mixed-integer program.
   * It is typically used to create oracles for optimizing over it,
   * e.g., using \ref SCIPOracle and \ref ExactSCIPOracle.
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

    MixedIntegerProgram(Space& space, SCIP* scip);
#endif

    /**
     * \brief Destructor.
     */

    virtual ~MixedIntegerProgram();

    /**
     * \brief Returns the space.
     *
     * Returns a const reference to the space.
     */

    inline const Space& space() const
    {
      return _space;
    }

    /**
     * \brief Returns the number of columns.
     *
     * Returns the number of columns.
     */

    inline std::size_t numColumns() const
    {
      return _space.dimension();
    }

    /**
     * \brief Returns the number of rows.
     *
     * Returns the number of rows.
     */

    inline std::size_t numRows() const
    {
      return _rowNames.size();
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
     * \brief Returns the name of the given row.
     *
     * Returns the name of the given \c row.
     */

    inline const std::string& rowName(std::size_t row) const
    {
      return _rowNames[row];
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

    bool checkPointRows(const soplex::SVectorRational* point);

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

    bool checkRayRows(const soplex::SVectorRational* ray);

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

    void setFace(Face* newFace);

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

    void getFixedVariableEquations(soplex::LPRowSetRational& rows,
      std::vector<std::string>* names = NULL);

  protected:
    const Space& _space; // Space with column names.
    soplex::LPColSetRational _columns; // Columns
    soplex::LPRowSetRational _rows; // Rows
    std::vector<std::string> _rowNames; // Row names
    std::vector<bool> _integrality; // Integrality constraints
    soplex::DVectorRational _worker; // Temporary dense rational vector.
    Face* _currentFace; // Currently active face.
  };

  /**
   * \brief An oracle that postprocesses floating-point solutions from an approximate oracle.
   *
   * An oracle that can postprocess floating-point solutions
   * of a mixed-integer program to get correct rational ones.
   */

  class MixedIntegerProgramCorrectorOracle: public OracleBase
  {
  public:

    /**
     * \brief Constructs an oracle that postprocesses floating-point solutions to make them exact.
     *
     * Constructs an oracle that postprocesses floating-point solutions to make them exact for the
     * given \c mip. It calls the \c approximateOracle and then solves an LP to recompute the
     * continuous part.
     *
     * \param name              Name of the oracle.
     * \param mip               Associated mixed-integer program.
     * \param approximateOracle Approximate oracle whose solutions are postprocessed.
     */

    MixedIntegerProgramCorrectorOracle(const std::string& name, MixedIntegerProgram& mip,
      OracleBase* approximateOracle);

    /**
     * \brief Constructs a heuristic oracle that postprocesses floating-point solutions.
     *
     * Constructs an oracle that postprocesses floating-point solutions to make them exact for the
     * given \c mip. It calls the \c approximateOracle and then solves an LP to recompute the
     * continuous part.
     *
     * \param name              Name of the oracle.
     * \param mip               Associated mixed-integer program.
     * \param approximateOracle Approximate oracle whose solutions are postprocessed.
     * \param nextOracle        Next associated oracle.
     */

    MixedIntegerProgramCorrectorOracle(const std::string& name, MixedIntegerProgram& mip,
      OracleBase* approximateOracle, OracleBase* nextOracle);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~MixedIntegerProgramCorrectorOracle();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation adds a corresponding equation to the LP that is used for the
     * postprocessing.
     */

    virtual void setFace(Face* newFace = NULL);

    /**
     * \brief Runs the oracle to maximize the dense rational \p objective.
     *
     * Runs the optimization oracle to maximize the given dense rational \p objective
     * over the current face \f$ F \f$ (see setFace()) and returns \p result.
     * If \p maxHeuristic is less than thisHeuristic() or if the objective value
     * requested by \p objectiveBound is not exceeded, then the call must be forwarded to the
     * next oracle.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param maxHeuristic   Requested maximum heuristic level.
     * \param minHeuristic   Requested minimum heuristic level.
     *
     * This implementation forwards the call to the oracle and then solves an LP for every
     * solution that is returned.
     */

    virtual void maximize(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound = ObjectiveBound(),
      std::size_t maxHeuristic = std::numeric_limits<std::size_t>::max(),
      std::size_t minHeuristic = 0);

  protected:

    /**
     * \brief Initializes the LP and solver.
     *
     * Initializes the LP and solver.
     */

    void initializeLP();

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

    soplex::DSVectorRational* correctDirection(const soplex::SVectorRational* direction,
      const soplex::VectorRational& objective);

  protected:

    MixedIntegerProgram& _mip; // Associated mixed-integer program.
    OracleBase* _approximateOracle; // Approximate oracle.
    soplex::SoPlex _spx; // LP solver
    soplex::DVectorRational _denseVector; // Temporary dense rational vector.
  };

} /* namespace ipo */

#endif /* IPO_MIXED_INTEGER_PROGRAM_H_ */
