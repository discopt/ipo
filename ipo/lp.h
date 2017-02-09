#ifndef IPO_LP_H_
#define IPO_LP_H_

#include <memory>

#include "common.h"

#include "rational.h"
#include "linear_constraint.h"
#include "space.h"

struct Scip;
typedef struct Scip SCIP;

namespace ipo {

  /**
   * \brief A linear set, i.e., an explicitly given polyhedron.
   *
   * A linear set, i.e., an explicitly given polyhedron.
   */

  class LinearSet
  {
  public:
    /**
     * \brief Creates a linear set.
     *
     * Creates a linear set from given data.
     */

    LinearSet(const std::vector<Rational>& lowerBounds, const std::vector<Rational>& upperBounds,
      const std::vector<LinearConstraint>& rows, const std::vector<std::string>& variableNames = std::vector<std::string>(),
      const std::vector<std::string>& rowNames = std::vector<std::string>());

    /**
     * \brief Creates a linear set as a copy of another one.
     * 
     * Creates a linear set as a copy of another one.
     */
    
    LinearSet(const LinearSet& other);

#ifdef IPO_WITH_SCIP

    /**
     * \brief Creates a linear set from a SCIP instance.
     * 
     * Creates a linear set from a SCIP instance.
     */

    LinearSet(SCIP* scip);

    /**
     * \brief Creates a linear set from a SCIP instance read from \p fileName.
     * 
     * Creates a linear set from a SCIP instance read from \p fileName.
     */

    LinearSet(const std::string& fileName);

    /**
     * \brief Constructs the linear set from a SCIP instance.
     * 
     * Constructs the linear set from a SCIP instance.
     */

    void constructLinearSet(SCIP* scip);

    /**
     * \brief Reads \p fileName into a SCIP instance.
     * 
     * Reads \p fileName into a SCIP instance.
     */

    SCIP* readSCIP(const std::string& fileName);

    /**
     * \brief Frees the SCIP instance \p scip.
     * 
     * Frees the SCIP instance \p scip.
     */

    void freeSCIP(SCIP* scip);

#endif /* IPO_WITH_SCIP */
    

    /**
     * \brief Destructor.
     */

    virtual ~LinearSet();

    /**
     * \brief Returns the ambient space.
     *
     * Returns a const-reference to the ambient space.
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

    inline std::size_t numVariables() const
    {
      return _space.dimension();
    }

    /**
     * \brief Returns the number of row constraints.
     *
     * Returns the number of row constraints.
     */

    inline std::size_t numRows() const
    {
      return _rowNames.size();
    }

    /**
     * \brief Returns the lower bound of the \p variableIndex variable.
     * 
     * Returns the lower bound of the \p variableIndex variable.
     */

    inline const Rational& lowerBound(std::size_t variableIndex) const
    {
      return _solver.lowerRational(variableIndex);
    }

    /**
     * \brief Returns the lower bound of the \p variableIndex variable as a LinearConstraint.
     * 
     * Returns the lower bound of the \p variableIndex variable as a LinearConstraint.
     */

    LinearConstraint lowerBoundConstraint(std::size_t variableIndex) const;

    /**
     * \brief Returns the upper bound of the \p variableIndex variable.
     * 
     * Returns the upper bound of the \p variableIndex variable.
     */

    inline const Rational& upperBound(std::size_t variableIndex) const
    {
      return _solver.upperRational(variableIndex);
    }

    /**
     * \brief Returns the upper bound of the \p variableIndex variable as a LinearConstraint.
     * 
     * Returns the upperbound of the \p variableIndex variable as a LinearConstraint.
     */

    LinearConstraint upperBoundConstraint(std::size_t variableIndex) const;

    /**
     * \brief Returns the name of the \p variableIndex variable.
     * 
     * Returns the name of the \p variableIndex variable.
     */

    inline const std::string& variableName(std::size_t variableIndex) const
    {
      return _space[variableIndex];
    }

    /**
     * \brief Returns the name of the \p rowIndex row.
     * 
     * Returns the name of the \p rowIndex row.
     */

    inline const std::string& rowName(std::size_t rowIndex) const
    {
      return _rowNames[rowIndex];
    }

    /**
     * \brief Returns the row index by \p rowIndex.
     * 
     * Returns the row index by \p rowIndex.
     */

    LinearConstraint rowConstraint(std::size_t rowIndex) const;

    /**
     * \brief Returns the row index by \p rowIndex.
     * 
     * Returns the row index by \p rowIndex.
     */

    void getConstraints(std::vector<LinearConstraint>& constraints, bool excludeEquations = false, 
      bool includeBounds = true, bool includeRows = true) const;

    /**
     * \brief Changes the lower bound of the specified variable.
     * 
     * Changes the lower bound of the specified variable.
     */

    void changeLowerBound(std::size_t variableIndex, const Rational& newLower);

    /**
     * \brief Changes the upper bound of the specified variable.
     * 
     * Changes the upper bound of the specified variable.
     */

    void changeUpperBound(std::size_t variableIndex, const Rational& newUpper);

    /**
     * \brief Changes the lower and upper bounds of the specified variable.
     * 
     * Changes the lower and upper bounds of the specified variable.
     */

    void changeBounds(std::size_t variableIndex, const Rational& newLower, const Rational& newUpper);

    /**
     * \brief Adds the given linear constraint.
     * 
     * Adds the given linear constraint. Returns the index of the new row.
     */

    std::size_t addConstraint(const LinearConstraint& constraint, const std::string& rowName = "");

    /**
     * \brief Removes a single row constraint.
     * 
     * Removes the specified row constraint.
     */

    void removeConstraint(std::size_t rowIndex);

    /**
     * \brief Removes the last row constraints until \p newNumRows are left.
     * 
     * Removes the last row constraints until \p newNumRows are left.
     */

    void removeLastConstraints(std::size_t newNumRows);

  protected:
    Space _space; // Space with variable names.
    soplex::SoPlex _solver;
    std::vector<std::string> _rowNames; // Row constraints' names.
  };

  class MixedIntegerLinearSet : public LinearSet
  {
  public:
    /**
     * \brief Creates a mixed-integer linear set.
     *
     * Creates a mixed-integer linear set.
     */

    MixedIntegerLinearSet(const std::vector<bool>& integrality, const std::vector<Rational>& lowerBounds,
      const std::vector<Rational>& upperBounds, const std::vector<LinearConstraint>& rows,
      const std::vector<std::string>& variableNames = std::vector<std::string>(), 
      const std::vector<std::string>& rowNames = std::vector<std::string>());

    /**
     * \brief Creates a mixed-integer linear set as a copy of another mixed-integer linear set.
     *
     * Creates a mixed-integer linear set as a copy of another mixed-integer linear set.
     */

    MixedIntegerLinearSet(const MixedIntegerLinearSet& other);

    /**
     * \brief Creates a mixed-integer linear set from a \ref LinearSet.
     *
     * Creates a mixed-integer linear set from a \ref LinearSet.
     */

    MixedIntegerLinearSet(const LinearSet& linearSet);

#ifdef IPO_WITH_SCIP

    /**
     * \brief Creates a mixed-integer linear set from a SCIP instance.
     * 
     * Creates a mixed-integer linear set from a SCIP instance.
     */

    MixedIntegerLinearSet(SCIP* scip);

    /**
     * \brief Creates a mixed-integer linear set from a SCIP instance read from \p fileName.
     * 
     * Creates a mixed-integer linear set from a SCIP instance read from \p fileName.
     */

    MixedIntegerLinearSet(const std::string& fileName);

    /**
     * \brief Reads integrality information from \p scip instance.
     * 
     * Reads integrality information from \p scip instance.
     */

    void constructIntegrality(SCIP *scip);

#endif /* IPO_WITH_SCIP */
    
    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    virtual ~MixedIntegerLinearSet();

    /**
     * \brief Returns if the \p variableIndex variable is integer.
     * 
     * Returns if the \p variableIndex variable is integer.
     */

    inline bool isIntegral(std::size_t variableIndex) const
    {
      assert(variableIndex < numVariables());
      return _integrality[variableIndex];
    }

  protected:
    std::vector<bool> _integrality;
  };

  class LinearProgram: public LinearSet
  {
  public:
    enum Result
    {
      INFEASIBLE = -1,
      OPTIMAL = 0,
      UNBOUNDED = 1
    };

    /**
     * \brief Creates a linear program.
     *
     * Creates a linear program.
     */

    LinearProgram(const std::vector<Rational>& objective, const std::vector<Rational>& lowerBounds,
      const std::vector<Rational>& upperBounds, const std::vector<LinearConstraint>& rows,
      const std::vector<std::string>& variableNames = std::vector<std::string>(),
      const std::vector<std::string >& rowNames = std::vector<std::string>());

    /**
     * \brief Creates a linear program as a copy of another linear program.
     *
     * Creates a linear program as a copy of another linear program.
     */

    LinearProgram(const LinearProgram& linearProgram);

    /**
     * \brief Creates a linear program from a \ref LinearSet.
     *
     * Creates a linear program from a \ref LinearSet.
     */

    LinearProgram(const LinearSet& linearSet);

#ifdef IPO_WITH_SCIP

    /**
     * \brief Creates a linear program from a SCIP instance.
     * 
     * Creates a linear program from a SCIP instance.
     */

    LinearProgram(SCIP* scip);
    
    /**
     * \brief Creates a linear program from a SCIP instance read from \p fileName.
     * 
     * Creates a linear program from a SCIP instance read from \p fileName.
     */

    LinearProgram(const std::string& fileName);

    /**
     * \brief Reads integrality information from \p scip instance.
     * 
     * Reads integrality information from \p scip instance.
     */

    void constructObjective(SCIP *scip);

#endif /* IPO_WITH_SCIP */

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    virtual ~LinearProgram();

    /**
     * \brief Replaces the current objective vector by \p objective.
     *
     * Replaces the current objective vector by \p objective.
     */

    void changeObjective(const std::vector<Rational>& objective);

    /**
     * \brief Solves the LP.
     *
     * Solves the LP and returns the result. If the LP is unbounded, \p vector contains an unbounded ray,
     * and otherwise, if feasible, contains an optimum solution.
     */

    Result solve(Vector& vector, Rational& objectiveValue);
  };

} /* namespace ipo */

#endif /* IPO_LP_H_ */
