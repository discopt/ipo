#pragma once

#include <vector>

#include <ipo/config.hpp>
#include <ipo/export.hpp>

namespace ipo
{

  enum LPSense
  {
    MINIMIZE = -1,
    MAXIMIZE = 1
  };

  enum LPStatus
  {
    OPTIMAL = 1,          ///< LP is solved to optimality.
    UNBOUNDED = 2,        ///< LP is unbounded.
    INFEASIBLE = 3,       ///< LP is infeasible.
    TIME_LIMIT = 4,       ///< Aborted due to time limit.
    CYCLING = 5,          ///< Aborted due to detection of cycling.
    ITERATION_LIMIT = 6,  ///< Aborted due to iteration limit.
    NUMERICS = 7,         ///< Aborted due to numerical difficulties.
  };

  std::ostream& operator<<(std::ostream& stream, LPStatus status);

  typedef std::size_t LPRow;
  typedef std::size_t LPColumn;

  template <typename NumberType, bool RowsRemovable = false, bool ColumnsRemovable = false>
  class LP
  {
  public:
    typedef NumberType Number;

    static const Number& plusInfinity();

    static const Number& minusInfinity();

    LP();

    ~LP();

    std::size_t numRows() const;

    std::size_t numColumns() const;

    LPStatus status() const;

    Number getObjectiveValue() const;

    bool hasPrimalSolution() const;

    Number getPrimalValue(LPColumn column) const;

    std::vector<Number> getPrimalSolution(const std::vector<LPColumn>& columns = std::vector<LPColumn>()) const;

    bool hasDualSolution() const;

    Number getDualValue(LPRow row) const;

    void setSense(LPSense newSense);
    
    LPColumn addColumn(
      const Number& lowerBound,
      const Number& upperBound,
      const Number& objectiveCoefficient,
      const std::string& name = "",
      LPColumn unusedColumn = std::numeric_limits<std::size_t>::max());

    LPRow addRow(
      const Number& lhs,
      std::size_t numNonzeros,
      const LPColumn* nonzeroVariables,
      const Number* nonzeroCoefficients,
      const Number& rhs,
      const std::string& name = "",
      LPRow unusedRow = std::numeric_limits<std::size_t>::max());

    void update();

    void changeUpper(LPColumn column, const Number& newUpperBound);

    void changeLower(LPColumn column, const Number& newLowerBound);

    void changeBounds(LPColumn column, const Number& newLowerBound, const Number& newUpperBound);
    
    void changeObjective(LPColumn column, const Number& newObjectiveCoefficient);

    void write(const std::string& fileName) const;

    LPStatus solve();

  private:
    void* _implementation;
  };

} /* namespace ipo */
