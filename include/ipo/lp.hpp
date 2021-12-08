#pragma once

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
    UNKNOWN = 0,
    OPTIMAL = 1,
    INFEASIBLE = 2,
    UNBOUNDED = 3,
    TIMEOUT = 4
  };

  typedef std::size_t LPRow;
  typedef std::size_t LPColumn;

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  class LP
  {
  public:
    typedef NumberType Number;

    LP();

    ~LP();

    std::size_t numRows() const;

    std::size_t numColumns() const;

    LPStatus status() const;

    double getObjectiveValue() const;

    bool hasPrimalSolution() const;

    double getPrimalValue(LPColumn column) const;

    bool hasDualSolution() const;

    double getDualValue(LPRow row) const;

    void setSense(LPSense newSense);
    
    LPColumn addColumn(
      double lowerBound,
      double upperBound,
      double objectiveCoefficient,
      const std::string& name = "",
      LPColumn unusedColumn = std::numeric_limits<std::size_t>::max());

    LPRow addRow(
      double lhs,
      std::size_t numNonzeros,
      const LPColumn* nonzeroVariables,
      const double* nonzeroCoefficients,
      double rhs,
      const std::string& name = "",
      LPRow unusedRow = std::numeric_limits<std::size_t>::max());

    void update();

    void changeObjective(LPColumn column, double newObjectiveCoefficient);

    LPStatus solve();

  private:
    void* _implementation;
  };

} /* namespace ipo */
