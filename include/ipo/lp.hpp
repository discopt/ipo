#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/arithmetic.hpp>

#include <vector>
#include <iostream>

#if defined(IPO_WITH_DOUBLE_LP) || defined(IPO_WITH_RATIONAL_LP)

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

  enum LPBasisStatus
  {
    BASIC = 0,
    NONBASIC_LOWER = 1,
    NONBASIC_UPPER = 2,
    NONBASIC_ZERO = 3
  };

  /**
   * \brief Writes LP status to \p stream.
   */

  std::ostream& operator<<(std::ostream& stream, LPStatus status);

  template <typename NumberType>
  class LP;

  template <typename NumberType>
  class LPImplementation;

  /**
   * \brief Unique identifier for a row or a column.
   */

  class LPKey
  {
  public:
    /**
     * \brief Creates an invalid key.
     */

    LPKey();

    /**
     * \brief Copy constructor.
     */

    LPKey(const LPKey& other);

    /**
     * \brief Assignment operator.
     */

    LPKey& operator=(const LPKey& other);

    /**
     * \brief Returns whether this key is valid.
     */

    bool isValid() const;

  private:
    /**
     * \brief Creates a key with a specific internal id.
     */

    LPKey(long id);

    /**
     * \brief Retrieves the internal id.
     */

    long id() const;

    friend class LP<double>;
    friend class LPImplementation<double>;

#if defined(IPO_WITH_RATIONAL_LP)
    friend class LP<rational>;
    friend class LPImplementation<rational>;
#endif /* IPO_WITH_RATIONAL_LP */

    /// Internal id.
    long _id;
  };

  template <typename NumberType>
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

    bool hasPrimalRay() const;

    Number getPrimalValue(int column) const;

    std::vector<Number> getPrimalSolution(const std::vector<int>& columns = std::vector<int>()) const;

    std::vector<Number> getPrimalRay(const std::vector<int>& columns = std::vector<int>()) const;

    bool hasDualSolution() const;

    Number getDualValue(int row) const;

    void setSense(LPSense newSense);

    LPKey addColumn(
      const Number& lowerBound,
      const Number& upperBound,
      const Number& objectiveCoefficient,
      const std::string& name = "");

    LPKey addRow(const Number& lhs, std::size_t numNonzeros, const std::size_t* nonzeroVariables,
      const Number* nonzeroCoefficients, const Number& rhs, const std::string& name = "");

    void update();

    void changeUpper(int column, const Number& newUpperBound);

    void changeLower(int column, const Number& newLowerBound);

    void changeBounds(int column, const Number& newLowerBound, const Number& newUpperBound);
    
    void changeObjective(int column, const Number& newObjectiveCoefficient);

    void changeRow(LPKey rowKey, const Number& lhs, std::size_t numNonzeros, const int* nonzeroVariables,
      const Number* nonzeroCoefficients, const Number& rhs);

    void write(const std::string& fileName) const;

    LPStatus solve(bool extreme = false);

    double getSolveTime() const;

  private:

    void* _implementation;
  };

} /* namespace ipo */

#endif /* IPO_WITH_DOUBLE_LP || IPO_WITH_RATIONAL_LP */
