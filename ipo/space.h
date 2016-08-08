#ifndef IPO_SPACE_H_
#define IPO_SPACE_H_

#include <vector>
#include <string>

#include "common.h"
#include "rational.h"
#include "vectors.h"

namespace ipo {

  /**
   * \brief Defines the ambient space for oracles.
   *
   * Defines the ambient space for one or more oracle(s), including the names of variables.
   * Contains relevant methods for printing linear forms, vectors, inequalities and equations.
   */

  class Space
  {
  public:
    /**
     * \brief Constructs a space without variables.
     *
     * Constructs the 0-dimensional space. Use \ref addVariable() to increase its dimension.
     */

    Space();

    /**
     * \brief Constructs a space with given \c variables.
     *
     * Constructs a space with given \c variables.
     */

    Space(const std::vector<std::string>& variables);

    /**
     * \brief Copy constructor.
     *
     * Copy constructor.
     */

    Space(const Space& other);

    /**
     * \brief Destructor.
     *
     */

    ~Space();

    /**
     * \brief Adds a variable with given \c name.
     *
     * Adds a variable with given \c name. Does only ensure uniqueness of this variable name
     * in debug mode.
     */

    void addVariable(const std::string& name);

    /**
     * \brief Returns the dimension of the space.
     *
     * Returns the dimension of the space, i.e., the number of variables.
     */

    inline std::size_t dimension() const
    {
      return _variables.size();
    }

    /**
     * \brief Returns the name of the variable indexed by \c var.
     *
     * Returns the name of the variable indexed by \c var.
     */

    inline const std::string& operator[](std::size_t var) const
    {
      return _variables[var];
    }

    /**
     * \brief Returns the name of the variable indexed by \c var.
     *
     * Returns the name of the variable indexed by \c var.
     */

    inline const std::string& variable(std::size_t var) const
    {
      return _variables[var];
    }

    /**
     * \brief Prints the linear form with these \c coefficients to \c stream.
     *
     * Prints the linear form with theses \c coefficients to \c stream using the variable names.
     * Does not emit a newline character.
     */

    void printLinearForm(std::ostream& stream, const soplex::SVectorRational* coefficients) const;

    /**
     * \brief Prints the linear form with these \c coefficients to \c stream.
     *
     * Prints the linear form with theses \c coefficients to \c stream using the variable names.
     * Does not emit a newline character.
     */

    void printLinearForm(std::ostream& stream, const SparseVector& coefficients) const;

    /**
     * \brief Prints the linear form with these \c coefficients to \c stream.
     *
     * Prints the linear form with theses \c coefficients to \c stream using the variable names.
     * Does not emit a newline character.
     */

    void printLinearForm(std::ostream& stream, const DenseVector* coefficients) const;

    /**
     * \brief Prints the \c row (inequality / equation) to \c stream.
     *
     * Prints the \c row (inequality / equation) to \c stream using the variable names. Does not
     * emit a newline character.
     */

    void printRow(std::ostream& stream, const soplex::LPRowRational& row) const;

    /**
     * \brief Prints the row (inequality / equation) \c index of the set of \c rows to \c stream.
     *
     * Prints the row (inequality / equation) \c index of the set of \c rows to \c stream using the
     * variable names. Does not emit a newline character.
     */

    void printRow(std::ostream& stream, const soplex::LPRowSetRational& rows, std::size_t index)
      const;

    /**
     * \brief Prints all \c rows (inequality / equation) to \c stream.
     *
     * Prints all \c rows (inequality / equation) to \c stream using the variable names. Emits a
     * newline character after each row.
     */

    void printRows(std::ostream& stream, const soplex::LPRowSetRational& rows) const;

    /**
     * \brief Prints \c vector to \c stream.
     *
     * Prints \c vector to \c stream using the variable names of the oracle.
     * Does not emit a newline character.
     */

    void printVector(std::ostream& stream, const soplex::SVectorRational* vector) const;

    /**
     * \brief Prints \c vector to \c stream.
     *
     * Prints \c vector to \c stream using the variable names of the oracle.
     * Does not emit a newline character.
     */

    void printVector(std::ostream& stream, const SparseVector& vector) const;

    /**
     * \brief Returns \c true iff spaces are equal.
     *
     * Returns \c true iff spaces have the same variables in the same order.
     */

    bool operator==(const Space& other) const;

    /**
     * \brief Returns \c true iff spaces are not equal.
     *
     * Returns \c true iff spaces are not equal (\sa operator==()).
     */

    inline bool operator!=(const Space& other) const
    {
      return !(*this == other);
    }

  protected:
    /*
     * \brief Implementation of row printing.
     *
     * Implementation of row printing.
     */

    void printRow(std::ostream& stream, const soplex::Rational* lhs, const soplex::Rational* rhs,
      const soplex::SVectorRational& vector) const;

    std::vector<std::string> _variables; // Variable names.
  };

}


#endif /* IPO_SPACE_H_ */
