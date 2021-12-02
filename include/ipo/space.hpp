#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <vector>
#include <string>

#include <ipo/sparse_vector.hpp>
#include <ipo/constraint.hpp>

namespace ipo
{
  /**
   * \brief Ambient space of certain dimension.
   *
   * Ambient space of certain dimension. Manages variable names.
   */

  class Space
  {
  public:
    /**
     * \brief Constructs a space of dimension 0.
     *
     * Constructs a space of dimension 0. Use \ref addVariable to add variables.
     */

    IPO_EXPORT
    Space();

    /**
     * \brief Constructs a space with given variable names.
     */

    IPO_EXPORT
    Space(std::vector<std::string>&& variableNames);

    /**
     * \brief Constructs a space with given variable names.
     */

    IPO_EXPORT
    Space(const std::vector<std::string>& variableNames);

    /**
     * \brief Comparison for equality.
     */

    IPO_EXPORT
    bool operator==(const Space& other) const;

    /**
     * \brief Comparison for non-equality.
     */

    IPO_EXPORT
    inline bool operator!=(const Space& other) const
    {
      return !(*this == other);
    }

    /**
     * \brief Returns the dimension of the space.
     */

    IPO_EXPORT
    inline std::size_t dimension() const
    {
      return _variableNames.size();
    }

    /**
     * \brief Returns the name of the given \p variableIndex.
     */

    IPO_EXPORT
    inline const std::string& variable(std::size_t variableIndex) const
    {
      return _variableNames[variableIndex];
    }

    /**
     * \brief Adds a variable with the \p name.
     */

    IPO_EXPORT
    void addVariable(const std::string& name);

    /**
     * \brief Prints a floating-point vector to \p stream.
     *
     * Prints a floating-point vector to given \p stream, using stored variable names. Nonzeros are
     * delimited by a comma.
     */

    IPO_EXPORT
    void printVector(std::ostream& str, const sparse_vector<double>& vector,
      bool rounded = false) const;

    IPO_EXPORT
    inline void printVector(std::ostream& str,
      const std::shared_ptr<sparse_vector<double>>& vector, bool rounded = false) const
    {
      return printVector(str, *vector, rounded);
    }

    IPO_EXPORT
    std::string printVector(const sparse_vector<double>& vector, bool rounded = false) const;

    IPO_EXPORT
    inline std::string printVector(const std::shared_ptr<sparse_vector<double>>& vector,
      bool rounded = false) const
    {
      return printVector(*vector, rounded);
    }

    IPO_EXPORT
    void printLinearForm(std::ostream& str, const sparse_vector<double>& vector,
      bool rounded = false) const;

    IPO_EXPORT
    void printLinearForm(std::ostream& str,
      const std::shared_ptr<sparse_vector<double>>& vector, bool rounded = false) const
    {
      return printLinearForm(str, *vector, rounded);
    }

    IPO_EXPORT
    std::string printLinearForm(const sparse_vector<double>& vector, bool rounded = false) const;

    IPO_EXPORT
    inline std::string printLinearForm(const std::shared_ptr<sparse_vector<double>>& vector,
      bool rounded = false) const
    {
      return printLinearForm(*vector, rounded);
    }

    IPO_EXPORT
    void printConstraint(std::ostream& str, const Constraint<double>& constraint,
      bool rounded = false) const;

    IPO_EXPORT
    std::string printConstraint(const Constraint<double>& constraint, bool rounded = false) const;

#if defined(IPO_WITH_GMP)
    
    IPO_EXPORT
    void printVector(std::ostream& str, const sparse_vector<mpq_class>& vector,
      bool rounded = false) const;

    IPO_EXPORT
    inline void printVector(std::ostream& str,
      const std::shared_ptr<sparse_vector<mpq_class>>& vector, bool rounded = false) const
    {
      return printVector(str, *vector, rounded);
    }

    IPO_EXPORT
    std::string printVector(const sparse_vector<mpq_class>& vector, bool rounded = false) const;

    IPO_EXPORT
    inline std::string printVector(const std::shared_ptr<sparse_vector<mpq_class>>& vector,
      bool rounded = false) const
    {
      return printVector(*vector, rounded);
    }

    IPO_EXPORT
    void printLinearForm(std::ostream& str, const sparse_vector<mpq_class>& vector,
      bool rounded = false) const;

    IPO_EXPORT
    inline void printLinearForm(std::ostream& str,
      const std::shared_ptr<sparse_vector<mpq_class>>& vector, bool rounded = false) const
    {
      return printLinearForm(str, *vector, rounded);
    }

    IPO_EXPORT
    std::string printLinearForm(const sparse_vector<mpq_class>& vector, bool rounded = false) const;

    IPO_EXPORT
    inline std::string printLinearForm(const std::shared_ptr<sparse_vector<mpq_class>>& vector,
      bool rounded = false) const
    {
      return printLinearForm(*vector, rounded);
    }

    IPO_EXPORT
    void printConstraint(std::ostream& str, const Constraint<mpq_class>& constraint,
      bool rounded = false) const;

    IPO_EXPORT
    std::string printConstraint(const Constraint<mpq_class>& constraint,
      bool rounded = false) const;

#endif /* IPO_WITH_GMP */

  private:
    std::vector<std::string> _variableNames;
  };

}
