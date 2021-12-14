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

    Space();

    /**
     * \brief Constructs a space with given variable names.
     */

    Space(std::vector<std::string>&& variableNames);

    /**
     * \brief Constructs a space with given variable names.
     */

    Space(const std::vector<std::string>& variableNames);

    /**
     * \brief Comparison for equality.
     */

    bool operator==(const Space& other) const;

    /**
     * \brief Comparison for non-equality.
     */

    inline bool operator!=(const Space& other) const
    {
      return !(*this == other);
    }

    /**
     * \brief Returns the dimension of the space.
     */

    inline std::size_t dimension() const
    {
      return _variableNames.size();
    }

    /**
     * \brief Returns the name of the given \p variableIndex.
     */

    inline const std::string& variable(std::size_t variableIndex) const
    {
      return _variableNames[variableIndex];
    }

    /**
     * \brief Adds a variable with the \p name.
     */

    void addVariable(const std::string& name);

    /**
     * \brief Prints a floating-point vector to \p stream.
     *
     * Prints a floating-point vector to given \p stream, using stored variable names. Nonzeros are
     * delimited by a comma.
     */

    void printVector(std::ostream& str, const sparse_vector<double>& vector,
      bool rounded = false) const;

    inline void printVector(std::ostream& str,
      const std::shared_ptr<sparse_vector<double>>& vector, bool rounded = false) const
    {
      return printVector(str, *vector, rounded);
    }

    std::string printVector(const sparse_vector<double>& vector, bool rounded = false) const;

    inline std::string printVector(const std::shared_ptr<sparse_vector<double>>& vector,
      bool rounded = false) const
    {
      return printVector(*vector, rounded);
    }

    void printLinearForm(std::ostream& str, const sparse_vector<double>& vector,
      bool rounded = false) const;

    void printLinearForm(std::ostream& str,
      const std::shared_ptr<sparse_vector<double>>& vector, bool rounded = false) const
    {
      return printLinearForm(str, *vector, rounded);
    }

    std::string printLinearForm(const sparse_vector<double>& vector, bool rounded = false) const;

    inline std::string printLinearForm(const std::shared_ptr<sparse_vector<double>>& vector,
      bool rounded = false) const
    {
      return printLinearForm(*vector, rounded);
    }

    void printConstraint(std::ostream& str, const Constraint<double>& constraint,
      bool rounded = false) const;

    std::string printConstraint(const Constraint<double>& constraint, bool rounded = false) const;

#if defined(IPO_RATIONAL)
    
    void printVector(std::ostream& str, const sparse_vector<rational>& vector,
      bool rounded = false) const;

    inline void printVector(std::ostream& str,
      const std::shared_ptr<sparse_vector<rational>>& vector, bool rounded = false) const
    {
      return printVector(str, *vector, rounded);
    }

    std::string printVector(const sparse_vector<rational>& vector, bool rounded = false) const;

    inline std::string printVector(const std::shared_ptr<sparse_vector<rational>>& vector,
      bool rounded = false) const
    {
      return printVector(*vector, rounded);
    }

    void printLinearForm(std::ostream& str, const sparse_vector<rational>& vector,
      bool rounded = false) const;

    inline void printLinearForm(std::ostream& str,
      const std::shared_ptr<sparse_vector<rational>>& vector, bool rounded = false) const
    {
      return printLinearForm(str, *vector, rounded);
    }

    std::string printLinearForm(const sparse_vector<rational>& vector, bool rounded = false) const;

    inline std::string printLinearForm(const std::shared_ptr<sparse_vector<rational>>& vector,
      bool rounded = false) const
    {
      return printLinearForm(*vector, rounded);
    }

    void printConstraint(std::ostream& str, const Constraint<rational>& constraint,
      bool rounded = false) const;

    std::string printConstraint(const Constraint<rational>& constraint,
      bool rounded = false) const;

#endif /* IPO_RATIONAL */

  private:
    std::vector<std::string> _variableNames;
  };

}
