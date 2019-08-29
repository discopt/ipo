#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <vector>
#include <string>

#ifdef IPO_WITH_GMP
#include <gmpxx.h>
#endif

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
     * \brief Returns the name of the given \p variable.
     */

    IPO_EXPORT
    inline const std::string& operator[](std::size_t variable) const
    {
      return _variableNames[variable];
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
    void printVector(std::ostream& stream, std::size_t numNonzeros,
        std::size_t* nonzeroCoordinates, double* nonzeroValues) const;
    
    /**
     * \brief Prints a floating-point vector into a string.
     *
     * Prints a floating-point vector into a string, using stored variable names. Nonzeros are 
     * delimited by a comma.
     */

    IPO_EXPORT
    std::string vectorToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
      double* nonzeroValues) const;

    /**
     * \brief Prints a floating-point linear form to \p stream.
     *
     * Prints a floating-point linear form to given \p stream, using stored variable names.
     * Nonzeros are delimited by a comma.
     */

    IPO_EXPORT
    void printLinearForm(std::ostream& stream, std::size_t numNonzeros,
      std::size_t* nonzeroColumns, double* nonzeroValues) const;

    /**
     * \brief Prints a floating-point linear form into a string.
     *
     * Prints a floating-point linear form into a string, using stored variable names. Nonzeros
     * are delimited by a comma.
     */

    IPO_EXPORT
    std::string linearFormToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
      double* nonzeroValues) const;

#if defined(IPO_WITH_GMP)
    /**
     * \brief Prints a rational vector to \p stream.
     *
     * Prints a rational vector to given \p stream, using stored variable names. Nonzeros are 
     * delimited by a comma.
     */

    IPO_EXPORT
    void printVector(std::ostream& stream, std::size_t numNonzeros,
      std::size_t* nonzeroCoordinates, mpq_class* nonzeroValues) const;

    /**
     * \brief Prints a rational vector into a string.
     *
     * Prints a rational vector into a string, using stored variable names. Nonzeros are 
     * delimited by a comma.
     */

    IPO_EXPORT
    std::string vectorToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
      mpq_class* nonzeroValues) const;

    /**
     * \brief Prints a rational linear form to \p stream.
     *
     * Prints a rational linear form to given \p stream, using stored variable names.
     * Nonzeros are delimited by a comma.
     */

    IPO_EXPORT
    void printLinearForm(std::ostream& stream, std::size_t numNonzeros,
      std::size_t* nonzeroCoordinates, mpq_class* nonzeroValues) const;

    /**
     * \brief Prints a rational linear form into a string.
     *
     * Prints a rational linear form into a string, using stored variable names. Nonzeros
     * are delimited by a comma.
     */

    IPO_EXPORT
    std::string linearFormToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
      mpq_class* nonzeroValues) const;

#endif /* IPO_WITH_GMP */

  private:
    std::vector<std::string> _variableNames;
  };

}
