#ifndef IPO_SPACE_H_
#define IPO_SPACE_H_

#include <vector>
#include <string>

#include "common.h"
#include "rational.h"
#include "vectors.h"
#include "linear_constraint.h"

namespace ipo {

  class Space;

  /**
   * \brief Storage class of a \ref Space.
   * 
   * Storage class of a \ref Space. Stores the variable names and a reference counter.
   */

  class SpaceData
  {
  public:
    /**
     * \brief Constructs without variables. Use addVariable() to add them.
     * 
     * Constructs without variables. Use addVariable() to add them.
     */

    SpaceData();

    /**
     * \brief Constructs with given \p variableNames.
     * 
     * Constructs with given \p variableNames. Use addVariable() to add further ones.
     */

    SpaceData(const std::vector<std::string>& variableNames);

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    ~SpaceData();

    /**
     * \brief Returns true iff all variables coincide.
     * 
     * Returns true iff all variables coincide.
     */

    bool operator==(const SpaceData& other) const;

    /**
     * \brief Returns true iff not all variables coincide.
     * 
     * Returns true iff not all variables coincide.
     */

    inline bool operator!=(const SpaceData& other) const
    {
      return !(*this == other);
    }

    /**
     * \brief Returns the dimension of the space.
     * 
     * Returns the dimension of the space.
     */

    inline std::size_t dimension() const
    {
      return _variableNames.size();
    }

    /**
     * \brief Returns the name of the given \p variable.
     * 
     * Returns a const-reference to the name of the given \p variable.
     */

    inline const std::string& variableName(std::size_t variable) const
    {
      return _variableNames[variable];
    }

    /**
     * \brief Adds a variable with given \p name.
     * 
     * Adds a variable with given \p name to the space.
     */
    
    void addVariable(const std::string& name);

    /**
     * \brief Increases the usage counter by 1.
     * 
     * Increases the usage counter by 1.
     */

    inline void markUsed()
    {
      _usage++;
    }

    /**
     * \brief Decreases the usage counter by 1.
     * 
     * Decreases the usage counter by 1 and frees it if the latter reached 0.
     */
    
    void unmarkUsed();

    friend class Space;

  protected:
    std::vector<std::string> _variableNames; // The variable names.
    std::size_t _usage; // Reference counter.
  };

  /**
   * \brief Defines the ambient space for polyhedra.
   *
   * Defines the ambient space for polyhedra defined by more oracles and stores the names of variables. Contains relevant methods 
   * for printing vectors, linear forms and linear constraints. Implementation is by reference counting, i.e., copying is cheap
   * and ownership is managed automatically.
   */

  class Space
  {
  public:
    /**
     * \brief Constructs the empty space.
     * 
     * Constructs the empty space.
     */

    Space();

    /**
     * \brief Constructs from \ref SpaceData object.
     * 
     * Constructs from \ref SpaceData object.
     */

    inline Space(SpaceData* data)
      : _data(data)
    {
      _data->_usage++;
    }

    /**
     * \brief Copy constructor.
     * 
     * Copy constructor.
     */

    inline Space(const Space& other)
      : _data(other._data)
    {
      _data->markUsed();
    }

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    ~Space()
    {
      _data->unmarkUsed();
    }

    /**
     * \brief Assignment operator.
     * 
     * Assignment operator.
     */

    inline Space& operator=(const Space& other)
    {
      if (_data != other._data)
      {
        _data->unmarkUsed();
        _data = other._data;
        _data->markUsed();
      }
      return *this;
    }

    /**
     * \brief Returns true if spaces are equal.
     * 
     * Returns true if spaces are equal, including variable names.
     */

    inline bool operator==(const Space& other) const
    {
      if (_data == other._data)
        return true;
      return *_data == *other._data;
    }
    
    /**
     * \brief Returns true if spaces are not equal.
     * 
     * Returns true if spaces are not equal, also considering variable names.
     */

    inline bool operator!=(const Space& other) const
    {
      return !(*this == other);
    }

    /**
     * \brief Returns the dimension of the space.
     * 
     * Returns the dimension of the space.
     */

    inline std::size_t dimension() const
    {
      return _data->dimension();
    }

    /**
     * \brief Returns the name of the given \p variable.
     * 
     * Returns a const-reference to the name of the given \p variable.
     */

    inline const std::string& operator[](std::size_t variable) const
    {
      return _data->variableName(variable);
    }

    /**
     * \brief Prints given \p vector to given \p stream.
     * 
     * Prints given \p vector to given \p stream, using stored variable names. Nonzeros are delimited by a comma.
     */

    void printVector(std::ostream& stream, const Vector& vector) const;

    /**
     * \brief Prints given \p linearForm to given \p stream.
     * 
     * Prints given \p linearForm to given \p stream, using stored variable names. Nonzeros are delimited by the sign of the next
     * entry.
     */

    void printLinearForm(std::ostream& stream, const Vector& linearForm) const;

    /**
     * \brief Prints given linear \p constraint to given \p stream.
     * 
     * Prints given linear \p constraint to given \p stream, using stored variable names. Nonzeros are delimited by the sign of 
     * the next entry.
     */

    void printLinearConstraint(std::ostream& stream, const LinearConstraint& constraint) const;

  protected:
    SpaceData* _data; // Pointer to data.
  };

}


#endif /* IPO_SPACE_H_ */
