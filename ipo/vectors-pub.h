#ifndef IPO_VECTORS_H_
#define IPO_VECTORS_H_

#include <set>

#include "common.h"
#include "rational.h"

namespace ipo {

  class ReferenceCountedVector;
  class MutableVector;
  class Vector;

  /**
   * \brief Internal representation for IPO vectors.
   *
   * Internal representation for IPO vectors. Stores only the nonzeros of the vector in exact and approximate form. Has fields
   * for the mutable and immutable usage that can be used by \ref MutableVector and \ref Vector.
   */

  class VectorData
  {
  private:
    /**
     * \brief Nonzero data.
     *
     * Nonzero data.
     */

    struct Nonzero
    {
      std::size_t index; // Coordinate of this nonzero entry, starting from 0.
      Rational value; // Exact value of this nonzero entry.
      double approximation; // Approximate value of this nonzero entry.

      /**
       * \brief Creates a nonzero entry.
       *
       * Creates a nonzero entry.
       */

      Nonzero(std::size_t index, const Rational& value);

      /**
       * \brief Destructor.
       *
       * Destructor.
       */

      ~Nonzero();

      /**
       * \brief Compares two nonzeros, based on their index.
       *
       * Compares two nonzeros, based on their index. Returns true iff this index is smaller.
       */

      bool operator<(const Nonzero& other) const;
    };

  public:
    /**
     * \brief Constructs a zero vector with a default initial memory amount.
     *
     * Constructs a zero vector with a default initial memory amount.
     */

    VectorData();

    /**
     * \brief Constructs a zero vector with given initial memory amount.
     *
     * Constructs a zero vector with given initial memory amount.
     *
     * \param initialMemory Number of nonzeros that the user expects.
     */

    VectorData(std::size_t initialMemory);

    /**
     *\brief Constructs a vector from another one by copying.
     *
     * Constructs a vector from another one by copying.
     */

    VectorData(const VectorData& other);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    ~VectorData();

    /**
     * \brief Copys another vector to this one.
     *
     * Copys another vector to this one, overwriting the previous contents.
     */

    VectorData& operator=(const VectorData& other);

    /**
     * \brief Returns the number of nonzeros.
     *
     * Returns the number of nonzeros.
     */

    inline const std::size_t size() const
    {
      return _nonzeros.size();
    }

    /**
     * \brief Returns the index of the nonzero specified by \p position.
     *
     * Returns the index of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Index of the nonzero.
     */

    inline const std::size_t index(std::size_t position) const
    {
      assert(position < _nonzeros.size());
      return _nonzeros[position].index;
    }

    /**
     * \brief Returns the value of the nonzero specified by \p position.
     *
     * Returns a const-reference to the value of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Value of the nonzero.
     */

    inline const Rational& value(std::size_t position) const
    {
      assert(position < _nonzeros.size());
      return _nonzeros[position].value;
    }

    /**
     * \brief Returns the value of the nonzero specified by \p position.
     *
     * Returns a reference to the value of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Reference to the value of the nonzero.
     */

    inline Rational& value(std::size_t position)
    {
      assert(position < _nonzeros.size());
      return _nonzeros[position].value;
    }

    /**
     * \brief Returns the approximate value of the nonzero specified by \p position.
     *
     * Returns the approximate value of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Approximate value of the nonzero.
     */

    inline const double approximation(std::size_t position) const
    {
      assert(position < _nonzeros.size());
      return _nonzeros[position].approximation;
    }

    /**
     * \brief Tests for equality with \c other vector.
     *
     * Tests for equality with \c other vector.
     */

    bool operator==(const VectorData& other) const;

    /**
     * \brief Fast ordering comparison.
     *
     * Fast ordering comparison. First compares the number of nonzeros, then the indices and approximate values and finally the
     * exact values if all the previous properties were equal.
     */

    bool operator<(const VectorData& other) const;

    /**
     * \brief Adds a nonzero to the vector.
     *
     * Adds a nonzero to the vector. Resizes the memory if required.
     */

    void add(std::size_t index, const Rational& value);

    /**
     * \brief Sorts the nonzero entries by index.
     *
     * Sorts the nonzero entries by index.
     */

    void sort();

    /**
     * \brief Checks if the nonzero entries are sorted.
     *
     * Checks if the nonzero entries are sorted, returning true iff this is the case.
     */

    bool isSorted() const;

    /**
     * \brief Deletes this object if both usage counters are zero.
     *
     * Deletes this object if both usage counters are zero.
     */

    void deleteIfUnused();

    /**
     * \brief Checks if this vector is mutable.
     *
     * Checks if this vector is mutable, returning true iff the immutable-usage counter is zero.
     */

    inline bool isMutable() const
    {
      return _immutableUsage == 0;
    }

    friend class ReferenceCountedVector;
    friend class MutableVector;
    friend class Vector;

  private:
    std::vector<Nonzero> _nonzeros;
    std::size_t _mutableUsage;
    std::size_t _immutableUsage;
  };

  /**
   * \brief Reference-counted vector.
   *
   * Reference-counted vector. It manages a pointer to a \ref VectorData object.
   */

  class ReferenceCountedVector
  {
  public:
    /**
     * \brief Constructs the zero vector.
     *
     * Constructs the zero vector.
     */

    ReferenceCountedVector();

    /**
     * \brief Constructs a vector from the \p data.
     *
     * Constructs a vector from the \p data.
     */

    inline ReferenceCountedVector(VectorData* data)
      : _data(data)
    {
      assert(_data != NULL);
    }

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    inline ~ReferenceCountedVector()
    {
      _data->deleteIfUnused();
    }

    /**
     * \brief Returns the number of nonzeros.
     *
     * Returns the number of nonzeros.
     */

    inline std::size_t size() const
    {
      return _data->size();
    }

    /**
     * \brief Returns the index of the nonzero specified by \p position.
     *
     * Returns the index of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Index of the nonzero.
     */

    inline std::size_t index(std::size_t position) const
    {
      return _data->index(position);
    }

    /**
     * \brief Returns the value of the nonzero specified by \p position.
     *
     * Returns the value of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Value of the nonzero.
     */

    inline const Rational& value(std::size_t position) const
    {
      return _data->value(position);
    }

    /**
     * \brief Returns the approximate value of the nonzero specified by \p position.
     *
     * Returns the approximate value of the nonzero specified by \p position.
     *
     * \param position Specifies the the nonzero, starting with 0.
     * \returns        Approximate value of the nonzero.
     */

    inline double approximation(std::size_t position) const
    {
      return _data->approximation(position);
    }

    /**
     * \brief Tests for equality with \c other vector.
     *
     * Tests for equality with \c other vector.
     */

    inline bool operator==(const ReferenceCountedVector& other) const
    {
      if (_data == other._data)
        return true;

      return (*_data) == (*other._data);
    }

    /**
     * \brief Tests for inequality with \c other vector.
     *
     * Tests for inequality with \c other vector.
     */

    inline bool operator!=(const ReferenceCountedVector& other) const
    {
      return !(*this == other);
    }

    /**
     * \brief Fast ordering comparison.
     *
     * Fast ordering comparison. First compares the number of nonzeros, then the indices and approximate values and finally the
     * exact values if all the previous properties were equal.
     */

    inline bool operator<(const ReferenceCountedVector& other) const
    {
      if (_data == other._data)
        return false;

      return (*_data) < (*other._data);
    }

    /**
     * \brief Checks if the nonzero entries are sorted.
     *
     * Checks if the nonzero entries are sorted, returning true iff this is the case.
     */

    inline bool isSorted() const
    {
      return _data->isSorted();
    }

    friend MutableVector;
    friend Vector;

  protected:
    VectorData* _data;
  };

  /**
   * \brief Reference-counted mutable vector.
   *
   * Reference-counted mutable vector. It manages a pointer to a \ref VectorData object. The modification methods throw exceptions
   * if called when the \ref VectorData object is referenced by a \ref Vector object since the latter expects that it is not
   * changed.
   */

  class MutableVector : public ReferenceCountedVector
  {
  public:
    /**
     * \brief Constructs the zero vector.
     *
     * Constructs the zero vector.
     */

    inline MutableVector()
      : ReferenceCountedVector()
    {
      _data->_mutableUsage++;
    }

    /**
     * \brief Constructs a vector from the \p data.
     *
     * Constructs a vector from the \p data.
     */

    inline MutableVector(VectorData* data)
      : ReferenceCountedVector(data)
    {
      _data->_mutableUsage++;
    }

    /**
     * \brief Constructs a vector from \p other mutable vector without copying it.
     *
     * Constructs a vector from \p other mutable vector without copying it.
     */

    inline MutableVector(const MutableVector& other)
      : ReferenceCountedVector(other._data)
    {
      _data->_mutableUsage++;
    }

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    inline ~MutableVector()
    {
      assert(_data->_mutableUsage > 0);
      _data->_mutableUsage--;
    }

    /**
     * \brief Assignment from \c other mutable vector.
     *
     * Assignment from \c other mutable vector.
     */

    MutableVector& operator=(const MutableVector& other)
    {
      assert(_data->_mutableUsage > 0);
      _data->_mutableUsage--;
      _data->deleteIfUnused();

      _data = other._data;
      _data->_mutableUsage++;
    }

    /**
     * \brief Checks if this vector is mutable.
     *
     * Checks if this vector is mutable, returning true iff the immutable-usage counter is zero.
     */

    inline bool isMutable() const
    {
      return _data->isMutable();
    }

    /**
     * \brief Adds a nonzero to the vector.
     *
     * Adds a nonzero to the vector. Resizes the memory if required. Throws an exception if the underlying \ref VectorData object
     * is referenced by a \ref Vector since the latter expects it not to be changed.
     */

    inline void add(std::size_t index, const Rational& value)
    {
      if (!isMutable())
        throw std::runtime_error("MutableVector::add not allowed since underlying vector referenced by Vector.");

      _data->add(index, value);
    }

    /**
     * \brief Sorts the nonzero entries by index.
     *
     * Sorts the nonzero entries by index. Throws an exception if the underlying \ref VectorData object is referenced by
     * a \ref Vector since the latter expects it not to be changed.
     */

    inline void sort()
    {
      if (!isMutable())
        throw std::runtime_error("MutableVector::sort not allowed since underlying vector referenced by Vector.");

      _data->sort();
    }

    /**
     * \brief Swaps two mutable vector objects efficiently.
     *
     * Swaps two mutable vector objects efficiently.
     */

    void swap(MutableVector& other)
    {
      std::swap(_data, other._data);
    }
  };


  /**
   * \brief Reference-counted immutable vector.
   *
   * Reference-counted immutable vector. It manages a pointer to a \ref VectorData object.
   */

  class Vector : public ReferenceCountedVector
  {
  public:
    /**
     * \brief Constructs the zero vector.
     *
     * Constructs the zero vector.
     */

    inline Vector()
      : ReferenceCountedVector()
    {
      _data->_immutableUsage++;
    }

    /**
     * \brief Constructs a vector from the \p data.
     *
     * Constructs a vector from the \p data.
     */

    inline Vector(VectorData* data)
      : ReferenceCountedVector(data)
    {
      _data->_immutableUsage++;
    }

    /**
     * \brief Constructs a vector from \p other immutable vector without copying it.
     *
     * Constructs a vector from \p other immutable vector without copying it.
     */

    inline Vector(const Vector& other)
      : ReferenceCountedVector(other._data)
    {
      _data->_immutableUsage++;
    }

    /**
     * \brief Constructs a vector from mutable vector \p other without copying it.
     *
     * Constructs a vector from mutable vector \p other without copying it.
     */

    inline Vector(const MutableVector& other)
      : ReferenceCountedVector(other._data)
    {
      _data->_immutableUsage++;
    }

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    inline ~Vector()
    {
      assert(_data->_immutableUsage > 0);
      _data->_immutableUsage--;
    }

    /**
     * \brief Assignment from \c other immutable vector.
     *
     * Assignment from \c other immutable vector.
     */

    Vector& operator=(const Vector& other)
    {
      assert(_data->_immutableUsage > 0);
      if (other._data != _data)
      {
        _data->_immutableUsage--;
        _data->deleteIfUnused();

        _data = other._data;
        _data->_immutableUsage++;
      }
      return *this;
    }

    /**
     * \brief Assignment from mutable vector \c other.
     *
     * Assignment from mutable vector \c other.
     */

    Vector& operator=(const MutableVector& other)
    {
      assert(_data->_immutableUsage > 0);
      _data->_immutableUsage--;
      _data->deleteIfUnused();

      _data = other._data;
      _data->_immutableUsage++;
      return *this;
    }

    void swap(Vector& other)
    {
      std::swap(_data, other._data);
    }
  };

  MutableVector zeroVector();
  MutableVector unitVector(std::size_t index);

  Rational operator*(const ReferenceCountedVector& a, const ReferenceCountedVector& b);
  Vector operator+(const Vector& a, const Vector& b);
  Vector operator-(const Vector& a, const Vector& b);
  MutableVector operator-(const ReferenceCountedVector& source);

  struct InnerDescription
  {
    std::vector<Vector> points;
    std::vector<Vector> rays;
  };

} /* namespace ipo */

#endif /* IPO_VECTORS_H_ */
