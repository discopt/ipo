#ifndef IPO_UNIQUE_RATIONAL_VECTORS_H_
#define IPO_UNIQUE_RATIONAL_VECTORS_H_

#include <list>
#include <map>
#include <vector>

#include "ipo.h"
#include "spx_gmp.h"

namespace ipo {

  /// A point.

  typedef soplex::DSVectorRational Point;

  /// An (unbounded) direction.

  typedef soplex::DSVectorRational Direction;

  /**
   * \brief Computes the bitsize of a vector.
   *
   * Computes the bitsize of a vector.
   */

  std::size_t bitSize(const soplex::SVectorRational& vector);

  /// A index subset.

  typedef std::vector<std::size_t> VectorSubset;

  // TODO: Remove inheritance and filtered version: Deprecated due to CacheOracle.
  
  /**
   * \brief Defines a container for unique sparse rational vectors.
   *
   * A container for sparse rational vectors that ensures uniqueness.
   * It also stores floating-point versions.
   */

  class UniqueRationalVectorsBase
  {
  protected:
    typedef std::size_t Hash;

  public:
    /**
     * \brief Constructor.
     */

    UniqueRationalVectorsBase(std::size_t numVariables);

    /**
     * \brief Destructor.
     */

    virtual ~UniqueRationalVectorsBase();

    /**
     * \brief Inserts a copy of a vector.
     *
     * Inserts a copy of a vector.
     */

    virtual std::size_t insertCopy(const Point* vector) = 0;

    /**
     * \brief Inserts a vector and frees it, if duplicate.
     *
     * Inserts a vector and frees it, if duplicate.
     */

    virtual std::size_t insertFree(Point const* vector) = 0;

    /**
     * \brief Returns the index of first vector.
     *
     * Returns the index of first vector.
     */

    virtual std::size_t first() = 0;

    /**
     * \brief Returns the index of the vector coming after \c index.
     *
     * Returns the index of the vector coming after \c index.
     * If \c index is the last vector, it returns
     * some index greater than or equal to \ref size().
     */

    virtual std::size_t next(std::size_t index) = 0;

    /**
     * \brief Returns an upper bound on the indices for the vectors.
     *
     * Returns an upper bound on the indices for the vectors.
     */

    virtual std::size_t size() const = 0;

    /**
     * \brief Returns the vector of a given index.
     *
     * Returns a const reference to the vector at \c index.
     */

    inline const soplex::DSVectorRational* operator[](std::size_t index) const
    {
      return vector(index);
    }

    /**
     * \brief Returns the vector of a given index.
     *
     * Returns a const reference to the vector at \c index.
     */

    virtual const soplex::DSVectorRational* vector(std::size_t index) const = 0;

    /**
     * \brief Returns the floating-point vector of a given index.
     *
     * Returns a const reference to the floating-point version of the vector at \c index.
     */

    virtual const soplex::DSVectorReal* approximation(std::size_t index) const = 0;

    /**
     * \brief Returns the number of nonzeros of vector at \c index.
     *
     * Returns the number of nonzeros of vector at \c index.
     */

    virtual std::size_t nonzeros(std::size_t index) const = 0;

    /**
     * \brief Returns the encoding length of the vector at \c index.
     *
     * Returns the encoding length of the vector at \c index.
     */

    virtual std::size_t bitSize(std::size_t index) const = 0;

    /**
     * \brief Returns the number of vectors by counting them.
     *
     * Returns the number of vectors by counting them.
     */

    inline std::size_t count()
    {
      std::size_t cnt = 0;
      for (std::size_t index = first(); index < size(); index = next(index))
        ++cnt;
      return cnt;
    }

    /**
     * \brief Returns the ambient dimension of the vectors.
     */

    inline std::size_t numVariables() const
    {
      return _numVariables;
    }

  protected:
    /**
     * \brief Computes the hash of a vector.
     *
     * Computes the hash of a \c vector.
     */

    Hash computeHash(const soplex::DSVectorRational* vector) const;

    std::size_t _numVariables; // Ambient dimension.
  };

  /**
   * \brief A container for unique sparse rational vectors.
   *
   * A container for unique sparse rational vectors.
   * It is the default implementation of \c UniqueRationalVectorsBase.
   */

  class UniqueRationalVectors: public UniqueRationalVectorsBase
  {
  protected:
    typedef std::list<std::size_t> IndexList;
    typedef std::map<Hash, IndexList> HashMap;

  public:
    /**
     * \brief Constructor.
     */

    UniqueRationalVectors(std::size_t numVariables);

    /**
     * \brief Destructor.
     *
     * Destructor.
     * It frees all stored vectors.
     */

    virtual ~UniqueRationalVectors();

    /**
     * \brief Inserts a copy of a vector.
     *
     * Inserts a copy of a vector.
     */

    virtual std::size_t insertCopy(const Point* vector);

    /**
     * \brief Inserts a vector and frees it, if duplicate.
     *
     * Inserts a vector and frees it, if duplicate.
     */

    virtual std::size_t insertFree(Point const* vector);

    /**
     * \brief Returns the index of first vector.
     *
     * Returns the index of first vector.
     */

    virtual std::size_t first();

    /**
     * \brief Returns the index of the vector coming after \c index.
     *
     * Returns the index of the vector coming after \c index.
     * If \c index is the last vector, it returns
     * some index greater than or equal to \ref size().
     */

    virtual std::size_t next(std::size_t index);

    /**
     * \brief Returns an upper bound on the indices for the vectors.
     *
     * Returns an upper bound on the indices for the vectors.
     */

    virtual std::size_t size() const;

    /**
     * \brief Returns the vector of a given index.
     *
     * Returns a const reference to the vector at \c index.
     */

    virtual const soplex::DSVectorRational* vector(std::size_t index) const;

    /**
     * \brief Returns the number of nonzeros of vector at \c index.
     *
     * Returns the number of nonzeros of vector at \c index.
     */

    virtual const soplex::DSVectorReal* approximation(std::size_t index) const;

    /**
     * \brief Returns the number of nonzeros of vector at \c index.
     *
     * Returns the number of nonzeros of vector at \c index.
     */

    virtual std::size_t nonzeros(std::size_t index) const;

    /**
     * \brief Returns the encoding length of the vector at \c index.
     *
     * Returns the encoding length of the vector at \c index.
     */

    virtual std::size_t bitSize(std::size_t index) const;

    /**
     * \brief Returns a writable version of the vector at \c index.
     *
     * Returns a writable version of the vector at \c index.
     */

    soplex::DSVectorRational const* get(std::size_t index);

    /**
     * \brief Clears the container without freeing the vectors.
     *
     * Clears the container without freeing the vectors.
     * The approximate versions are freed, though.
     */

    void extractAll();

  protected:
    /**
     * \brief Finds a given vector using a hash map.
     *
     * Finds a given vector using a hash map.
     */

    bool find(const Point* vector, Hash& hash, std::size_t& index, HashMap::iterator& hashMapIter);

    std::vector<Point const*> _vectors; // The rational sparse vectors.
    std::vector<soplex::DSVectorReal*> _approximations; // The floating-point versions of the vectors.
    std::vector<std::size_t> _bitSizes; // The bitsizes of the vectors.
    HashMap _hashMap; // The hash map to speed-up searching for duplicates.
  };

  /**
   * \brief A filtered view onto another \ref UniqueRationalVectorsBase container.
   *
   * References another \ref UniqueRationalVectorsBase container,
   * allowing the user to mark every vector as visible or invisible.
   * Iteration using \ref first() and \ref next()
   * will only return indices of visible vectors.
   */

  class FilteredUniqueRationalVectors: public UniqueRationalVectorsBase
  {
  protected:
    typedef std::list<std::size_t> IndexList;
    typedef std::map<Hash, IndexList> HashMap;

  public:
    /**
     * \brief Constructor referencing another container.
     *
     * Constructor referencing another container.
     */

    FilteredUniqueRationalVectors(UniqueRationalVectorsBase& base);

    /**
     * \brief Destructor.
     *
     * Destructor.
     * It does not free any vectors of the referenced container.
     */

    virtual ~FilteredUniqueRationalVectors();

    /**
     * \brief Inserts a copy of a vector into the referenced container.
     *
     * Inserts a copy of a vector into the referenced container
     * and marks it as visible.
     */

    virtual std::size_t insertCopy(const Point* vector);

    /**
     * \brief Inserts a vector and frees it, if duplicate.
     *
     * Inserts a vector into the referenced container and frees it, if duplicate.
     * If it is added, it is also marked as visible.
     */

    virtual std::size_t insertFree(Point const* vector);

    /**
     * \brief Returns the index of first visible vector.
     *
     * Returns the index of first visible vector.
     */

    virtual std::size_t first();

    /**
     * \brief Returns the index of the vector coming after \c index.
     *
     * Returns the index of the vector coming after \c index.
     * If \c index is the last vector, it returns
     * some index greater than or equal to \ref size().
     */

    virtual std::size_t next(std::size_t index);

    /**
     * \brief Returns an upper bound on the indices for the vectors.
     *
     * Returns an upper bound on the indices for the vectors.
     */

    virtual std::size_t size() const;

    /**
     * \brief Returns the vector of a given index.
     *
     * Returns a const reference to the vector at \c index.
     */

    virtual const soplex::DSVectorRational* vector(std::size_t index) const;

    /**
     * \brief Returns the number of nonzeros of vector at \c index.
     *
     * Returns the number of nonzeros of vector at \c index.
     */

    virtual const soplex::DSVectorReal* approximation(std::size_t index) const;

    /**
     * \brief Returns the number of nonzeros of vector at \c index.
     *
     * Returns the number of nonzeros of vector at \c index.
     */

    virtual std::size_t nonzeros(std::size_t index) const;

    /**
     * \brief Returns the encoding length of the vector at \c index.
     *
     * Returns the encoding length of the vector at \c index.
     */

    virtual std::size_t bitSize(std::size_t index) const;

    /**
     * \brief Sets the visibility of the vector at \c index.
     *
     * Sets the visibility of the vector at \c index.
     */

    void set(std::size_t index, bool enabled);

    /**
     * \brief Sets the visibility of all vectors.
     *
     * Sets the visibility of all vectors.
     */

    void setAll(bool enabled);

    /**
     * \brief Marks the vector at \c index as visible.
     *
     * Marks the vector at \c index as visible.
     */

    inline void enable(std::size_t index)
    {
      set(index, true);
    }

    /**
     * \brief Marks the vector at \c index as invisible.
     *
     * Marks the vector at \c index as invisible.
     */

    inline void disable(std::size_t index)
    {
      set(index, false);
    }

    /**
     * \brief Returns whether the vector at \c index is visible.
     *
     * Returns whether the vector at \c index is visible.
     */

    inline bool get(std::size_t index) const
    {
      return _enabled[index];
    }

  protected:
    /**
     * \brief Updates data structures.
     *
     * Updates data structures in case
     * the referenced container changed.
     */

    void update();

    UniqueRationalVectorsBase& _base; // Referenced container.
    std::vector<bool> _enabled; // Visibility of the vectors.
  };

} /* namespace ipo */

#endif /* IPO_UNIQUE_RATIONAL_VECTORS_H_ */
