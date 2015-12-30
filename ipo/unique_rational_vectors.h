#ifndef IPO_UNIQUE_RATIONAL_VECTORS_H_
#define IPO_UNIQUE_RATIONAL_VECTORS_H_

#include <list>
#include <map>
#include <vector>

#include "spx_gmp.h"

namespace ipo {

  typedef soplex::DSVectorRational Point;
  typedef soplex::DSVectorRational Ray;

  std::size_t bitSize(const soplex::SVectorRational& point);

  // TODO: deprecate this!
  void pointDifference(Point& result, const Point& a, const Point& b);

  typedef std::vector<std::size_t> VectorSubset;

  class UniqueRationalVectorsBase
  {
  protected:
    typedef std::size_t Hash;

  public:
    UniqueRationalVectorsBase(std::size_t numVariables);
    virtual ~UniqueRationalVectorsBase();

    virtual std::size_t insertCopy(const Point* vector) = 0;
    virtual std::size_t insertFree(Point* vector) = 0;
    virtual std::size_t first() = 0;
    virtual std::size_t next(std::size_t index) = 0;
    virtual std::size_t size() = 0;

    inline const soplex::DSVectorRational* operator[](std::size_t index) const
    {
      return vector(index);
    }

    virtual const soplex::DSVectorRational* vector(std::size_t index) const = 0;
    virtual const soplex::DSVectorReal* approximation(std::size_t index) const = 0;
    virtual std::size_t nonzeros(std::size_t index) const = 0;
    virtual std::size_t bitSize(std::size_t index) const = 0;

    inline std::size_t count()
    {
      std::size_t cnt = 0;
      for (std::size_t index = first(); index < size(); index = next(index))
        ++cnt;
      return cnt;
    }

    inline std::size_t numVariables() const
    {
      return _numVariables;
    }

  protected:
    Hash computeHash(const soplex::DSVectorRational* vector) const;

    std::size_t _numVariables;
  };

  class UniqueRationalVectors: public UniqueRationalVectorsBase
  {
  protected:
    typedef std::list<std::size_t> IndexList;
    typedef std::map<Hash, IndexList> HashMap;

  public:
    UniqueRationalVectors(std::size_t numVariables);
    virtual ~UniqueRationalVectors();

    virtual std::size_t insertCopy(const Point* vector);
    virtual std::size_t insertFree(Point* vector);
    virtual std::size_t first();
    virtual std::size_t next(std::size_t index);
    virtual std::size_t size();
    virtual const soplex::DSVectorRational* vector(std::size_t index) const;
    virtual const soplex::DSVectorReal* approximation(std::size_t index) const;
    virtual std::size_t nonzeros(std::size_t index) const;
    virtual std::size_t bitSize(std::size_t index) const;

    soplex::DSVectorRational* get(std::size_t index);
    void extractAll();

  protected:
    bool find(const Point* vector, Hash& hash, std::size_t& index, HashMap::iterator& hashMapIter);

    std::vector<Point*> _vectors;
    std::vector<soplex::DSVectorReal*> _approximations;
    std::vector<std::size_t> _bitSizes;
    HashMap _hashMap;
  };

  class FilteredUniqueRationalVectors: public UniqueRationalVectorsBase
  {
  protected:
    typedef std::list<std::size_t> IndexList;
    typedef std::map<Hash, IndexList> HashMap;

  public:
    FilteredUniqueRationalVectors(UniqueRationalVectorsBase& base);
    virtual ~FilteredUniqueRationalVectors();

    virtual std::size_t insertCopy(const Point* vector);
    virtual std::size_t insertFree(Point* vector);
    virtual std::size_t first();
    virtual std::size_t next(std::size_t index);
    virtual std::size_t size();
    virtual const soplex::DSVectorRational* vector(std::size_t index) const;
    virtual const soplex::DSVectorReal* approximation(std::size_t index) const;
    virtual std::size_t nonzeros(std::size_t index) const;
    virtual std::size_t bitSize(std::size_t index) const;

    void set(std::size_t index, bool enabled);
    void setAll(bool enabled);

    inline void enable(std::size_t index)
    {
      set(index, true);
    }

    inline void disable(std::size_t index)
    {
      set(index, false);
    }

    inline bool get(std::size_t index) const
    {
      return _enabled[index];
    }

  protected:
    void update();

    UniqueRationalVectorsBase& _base;
    std::vector<bool> _enabled;
  };

} /* namespace polycomb */

#endif /* IPO_UNIQUE_RATIONAL_VECTORS_H_ */
