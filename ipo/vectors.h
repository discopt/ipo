#ifndef IPO_VECTORS_H_
#define IPO_VECTORS_H_

#include <set>

#include "common.h"
#include "rational.h"
#include "spx_gmp.h"

namespace ipo {

  typedef soplex::DVectorRational DenseVector;
  typedef soplex::DVectorReal DenseVectorApproximation;

  class SparseVector
  {
  protected:
    struct Nonzero
    {
      std::size_t index;
      Rational value;
      double approximation;
    };

    struct Implementation
    {
      std::size_t usage;
      std::size_t size;
      std::size_t capacity;
      Nonzero* nonzeros;
    };

  public:
    SparseVector(std::size_t capacity);
    SparseVector(SparseVector& other);
    ~SparseVector();

    SparseVector& operator=(SparseVector& other);

    inline std::size_t size() const
    {
      return impl->size;
    }

    inline std::size_t index(std::size_t position) const
    {
      assert(position < impl->size);
      return impl->nonzeros[position].index;
    }

    inline const Rational& value(std::size_t position) const
    {
      assert(position < impl->size);
      return impl->nonzeros[position].value;
    }

    inline double approximation(std::size_t position) const
    {
      assert(position < impl->size);
      return impl->nonzeros[position].approximation;
    }

    bool operator<(const SparseVector& other) const;

    void add(std::size_t index, const Rational& value);

    bool isSorted() const;

    void swap(SparseVector& other);

  protected:
  public: // TODO: for debugging
    Implementation* impl;
  };

  inline void swap(SparseVector& u, SparseVector& v)
  {
    u.swap(v);
  }

  class UniqueSparseVectors
  {
  protected:
    struct Data
    {
      double hash;
      SparseVector vector;

      Data(double hash, SparseVector& vector);
      Data(const Data& other);
      ~Data();

      bool operator<(const Data& other) const;
    };

  public:
    struct Iterator
    {
      Iterator(std::set<Data>::iterator iter)
        : _iter(iter)
      {

      }

      inline Iterator operator++()
      {
        ++_iter;
        return Iterator(_iter);
      }

      inline SparseVector& operator*()
      {
        return const_cast<SparseVector&>(_iter->vector);
      }
      
    private:
      std::set<Data>::iterator _iter;
    };

    UniqueSparseVectors(std::size_t ambientDimension);
    ~UniqueSparseVectors();

    bool insert(SparseVector& vector);

  private:
    std::set<Data> _data;
    DenseVectorApproximation _hashVector;
  };
}

#endif /* IPO_VECTORS_H_ */
