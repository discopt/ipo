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

      bool operator<(const Nonzero& other) const;
    };

    struct Implementation
    {
      std::size_t usage;
      std::size_t size;
      std::size_t capacity;
      Nonzero* nonzeros;
    };

  public:
    SparseVector();
    SparseVector(std::size_t capacity);
    SparseVector(const SparseVector& other);
    ~SparseVector();

    SparseVector& operator=(const SparseVector& other);

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

    bool operator==(const SparseVector& other) const;

    bool operator<(const SparseVector& other) const;

    void add(std::size_t index, const Rational& value);

    void sort();

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

  SparseVector denseToSparseVector(const DenseVector& source);
  void assign(DenseVectorApproximation& target, const DenseVector& source);
  void assign(DenseVector& target, const DenseVector& source);
  void assign(DenseVector& target, const SparseVector& source);
  void assign(SparseVector& target, const DenseVector& source);
  void assign(SparseVector& target, const SparseVector& source);
  Rational scalarProduct(const DenseVector& a, const DenseVector& b);
  Rational scalarProduct(const DenseVector& a, const SparseVector& b);
  Rational scalarProduct(const SparseVector& a, const DenseVector& b);
  Rational scalarProduct(const SparseVector& a, const SparseVector& b);
  std::size_t differingIndex(const SparseVector& a, const SparseVector& b);
  void add(DenseVector& target, const SparseVector& source);
  void subtract(DenseVector& target, const SparseVector& source);

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
      
      inline bool operator==(const Iterator& other) const
      {
        return _iter == other._iter;
      }

      inline bool operator!=(const Iterator& other) const
      {
        return ! (*this == other);
      }
      
    private:
      std::set<Data>::iterator _iter;
    };

    UniqueSparseVectors(std::size_t ambientDimension);
    ~UniqueSparseVectors();

    inline std::size_t size() const
    {
      return _data.size();
    }

    void clear();
    
    bool insert(SparseVector& vector);

    inline Iterator begin()
    {
      return Iterator(_data.begin());
    }

    inline Iterator end()
    {
      return Iterator(_data.end());
    }

  private:
    std::set<Data> _data;
    DenseVectorApproximation _hashVector;
  };
}

#endif /* IPO_VECTORS_H_ */
