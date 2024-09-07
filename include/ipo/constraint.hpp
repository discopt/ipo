#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#include <ipo/sparse_vector.hpp>

#include <memory>
#include <limits>
#include <map>
#include <ostream>

namespace ipo
{
  enum class ConstraintType
  {
    EQUATION,
    LESS_OR_EQUAL,
    GREATER_OR_EQUAL,
    RANGED
  };

  template <typename T>
  class Constraint
  {
  public:
    Constraint(const T& lhs, std::shared_ptr<sparse_vector<T>> vector, const T& rhs,
      ConstraintType type)
      : _lhs(lhs), _vector(vector), _rhs(rhs), _type(type)
    {
      assert(_type != ConstraintType::RANGED || lhs <= rhs);
      assert(_type != ConstraintType::EQUATION || lhs == rhs);
    }

    Constraint(T&& lhs, std::shared_ptr<sparse_vector<T>> vector, T&& rhs, ConstraintType type)
      : _lhs(std::move(lhs)), _vector(vector), _rhs(std::move(rhs)), _type(type)
    {
      assert(_type != ConstraintType::RANGED || lhs <= rhs);
      assert(_type != ConstraintType::EQUATION || lhs == rhs);
    }

    Constraint(const T& lhs, std::shared_ptr<sparse_vector<T>> vector, const T& rhs)
      : _lhs(lhs), _vector(vector), _rhs(rhs), _type(lhs == rhs ? ConstraintType::EQUATION : ConstraintType::RANGED)
    {
      assert(lhs <= rhs);
    }

    Constraint(T&& lhs, std::shared_ptr<sparse_vector<T>> vector, T&& rhs)
      : _lhs(std::move(lhs)), _vector(std::move(vector)), _rhs(std::move(rhs)),
      _type(lhs == rhs ? ConstraintType::EQUATION : ConstraintType::RANGED)
    {
      assert(lhs <= rhs);
    }

    Constraint(const T& lhs, std::shared_ptr<sparse_vector<T>> vector)
      : _lhs(lhs), _vector(vector), _rhs(0), _type(ConstraintType::GREATER_OR_EQUAL)
    {

    }

    Constraint(T&& lhs, std::shared_ptr<sparse_vector<T>> vector)
      : _lhs(std::move(lhs)), _vector(std::move(vector)), _rhs(0),
      _type(ConstraintType::GREATER_OR_EQUAL)
    {

    }

    Constraint(std::shared_ptr<sparse_vector<T>> vector, const T& rhs)
      : _lhs(0), _vector(vector), _rhs(rhs), _type(ConstraintType::LESS_OR_EQUAL)
    {

    }

    Constraint(std::shared_ptr<sparse_vector<T>> vector, T&& rhs)
      : _lhs(0), _vector(std::move(vector)), _rhs(std::move(rhs)),
      _type(ConstraintType::LESS_OR_EQUAL)
    {

    }

    ~Constraint()
    {

    }

    inline bool operator==(const Constraint<T>& other) const
    {
      if (_vector != other._vector || _type != other._type)
        return false;
      if (hasLhs() && _lhs != other._lhs)
        return false;
      if (hasRhs() && _rhs != other._rhs)
        return false;
      return true;
    }

    inline bool operator!=(const Constraint<T>& other) const
    {
      return !(*this == other);
    }

    const T& lhs() const
    {
      return _lhs;
    }

    const sparse_vector<T>& vector() const
    {
      return *_vector;
    }

    std::shared_ptr<sparse_vector<T>> sharedVector()
    {
      return _vector;
    }

    ConstraintType type() const
    {
      return _type;
    }

    const T& rhs() const
    {
      return _rhs;
    }

    bool hasLhs() const
    {
      return _type != ConstraintType::LESS_OR_EQUAL;
    }

    bool hasRhs() const
    {
      return _type != ConstraintType::GREATER_OR_EQUAL;
    }

    bool isAlwaysSatisfied() const
    {
      return _vector->empty() && _lhs <= 0 && _rhs >= 0;
    }

    bool isNeverSatisfied() const
    {
      return (_vector->empty() && (_lhs > 0 || _rhs < 0)) || (_lhs > _rhs);
    }

    friend void scaleIntegral(Constraint<double>&);

#if defined(IPO_RATIONAL)

    friend void scaleIntegral(Constraint<rational>&);

#endif /* IPO_RATIONAL */

  protected:
    T _lhs;
    std::shared_ptr<sparse_vector<T>> _vector;
    T _rhs;
    ConstraintType _type;
  };

  template <typename T>
  Constraint<T> alwaysSatisfiedConstraint()
  {
    auto zero = std::make_shared<sparse_vector<T>>();
    return Constraint<T>(zero, T(1));
  }

  template <typename T>
  Constraint<T> neverSatisfiedConstraint()
  {
    auto zero = std::make_shared<sparse_vector<T>>();
    return Constraint<T>(zero, T(-1));
  }

  std::ostream& operator<<(std::ostream& stream, const Constraint<double>& constraint);

#if defined(IPO_RATIONAL)

  std::ostream& operator<<(std::ostream& stream, const Constraint<rational>& constraint);

#endif /* IPO_RATIONAL */

  void scaleIntegral(Constraint<double>& constraint);

#if defined(IPO_RATIONAL)

  void scaleIntegral(Constraint<rational>& constraint);

#endif /* IPO_RATIONAL */


  template <typename To, typename From>
  inline Constraint<To> convertConstraint(const Constraint<From>& from)
  {
    To lhs = convertNumber<To, From>(from.lhs());
    To rhs = convertNumber<To, From>(from.rhs());
    return Constraint<To>(lhs,
      std::make_shared<sparse_vector<To>>(convertSparseVector<To, From>(from.vector())), rhs, from.type());
  }

  template <typename Number>
  class ConstraintSet
  {
  protected:
    struct Payload
    {
      Constraint<Number> constraint; /// The actual constraint.
      double hash; /// The scalar product of the normal vector with the random vector used for hashing.
      double inverseNorm; /// Inverse of Euclidean norm of normal vector.

      Payload(const Constraint<Number>& cons, double hsh)
        : constraint(cons), hash(hsh)
      {
        inverseNorm = euclideanNorm(constraint.vector());
        if (inverseNorm > 0.0)
          inverseNorm = 1.0 / inverseNorm;
      }
    };

    struct const_iterator
    {
      typename std::vector<Payload>::const_iterator iterator;

      const_iterator(typename std::vector<Payload>::const_iterator iter)
        : iterator(iter)
      {

      }

      const Constraint<Number>& operator*() const
      {
        return iterator->constraint;
      }

      bool operator!=(const const_iterator& other) const
      {
        return iterator != other.iterator;
      }
      
      const_iterator operator++()
      {
        ++iterator;
        return *this;
      }
    };

    bool exists(const Constraint<Number>& constraint, double hashProduct);

    void force_push_back(const Constraint<Number>& constraint, double hashProduct)
    {
      Payload payload(constraint, hashProduct);
      _data.emplace_back(payload);
      _hashToIndex.insert(std::make_pair(hashProduct, _data.size() - 1));
    }

    size_t _n;
    std::vector<Payload> _data;
    double* _hashVector;
    std::multimap<double, std::size_t> _hashToIndex;

  public:
    ConstraintSet(std::size_t ambientDimension)
      : _n(ambientDimension)
    {
      _hashVector = generateRandomVectorSphere(ambientDimension);
    }

    ConstraintSet(std::size_t ambientDimension, const std::vector<Constraint<Number>>& constraints)
      : _n(ambientDimension)
    {
      _hashVector = generateRandomVectorSphere(ambientDimension);

      for (const auto& constraint : constraints)
        push_back(constraint);
    }

    ~ConstraintSet()
    {
      delete[] _hashVector;
    }

    const_iterator begin() const
    {
      return const_iterator(_data.begin());
    }

    const_iterator end() const
    {
      return const_iterator(_data.end());
    }

    bool push_back(const Constraint<Number>& constraint)
    {
      double product = 0.0;
      for (const auto& iter : constraint.vector())
        product += convertNumber<double>(iter.second) * _hashVector[iter.first];
      product = fabs(product);

      if (exists(constraint, product))
        return false;

      force_push_back(constraint, product);
      return true;
    }

    void force_push_back(const Constraint<Number>& constraint)
    {
      double product = 0.0;
      for (const auto& iter : constraint.vector())
        product += iter.second * _hashVector[iter.first];
      product = fabs(product);

      force_push_back(constraint, product);
    }

    bool exists(const Constraint<Number>& constraint)
    {
      double product = 0.0;
      for (const auto& iter : constraint.vector())
        product += convertNumber<double>(iter.second) * _hashVector[iter.first];
      product = fabs(product);

      return exists(constraint, product);
    }
  };
}
