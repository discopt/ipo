#pragma once

#include <ostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <cmath>

#include <ipo/arithmetic.hpp>

template <class T>
class sparse_vector
{
public:
  typedef std::size_t key_type;
  typedef T mapped_type;
  typedef std::pair<key_type, mapped_type> value_type;
  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

  sparse_vector()
  {

  }

  sparse_vector(const sparse_vector<T>& other)
    : _data(other._data)
  {

  }

  sparse_vector(sparse_vector<T>&& other)
    : _data(std::move(other._data))
  {

  }

  ~sparse_vector()
  {

  }

  sparse_vector<T>& operator=(const sparse_vector<T>& other)
  {
    _data = other._data;
    return *this;
  }

  sparse_vector<T>& operator=(sparse_vector<T>&& other)
  {
    _data = std::move(other._data);
    return *this;
  }

  sparse_vector(std::vector<value_type>&& other, bool sort)
    : _data(std::move(other))
  {
    struct Less
    {
      bool operator()(const value_type& a, const value_type& b) const
      {
        return a.first < b.first;
      }
    };

    if (sort)
      std::sort(_data.begin(), _data.end(), Less());
  }

  std::size_t size() const
  {
    return _data.size();
  }

  bool empty() const
  {
    return _data.empty();
  }

  void clear()
  {
    _data.clear();
  }

  void push_back(key_type key, const mapped_type& value)
  {
    assert(_data.empty() || key > _data.back().first);
    _data.push_back(std::make_pair(key, value));
  }

  void push_back(key_type key, mapped_type&& value)
  {
    assert(_data.empty() || key > _data.back().first);
    _data.push_back(std::make_pair(key, value));
  }

  iterator begin()
  {
    return _data.begin();
  }

  iterator end()
  {
    return _data.end();
  }

  const_iterator begin() const
  {
    return _data.begin();
  }

  const_iterator end() const
  {
    return _data.end();
  }

  iterator find(const key_type& key)
  {
#if !defined(NDEBUG)
    for (std::size_t i = 1; i < _data.size(); ++i)
      assert(_data[i-1].first < _data[i].first);
#endif /* !NDEBUG */

    std::size_t left = 0;
    std::size_t right = _data.size();
    while (left < right)
    {
      std::size_t middle = (left + right) / 2;
      std::size_t k = _data[middle].first;
      if (key < k)
        right = middle;
      else if (key > k)
        left = middle + 1;
      else
        return _data.begin() + middle;
    }
    return _data.end();
  }

  const_iterator find(const key_type& key) const
  {
#if !defined(NDEBUG)
    for (std::size_t i = 1; i < _data.size(); ++i)
      assert(_data[i-1].first < _data[i].first);
#endif /* !NDEBUG */
    
    std::size_t left = 0;
    std::size_t right = _data.size();
    while (left < right)
    {
      std::size_t middle = (left + right) / 2;
      std::size_t k = _data[middle].first;
      if (key < k)
        right = middle;
      else if (key > k)
        left = middle + 1;
      else
        return _data.begin() + middle;
    }
    return _data.end();
  }

  const mapped_type& find(key_type index, const mapped_type& backup) const
  {
    auto iter = find(index);
    return (iter != end()) ? iter->second : backup;
  }

  friend void swap(sparse_vector<T>& a, sparse_vector<T>& b)
  {
    std::swap(a._data, b._data);
  }

  std::size_t getInequalIndex(const sparse_vector<T>& other) const
  {
    // We first compare the support.

    sparse_vector<T>::const_iterator iter1 = begin();
    sparse_vector<T>::const_iterator iter2 = other.begin();
    while (true)
    {
      if (iter1 == end())
      {
        if (iter2 == other.end())
          return std::numeric_limits<std::size_t>::max();
        else
          return iter2->first;
      }
      else if (iter2 == other.end())
        return iter1->first;
      else
      {
        ++iter1;
        ++iter2;
      }
    }

    // We now know that they have the same support, so we have to compare values.

    iter1 = begin();
    iter2 = other.begin();
    while (iter1 != end())
    {
      assert(iter2 != end());
      assert(iter1->first == iter2->first);
      if (iter1->second != iter2->second)
        return iter1->first;
      ++iter1;
      ++iter2;
    }
    return std::numeric_limits<std::size_t>::max();
  }

protected:
  std::vector<value_type> _data;
};

template <typename T,typename U>
T operator*(const sparse_vector<T>& a, const sparse_vector<U>& b)
{
  T result(0);
  typename sparse_vector<T>::const_iterator ia = a.begin();
  if (ia == a.end())
    return result;
  typename sparse_vector<U>::const_iterator ib = b.begin();
  if (ib == b.end())
    return result;
  while (true)
  {
    if (ia->first < ib->first)
    {
      ++ia;
      if (ia == a.end())
        return result;
    }
    else if (ia->first > ib->first)
    {
      ++ib;
      if (ib == b.end())
        return result;
    }
    else
    {
      assert(ia->first == ib->first);
      result += ia->second * convertNumber<T,U>(ib->second);
      ++ia;
      if (ia == a.end())
        return result;
      ++ib;
      if (ib == b.end())
        return result;
    }
  }
  return result;
}

template <typename T, typename U>
T operator*(const sparse_vector<T>& a, const U* b)
{
  T result(0);
  for (auto& iter : a)
    result += iter.second * convertNumber<T>(b[iter.first]);
  return result;
}

template <typename T, typename U>
T operator*(const T* a, const sparse_vector<U>& b)
{
  T result(0);
  for (auto& iter : b)
    result += a[iter.first] * convertNumber<T>(iter.second);
  return result;
}

template <typename T>
std::ostream& operator<<(std::ostream& stream, const sparse_vector<T>& vector)
{
  stream << '[';
  bool first = true;
  for (auto pair : vector)
  {
    if (first)
      first = false;
    else
      stream << ',';
    stream << pair.first << ":" << pair.second;
  }
  return stream << ']';
}

template <typename T, typename U>
double squaredEuclideanDistance(const sparse_vector<T>& a, const sparse_vector<U>& b)
{
  double result = 0.0;
  typename sparse_vector<T>::const_iterator ia = a.begin();
  typename sparse_vector<U>::const_iterator ib = b.begin();
  if (ia != a.end() && ib != b.end())
  {
    while (true)
    {
      if (ia->first < ib->first)
      {
        double x = convertNumber<double>(ia->second);
        result += x*x;
        ++ia;
        if (ia == a.end())
          break;
      }
      else if (ia->first > ib->first)
      {
        double x = convertNumber<double>(ib->second);
        result += x*x;
        ++ib;
        if (ib == b.end())
          break;
      }
      else
      {
        double x = convertNumber<double, T>(ia->second) - convertNumber<double, U>(ib->second);
        result += x*x;
        ++ia;
        ++ib;
        if (ia == a.end() || ib == b.end())
          break;
      }
    }
  }
  for (; ia != a.end(); ++ia)
  {
    double x = convertNumber<double>(ia->second);
    result += x*x;
  }
  for (; ib != b.end(); ++ib)
  {
    double x = convertNumber<double>(ib->second);
    result += x*x;
  }
  return result;
}

template <typename T, typename U>
double euclideanDistance(const sparse_vector<T>& a, const sparse_vector<U>& b)
{
  return sqrt(squaredEuclideanDistance(a, b));
}

template <typename T>
double squaredEuclideanNorm(const sparse_vector<T>& vector)
{
  double result = 0.0;
  for (const auto& iter : vector)
  {
    double x = convertNumber<double>(iter.second);
    result += x*x;
  }
  return result;
}

template <typename T>
double euclideanNorm(const sparse_vector<T>& vector)
{
  return sqrt(squaredEuclideanNorm(vector));
}

template<typename To, typename From>
inline sparse_vector<To> convertTo(const sparse_vector<From>& vector)
{
  std::vector<std::pair<std::size_t, To>> entries;
  entries.reserve(vector.size());
  for (const auto& iter : vector)
    entries.push_back(std::make_pair(iter.first, convertNumber<To>(iter.second)));
  return sparse_vector<To>(std::move(entries), false);
}
