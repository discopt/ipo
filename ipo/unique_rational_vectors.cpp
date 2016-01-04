#include "unique_rational_vectors.h"

#include <limits>
#include <iostream>

using namespace soplex;

namespace ipo {

  std::size_t pointBitSize(const SVectorRational& vector)
  {
    std::size_t size = 0;
    for (int p = vector.size() - 1; p >= 0; --p)
      size += vector.value(p).sizeInBase(2);
    return size;
  }

  void pointDifference(Point& result, const Point& a, const Point& b)
  {
    result.clear();
    std::size_t aPos = 0;
    std::size_t bPos = 0;
    std::size_t aIndex, bIndex;
    while (aPos < a.size() && bPos < b.size())
    {
      aIndex = a.index(aPos);
      bIndex = b.index(bPos);
      if (aIndex < bIndex)
      {
        result.add(aIndex, a.value(aPos));
        aPos++;
      }
      else if (aIndex > bIndex)
      {
        result.add(bIndex, -b.value(bPos));
        bPos++;
      }
      else
      {
        if (a.value(aPos) != b.value(bPos))
          result.add(aIndex, a.value(aPos) - b.value(bPos));
        aPos++;
        bPos++;
      }
    }
    while (aPos < a.size())
    {
      assert(bPos == b.size());
      aIndex = a.index(aPos);
      result.add(aIndex, a.value(aPos));
      aPos++;
    }
    while (bPos < b.size())
    {
      assert(aPos == a.size());
      bIndex = b.index(bPos);
      result.add(bIndex, -b.value(bPos));
      bPos++;
    }
  }

  UniqueRationalVectorsBase::UniqueRationalVectorsBase(std::size_t numVariables) :
      _numVariables(numVariables)
  {

  }

  UniqueRationalVectorsBase::~UniqueRationalVectorsBase()
  {

  }

  UniqueRationalVectorsBase::Hash UniqueRationalVectorsBase::computeHash(const Point* vector) const
  {
    Hash hash = 0;
    for (std::size_t i = 0; i < vector->size(); ++i)
    {
      std::size_t var = vector->index(i);
      const soplex::Rational& val = vector->value(i);
      const mpq_t& q = val.getMpqRef();
      hash += (1773 * var + 1) * (mpz_get_si(mpq_numref(q)) * mpz_get_si(mpq_denref(q)));
    }
    return hash;
  }

  UniqueRationalVectors::UniqueRationalVectors(std::size_t numVariables) :
      UniqueRationalVectorsBase(numVariables)
  {

  }

  UniqueRationalVectors::~UniqueRationalVectors()
  {
    for (std::size_t i = 0; i < _vectors.size(); ++i)
    {
      delete _vectors[i];
      delete _approximations[i];
    }
  }

  std::size_t UniqueRationalVectors::insertCopy(const Point* vector)
  {
    Hash h;
    std::size_t index;
    HashMap::iterator hashMapIter;
    if (find(vector, h, index, hashMapIter))
      return index;
    else if (hashMapIter == _hashMap.end())
      hashMapIter = _hashMap.insert(std::make_pair(h, IndexList())).first;
    _vectors.push_back(new Point(*vector));
    DSVectorReal* approx = new DSVectorReal;
    for (int p = vector->size() - 1; p >= 0; --p)
      approx->add(vector->index(p), double(vector->value(p)));
    _approximations.push_back(approx);
    assert(index == _vectors.size() - 1);
    assert(index == _approximations.size() - 1);
    hashMapIter->second.insert(hashMapIter->second.end(), index);
    return index;
  }

  std::size_t UniqueRationalVectors::insertFree(Point* vector)
  {
    Point copy;
    copy = *vector;

    Hash h;
    std::size_t index;
    HashMap::iterator hashMapIter;
    if (find(vector, h, index, hashMapIter))
    {
      delete vector;
      return index;
    }
    else if (hashMapIter == _hashMap.end())
    {
      hashMapIter = _hashMap.insert(std::make_pair(h, IndexList())).first;
    }
    _vectors.push_back(vector);
    DSVectorReal* approx = new DSVectorReal;
    for (int p = vector->size() - 1; p >= 0; --p)
      approx->add(vector->index(p), double(vector->value(p)));
    _approximations.push_back(approx);
    assert(index == _vectors.size() - 1);
    assert(index == _approximations.size() - 1);
    hashMapIter->second.insert(hashMapIter->second.end(), index);

    return index;
  }

  std::size_t UniqueRationalVectors::first()
  {
    return 0;
  }

  std::size_t UniqueRationalVectors::next(std::size_t index)
  {
    return index + 1;
  }

  std::size_t UniqueRationalVectors::size()
  {
    return _vectors.size();
  }

  const DSVectorRational* UniqueRationalVectors::vector(std::size_t index) const
  {
    return _vectors[index];
  }

  const DSVectorReal* UniqueRationalVectors::approximation(std::size_t index) const
  {
    return _approximations[index];
  }

  std::size_t UniqueRationalVectors::nonzeros(std::size_t index) const
  {
    return _approximations[index]->size();
  }

  std::size_t UniqueRationalVectors::bitSize(std::size_t index) const
  {
    return _bitSizes[index];
  }

  soplex::DSVectorRational* UniqueRationalVectors::get(std::size_t index)
  {
    return _vectors[index];
  }

  void UniqueRationalVectors::extractAll()
  {
    _vectors.clear();
    for (std::size_t i = 0; i < _approximations.size(); ++i)
      delete _approximations[i];
    _approximations.clear();
    _hashMap.clear();
  }

  bool UniqueRationalVectors::find(const Point* vector, Hash& h, std::size_t& index, HashMap::iterator& hashMapIter)
  {
    h = computeHash(vector);
    hashMapIter = _hashMap.find(h);
    index = size();
    if (hashMapIter == _hashMap.end())
    {
      return false;
    }
    for (IndexList::const_iterator listIter = hashMapIter->second.begin(); listIter != hashMapIter->second.end();
        ++listIter)
    {
      const Point* p = _vectors[*listIter];
      if (p->size() != vector->size())
        continue;
      bool same = true;
      for (std::size_t i = 0; i < vector->size(); ++i)
      {
        if (vector->index(i) != p->index(i) || vector->value(i) != p->value(i))
        {
          same = false;
          break;
        }
      }
      if (same)
      {
        index = *listIter;
        return true;
      }
    }
    return false;
  }

  FilteredUniqueRationalVectors::FilteredUniqueRationalVectors(UniqueRationalVectorsBase& base) :
      UniqueRationalVectorsBase(base.numVariables()), _base(base)
  {
    _enabled.resize(base.size(), true);
  }

  FilteredUniqueRationalVectors::~FilteredUniqueRationalVectors()
  {

  }

  std::size_t FilteredUniqueRationalVectors::insertCopy(const Point* vector)
  {
    std::size_t index = _base.insertCopy(vector);
    update();
    return index;
  }

  std::size_t FilteredUniqueRationalVectors::insertFree(Point* vector)
  {
    std::size_t index = _base.insertFree(vector);
    update();
    return index;
  }

  std::size_t FilteredUniqueRationalVectors::first()
  {
    update();
    std::size_t index;
    for (index = _base.first(); index < _base.size(); index = _base.next(index))
    {
      if (_enabled[index])
        return index;
    }
    return std::numeric_limits<std::size_t>::max();
  }

  std::size_t FilteredUniqueRationalVectors::next(std::size_t index)
  {
    update();
    ++index;
    for (; index < _base.size(); index = _base.next(index))
    {
      if (_enabled[index])
        return index;
    }
    return std::numeric_limits<std::size_t>::max();
  }

  std::size_t FilteredUniqueRationalVectors::size()
  {
    return _base.size();
  }

  const soplex::DSVectorRational* FilteredUniqueRationalVectors::vector(std::size_t index) const
  {
    return _base.vector(index);
  }

  const DSVectorReal* FilteredUniqueRationalVectors::approximation(std::size_t index) const
  {
    return _base.approximation(index);
  }

  std::size_t FilteredUniqueRationalVectors::nonzeros(std::size_t index) const
  {
    return _base.nonzeros(index);
  }

  std::size_t FilteredUniqueRationalVectors::bitSize(std::size_t index) const
  {
    return _base.bitSize(index);
  }

  void FilteredUniqueRationalVectors::set(std::size_t index, bool enabled)
  {
    update();
    _enabled[index] = enabled;
  }

  void FilteredUniqueRationalVectors::setAll(bool enabled)
  {
    _enabled.clear();
    _enabled.resize(_base.size(), enabled);
  }

  void FilteredUniqueRationalVectors::update()
  {
    if (_enabled.size() < _base.size())
      _enabled.resize(_base.size(), true);
  }

//  std::ostream& HashedVectors::printStatistics(std::ostream& stream) const
//  {
//    stream << "PointSet has " << size() << " points in " << _hashMap.size()
//        << " buckets.\n";
//    std::size_t worstCardinality = 0;
//    const IndexList* worstList = NULL;
//    Hash worstHash = 0;
//    for (HashMap::const_iterator hashMapIter = _hashMap.begin();
//        hashMapIter != _hashMap.end(); ++hashMapIter)
//    {
//      if (hashMapIter->second.size() > worstCardinality)
//      {
//        worstList = &(hashMapIter->second);
//        worstCardinality = hashMapIter->second.size();
//        worstHash = hashMapIter->first;
//      }
//    }
//    if (worstList)
//    {
//      std::cerr << "Largest bucket has " << worstCardinality
//          << " points and hash " << worstHash << "\n";
//      for (IndexList::const_iterator listIter = worstList->begin();
//          listIter != worstList->end(); ++listIter)
//      {
//        stream << " Point in largest bucket:" << *_vectors[*listIter] << "\n";
//      }
//    }
//    return stream;
//  }

} /* namespace ipo */
