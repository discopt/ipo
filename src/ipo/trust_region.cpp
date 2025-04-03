#include "ipo/trust_region.hpp"

namespace ipo
{
  template <typename Number>
  ManhattanTrustRegion<Number>::ManhattanTrustRegion()
    : _center(nullptr), _size(0)
  {

  }

  template <typename Number>
  ManhattanTrustRegion<Number>::ManhattanTrustRegion(std::shared_ptr<sparse_vector<Number>> center, const Number& size)
    : _center(center), _size(size)
  {

  }

  template <typename Number>
  ManhattanTrustRegion<Number>::ManhattanTrustRegion(const ManhattanTrustRegion<Number>& other)
    : _center(other._center), _size(other._size)
  {

  }

  template <typename Number>
  ManhattanTrustRegion<Number>& ManhattanTrustRegion<Number>::operator=(const ManhattanTrustRegion<Number>& other)
  {
    this->_center = other._center;
    this->_size = other._size;
    return *this;
  }

  template <typename Number>
  ManhattanTrustRegion<Number>::~ManhattanTrustRegion()
  {

  }

  template class ManhattanTrustRegion<double>;

#if defined(IPO_WITH_RATIONAL)

  template class ManhattanTrustRegion<rational>;

#endif /* IPO_WITH_RATIONAL */

} /* namespace ipo */
