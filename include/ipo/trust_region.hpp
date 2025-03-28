#pragma once

#include <ipo/arithmetic.hpp>
#include <ipo/sparse_vector.hpp>

namespace ipo
{
  /**
   * \brief Information about a Manhattan-norm trust region.
   *
   * A trust region is parameterized by its center point \f$ c \in \mathbb{R}^n \f$ and its size \f$ s > 0 \f$ and is
   * defined as the set \f$ \{ x \in \mathbb{R}^n \mid ||x-c||_1 \leq s \} \f$ of points of Manhattan distance at
   * most \f$ s \f$ from \f$ c \f$.
   */

  template <typename NumberType>
  class ManhattanTrustRegion
  {
  public:
    /**
     * \brief Creates the whole space as trust region.
     */

    ManhattanTrustRegion();

    /**
     * \brief Creates a trust region with given \p center and \p size.
     */

    ManhattanTrustRegion(std::shared_ptr<sparse_vector<NumberType>> center, const NumberType& size);

    /**
     * \brief Copy constructor for a trust region.
     */

    ManhattanTrustRegion(const ManhattanTrustRegion<NumberType>& other);

    /**
     * \brief Assignment.
     */

    ManhattanTrustRegion<NumberType>& operator=(const ManhattanTrustRegion<NumberType>& other);

    /**
     * \brief Destructor.
     */

    virtual ~ManhattanTrustRegion();

    /**
     * \brief Returns whether the trust region is bounded, i.e., not the full space.
     */

    inline bool isBounded() const
    {
      return _center != nullptr;
    }

    /**
     * \brief Returns the center point \f$ c \f$ of the trust region.
     */

    inline const std::shared_ptr<sparse_vector<NumberType>> center() const
    {
      return _center;
    }

    /**
     * \brief Returns the size \f$ s \f$ of the trust region.
     */

    inline NumberType size() const
    {
      return _size;
    }

    /**
     * \brief Update the size \f$ s \f$ of the trust region.
     */

    inline void updateSize(const NumberType& newSize)
    {
      _size = newSize;
    }

  private:
    /// Center point of the trust region.
    std::shared_ptr<sparse_vector<NumberType>> _center;
    /// Maximum Manhattan distance to center.
    NumberType _size;
  };

} /* namespace ipo */

