#ifndef IPO_VECTOR_SPACE_GENERATORS_H_
#define IPO_VECTOR_SPACE_GENERATORS_H_

#include <vector>

#include "common.h"
#include "rational.h"

namespace ipo {

  class VectorSpaceGenerators
  {
  public:
    VectorSpaceGenerators();
    virtual ~VectorSpaceGenerators();

    void reset(std::size_t ambientDimension);
    std::size_t add(const soplex::SVectorRational& vector);
    std::size_t add(const soplex::SVectorRational& vector, bool dependent);
    std::size_t addLazy(const soplex::SVectorRational& vector);
    void flushLazy();
    bool isDependent(std::size_t index);
    bool isDependent(const soplex::SVectorRational& vector);
    bool isDependent(const soplex::VectorRational& vector);

  protected:
    void checkVector(std::size_t index);

    soplex::SoPlex* _spx;
    std::vector<std::size_t> _lazyVectors;
    bool _zeroRhs;
    soplex::DVectorRational _solution;
  };


} /* namespace ipo */

#endif /* IPO_VECTOR_SPACE_GENERATORS_H_ */
