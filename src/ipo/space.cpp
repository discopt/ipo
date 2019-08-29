#include <ipo/space.hpp>

#include <ostream>
#include <sstream>

namespace ipo
{
  Space::Space()
  {

  }

  Space::Space(const std::vector<std::string>& variableNames)
    : _variableNames(variableNames)
  {

  }

  bool Space::operator==(const Space& other) const
  {
    return _variableNames == other._variableNames;
  }

  void Space::addVariable(const std::string& name)
  {
    _variableNames.push_back(name);
  }

  void Space::printVector(std::ostream& stream, std::size_t numNonzeros,
    std::size_t* nonzeroCoordinates, double* nonzeroValues) const
  {
    if (numNonzeros == 0)
      stream << "()";
    else
    {
      stream << "(" << _variableNames[nonzeroCoordinates[0]] << "=" << nonzeroValues[0];
      for (std::size_t i = 1; i < numNonzeros; ++i)
        stream << "," << _variableNames[nonzeroCoordinates[i]] << "=" << nonzeroValues[i];
      stream << ")";
    }
  }

  std::string Space::vectorToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
    double* nonzeroValues) const
  {
    std::stringstream stream;
    printVector(stream, numNonzeros, nonzeroCoordinates, nonzeroValues);
    return stream.str();
  }

  void Space::printLinearForm(std::ostream& stream, std::size_t numNonzeros,
    std::size_t* nonzeroCoordinates, double* nonzeroValues) const
  {
    if (numNonzeros == 0)
        stream << "0";
    else
    {
      double x0 = nonzeroValues[0];
      if (x0 == -1)
        stream << "-";
      else if (x0 != 1)
        stream << x0;
      stream << _variableNames[nonzeroCoordinates[0]];
      for (std::size_t i = 1; i < numNonzeros; ++i)
      {
        double xi = nonzeroValues[i];
        if (xi == 1)
          stream << " + ";
        else if (xi == -1)
          stream << " - ";
        else if (xi > 0)
          stream << " + " << xi;
        else
          stream << " - " << (-xi);
        stream << _variableNames[nonzeroCoordinates[i]];
      }
    }
  }

  std::string Space::linearFormToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
    double* nonzeroValues) const
  {
    std::stringstream stream;
    printLinearForm(stream, numNonzeros, nonzeroCoordinates, nonzeroValues);
    return stream.str();
  }

#if defined(IPO_WITH_GMP)

  void Space::printVector(std::ostream& stream, std::size_t numNonzeros,
    std::size_t* nonzeroCoordinates, mpq_class* nonzeroValues) const
  {
    if (numNonzeros == 0)
      stream << "()";
    else
    {
      stream << "(" << _variableNames[nonzeroCoordinates[0]] << "=" << nonzeroValues[0];
      for (std::size_t i = 1; i < numNonzeros; ++i)
        stream << "," << _variableNames[nonzeroCoordinates[i]] << "=" << nonzeroValues[i];
      stream << ")";
    }
  }

  std::string Space::vectorToString(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
    mpq_class* nonzeroValues) const
  {
    std::stringstream stream;
    printVector(stream, numNonzeros, nonzeroCoordinates, nonzeroValues);
    return stream.str();
  }

  void Space::printLinearForm(std::ostream& stream, std::size_t numNonzeros,
    std::size_t* nonzeroCoordinates, mpq_class* nonzeroValues) const
  {
    if (numNonzeros == 0)
      stream << "0";
    else
    {
      const mpq_class& x0 = nonzeroValues[0];
      if (x0 == -1)
        stream << "-";
      else if (x0 != 1)
        stream << x0;
      stream << _variableNames[nonzeroCoordinates[0]];
      for (std::size_t i = 1; i < numNonzeros; ++i)
      {
        const mpq_class& xi = nonzeroValues[i];
        if (xi == 1)
          stream << " + ";
        else if (xi == -1)
          stream << " - ";
        else if (xi > 0)
          stream << " + " << xi;
        else
          stream << " - " << (-xi);
        stream << _variableNames[nonzeroCoordinates[i]];
      }
    }
  }

  std::string Space::linearFormToString(std::size_t numNonzeros, std::size_t* nonzeroColumns,
    mpq_class* nonzeroValues) const
  {
    std::stringstream stream;
    printLinearForm(stream, numNonzeros, nonzeroColumns, nonzeroValues);
    return stream.str();
  }

#endif /* IPO_WITH_GMP */

}
