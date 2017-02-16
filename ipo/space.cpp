#include "space.h"

using namespace soplex;

namespace ipo {

  static SpaceData* emptySpace = NULL;

  void freeStaticSpace()
  {
    if (emptySpace != NULL)
    {
      assert(emptySpace->usage() == 1);
      delete emptySpace;
      emptySpace = NULL;
    }
  }

  SpaceData::SpaceData()
    : _usage(0)
  {

  }

  SpaceData::SpaceData(const std::vector< std::string >& variableNames)
    : _variableNames(variableNames), _usage(0)
  {

  }

  SpaceData::~SpaceData()
  {

  }

  bool SpaceData::operator==(const SpaceData& other) const
  {
    return _variableNames == other._variableNames;
  }

  void SpaceData::addVariable(const std::string& name)
  {
#ifdef IPO_DEBUG
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      if (variableName(v) == name)
        throw std::runtime_error("Error while adding variable to SpaceData: Duplicate name.");
    }
#endif
    _variableNames.push_back(name);
  }

  void SpaceData::unmarkUsed()
  {
    _usage--;
    if (_usage == 0)
      delete this;
  }

  Space::Space()
  {
    if (emptySpace == NULL)
    {
      emptySpace = new SpaceData();
      emptySpace->markUsed();
    }
    _data = emptySpace;
    _data->markUsed();
  }

  void Space::printVector(std::ostream& stream, const Vector& vector) const
  {
    if (vector.size() == 0)
      stream << "()";
    else
    {
      assert(vector.index(vector.size() - 1) < dimension());

      stream << "(" << _data->variableName(vector.index(0)) << "=" << vector.value(0);
      for (std::size_t p = 1; p < vector.size(); ++p)
        stream << "," << _data->variableName(vector.index(p)) << "=" << vector.value(p);
      stream << ")";
    }
  }

  std::string Space::vectorToString(const Vector& vector) const
  {
    std::stringstream stream;
    printVector(stream, vector);
    return stream.str();
  }

  void Space::printLinearForm(std::ostream& stream, const Vector& linearForm) const
  {
    if (linearForm.size() == 0)
      stream << "0";
    else
    {
      assert(linearForm.index(linearForm.size() - 1) < dimension());

      const Rational& x0 = linearForm.value(0);
      if (x0 == -1)
        stream << "-";
      else if (x0 != 1)
        stream << x0;
      stream << _data->variableName(linearForm.index(0));
      for (std::size_t p = 1; p < linearForm.size(); ++p)
      {
        const Rational& xp = linearForm.value(p);
        if (xp == 1)
          stream << " + ";
        else if (xp == -1)
          stream << " - ";
        else if (xp > 0)
          stream << " + " << xp;
        else
          stream << " - " << (-xp);
        stream << _data->variableName(linearForm.index(p));
      }
    }
  }

  void Space::printLinearConstraint(std::ostream& stream, const LinearConstraint& constraint) const
  {
    printLinearForm(stream, constraint.normal());
    stream << ' ' << constraint.type() << "= " << constraint.rhs();
  }

  std::string Space::linearConstraintToString(const LinearConstraint& constraint) const
  {
    std::stringstream stream;
    printLinearConstraint(stream, constraint);
    return stream.str();
  }
}
