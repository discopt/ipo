#include "space.h"

using namespace soplex;

namespace ipo {

  Space::Space()
  {

  }

  Space::Space(const std::vector<std::string>& variables) : _variables(variables)
  {

  }

  Space::Space(const Space& other) : _variables(other._variables)
  {

  }

  Space::~Space()
  {

  }

  void Space::addVariable(const std::string& name)
  {
#ifdef IPO_DEBUG
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      if (variable(v) == name)
        throw std::runtime_error("Error while adding variable to Space: Duplicate name.");
    }
#endif
    _variables.push_back(name);
  }

  void Space::printLinearForm(std::ostream& stream, const SVectorRational* coefficients) const
  {
    bool first = true;
    for (std::size_t i = 0; i < coefficients->size(); ++i)
    {
      std::size_t v = coefficients->index(i);
      const Rational& x = coefficients->value(i);
      if (x < 0)
        stream << (first ? "-" : " - ");
      else if (!first)
        stream << " + ";
      if (x != 1 && x != -1)
      {
        if (x > 0)
          stream << x << ' ';
        else
          stream << (-x) << ' ';
      }
      stream << variable(v);
      first = false;
    }
  }

  void Space::printLinearForm(std::ostream& stream, const Vector& coefficients) const
  {
    bool first = true;
    for (std::size_t i = 0; i < coefficients.size(); ++i)
    {
      std::size_t v = coefficients.index(i);
      const Rational& x = coefficients.value(i);
      if (x < 0)
        stream << (first ? "-" : " - ");
      else if (!first)
        stream << " + ";
      if (x != 1 && x != -1)
      {
        if (x > 0)
          stream << x << ' ';
        else
          stream << (-x) << ' ';
      }
      stream << variable(v);
      first = false;
    }
  }


  void Space::printLinearForm(std::ostream& stream, const soplex::VectorRational* coefficients) const
  {
    bool first = true;
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      const Rational& x = (*coefficients)[v];
      if (x == 0)
        continue;
      if (x < 0)
        stream << (first ? "-" : " - ");
      else if (!first)
        stream << " + ";
      if (x != 1 && x != -1)
      {
        if (x > 0)
          stream << x << ' ';
        else
          stream << (-x) << ' ';
      }
      stream << variable(v);
      first = false;
    }
  }

  void Space::printRow(std::ostream& stream, const LPRowRational& row) const
  {
    const Rational& lhs = row.lhs();
    const Rational& rhs = row.rhs();
    printRow(stream, lhs > -infinity/2 ? &lhs : NULL,
      rhs < infinity / 2 ? &rhs : NULL, row.rowVector());
  }

  void Space::printRow(std::ostream& stream, const LPRowSetRational& rows, std::size_t index) const
  {
    assert(index < rows.num());
    const Rational& lhs = rows.lhs(index);
    const Rational& rhs = rows.rhs(index);
    printRow(stream, lhs > -infinity / 2 ? &lhs : NULL,
      rhs < infinity / 2 ? &rhs : NULL, rows.rowVector(index));
  }

  void Space::printRows(std::ostream& stream, const LPRowSetRational& rows) const
  {
    for (int i = 0; i < rows.num(); ++i)
    {
      printRow(stream, rows, i);
      stream << "\n";
    }
  }

  void Space::printRow(std::ostream& stream, const Rational* lhs, const Rational* rhs,
    const SVectorRational& vector) const
  {
    bool equation = lhs != NULL && rhs != NULL && *lhs == *rhs;
    if (lhs != NULL && !equation)
      stream << *lhs << " <= ";
    printLinearForm(stream, &vector);
    if (vector.size() == 0)
      stream << '0';
    if (rhs != NULL)
      stream << (equation ? " == " : " <= ") << *rhs;
  }

  void Space::printVector(std::ostream& stream, const SVectorRational* vector) const
  {
    bool delimit = false;
    for (std::size_t i = 0; i < vector->size(); ++i)
    {
      std::size_t v = vector->index(i);
      const soplex::Rational& x = vector->value(i);
      stream << (delimit ? ", " : "(") << variable(v) << "=" << x;
      delimit = true;
    }
    if (delimit)
      stream << ")";
    else
      stream << "()";
  }

  void Space::printVector(std::ostream& stream, const Vector& vector) const
  {
    bool delimit = false;
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
      std::size_t v = vector.index(i);
      const soplex::Rational& x = vector.value(i);
      stream << (delimit ? ", " : "(") << variable(v) << "=" << x;
      delimit = true;
    }
    if (delimit)
      stream << ")";
    else
      stream << "()";
  }

  bool Space::operator==(const Space& other) const
  {
    if (this == &other)
      return true;
    if (dimension() != other.dimension())
      return false;
    for (std::size_t v = 0; v < dimension(); ++v)
    {
      if (variable(v) != other.variable(v))
        return false;
    }
    return true;
  }

}
