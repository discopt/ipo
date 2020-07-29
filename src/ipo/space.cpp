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

  void Space::printVector(std::ostream& str, const sparse_vector<double>& vector,
    bool rounded) const
  {
    if (vector.empty())
    {
      str << "()";
      return;
    }
  
    sparse_vector<double>::const_iterator iter = vector.begin();
    assert(iter != vector.end());
    str << '(' << _variableNames[iter->first] << '=' << iter->second;
    ++iter;
    for (; iter != vector.end(); ++iter)
      str << ',' << _variableNames[iter->first] << '=' << iter->second;
    str << ')';
  }

  std::string Space::printVector(const sparse_vector<double>& vector, bool rounded) const
  {
    std::ostringstream str;
    printVector(str, vector, rounded);
    return str.str();
  }

  void Space::printLinearForm(std::ostream& str, const sparse_vector<double>& vector,
    bool rounded) const
  {
    if (vector.empty())
    {
      str << '0';
      return;
    }

    sparse_vector<double>::const_iterator iter = vector.begin();
    assert(iter != vector.end());
    double x0 = iter->second;
    if (x0 == -1)
      str << '-';
    else if (x0 != 1)
      str << x0;
    str << '<' << _variableNames[iter->first] << '>';
    ++iter;

    for (; iter != vector.end(); ++iter)
    {
      double xi = iter->second;
      if (xi == 1)
        str << " + ";
      else if (xi == -1)
        str << " - ";
      else if (xi > 0)
        str << " + " << xi;
      else
        str << " - " << (-xi);
      str << '<' << _variableNames[iter->first] << '>';
    }
  }

  std::string Space::printLinearForm(const sparse_vector<double>& vector, bool rounded) const
  {
    std::stringstream str;
    printLinearForm(str, vector, rounded);
    return str.str();
  }

  void Space::printConstraint(std::ostream& str, const Constraint<double>& constraint,
    bool rounded) const
  {
    switch (constraint.type())
    {
    case EQUATION:
      printLinearForm(str, constraint.vector());
      str << " == " << constraint.rhs();
    break;
    case LESS_OR_EQUAL:
      printLinearForm(str, constraint.vector());
      str << " <= " << constraint.rhs();
    break;
    case GREATER_OR_EQUAL:
      printLinearForm(str, constraint.vector());
      str << " >= " << constraint.lhs();
    break;
    case RANGED:
      str << constraint.lhs() << " <= ";
      printLinearForm(str, constraint.vector());
      str << " <= " << constraint.rhs();
    break;
    }
  }

  std::string Space::printConstraint(const Constraint<double>& constraint, bool rounded) const
  {
    std::stringstream str;
    printConstraint(str, constraint, rounded);
    return str.str();
  }

#if defined(IPO_WITH_GMP)

  void Space::printVector(std::ostream& str, const sparse_vector<mpq_class>& vector,
    bool rounded) const
  {
    if (vector.empty())
    {
      str << "()";
      return;
    }
  
    sparse_vector<mpq_class>::const_iterator iter = vector.begin();
    assert(iter != vector.end());

    if (rounded)
    {
      str << '(' << _variableNames[iter->first] << '=' << iter->second.get_d();
      ++iter;
      for (; iter != vector.end(); ++iter)
        str << ',' << _variableNames[iter->first] << '=' << iter->second.get_d();
      str << ')';
    }
    else
    {
      str << '(' << _variableNames[iter->first] << '=' << iter->second;
      ++iter;
      for (; iter != vector.end(); ++iter)
        str << ',' << _variableNames[iter->first] << '=' << iter->second;
      str << ')';
    }
  }   

  std::string Space::printVector(const sparse_vector<mpq_class>& vector, bool rounded) const
  {
    std::stringstream str;
    printVector(str, vector, rounded);
    return str.str();
  }

  void Space::printLinearForm(std::ostream& str, const sparse_vector<mpq_class>& vector,
      bool rounded) const
  {
    if (vector.empty())
    {
      str << '0';
      return;
    }

    sparse_vector<mpq_class>::const_iterator iter = vector.begin();
    assert(iter != vector.end());
    const mpq_class& x0 = iter->second;
    if (x0 == -1)
      str << '-';
    else if (x0 != 1)
    {
      if (rounded)
        str << x0.get_d();
      else
        str << x0;
    }
    str << '<' << _variableNames[iter->first] << '>';
    ++iter;

    if (rounded)
    {
      for (; iter != vector.end(); ++iter)
      {
        double xi = iter->second.get_d();
        if (xi == 1)
          str << " + ";
        else if (xi == -1)
          str << " - ";
        else if (xi > 0)
          str << " + " << xi;
        else
          str << " - " << (-xi);
        str << '<' << _variableNames[iter->first] << '>';
      }
    }
    else
    {
      for (; iter != vector.end(); ++iter)
      {
        const mpq_class& xi = iter->second;
        if (xi == 1)
          str << " + ";
        else if (xi == -1)
          str << " - ";
        else if (xi > 0)
          str << " + " << xi;
        else
          str << " - " << (-xi);
        str << '<' << _variableNames[iter->first] << '>';
      }
    }
  }

  std::string Space::printLinearForm(const sparse_vector<mpq_class>& vector, bool rounded) const
  {
    std::stringstream str;
    printLinearForm(str, vector, rounded);
    return str.str();
  }

  void Space::printConstraint(std::ostream& str, const Constraint<mpq_class>& constraint,
    bool rounded) const
  {
    switch (constraint.type())
    {
    case EQUATION:
      printLinearForm(str, constraint.vector(), rounded);
      str << " == " << constraint.rhs();
    break;
    case LESS_OR_EQUAL:
      printLinearForm(str, constraint.vector(), rounded);
      str << " <= " << constraint.rhs();
    break;
    case GREATER_OR_EQUAL:
      printLinearForm(str, constraint.vector(), rounded);
      str << " >= " << constraint.lhs();
    break;
    case RANGED:
      str << constraint.lhs() << " <= ";
      printLinearForm(str, constraint.vector(), rounded);
      str << " <= " << constraint.rhs();
    break;
    }
  }

  std::string Space::printConstraint(const Constraint<mpq_class>& constraint, bool rounded) const
  {
    std::stringstream str;
    printConstraint(str, constraint, rounded);
    return str.str();
  }

#endif /* IPO_WITH_GMP */

}
