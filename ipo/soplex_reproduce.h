#ifndef _SOPLEX_REPRODUCE_H
#define _SOPLEX_REPRODUCE_H

#include <string>

#include <rational.h>
#include <dsvector.h>
#include <soplex.h>

namespace soplex
{

class ReproSoPlex : public soplex::SoPlex
{
public:
  ReproSoPlex()
  {

  }

  ~ReproSoPlex()
  {

  }

  double intParam(const SoPlex::IntParam param)
  {
    return _spx.intParam(param);
  }

  Real realParam(const SoPlex::RealParam param)
  {
    return _spx.realParam(param);
  }

  void setIntParam(const SoPlex::IntParam param, const int value)
  {
    _stream << "spx.setIntParam(SoPlex::IntParam(" << int(param) << "), " << value << ");\n";
    _spx.setIntParam(param, value);
  }

  void setRealParam(const SoPlex::RealParam param, const Real value)
  {
    _stream << "spx.setRealParam(SoPlex::RealParam(" << int(param) << "), " << value << ");\n";
    _spx.setRealParam(param, value);
  }

  void setBoolParam(const SoPlex::BoolParam param, const bool value)
  {
    _stream << "spx.setBoolParam(SoPlex::BoolParam(" << int(param) << "), " << (value ? "true" : "false") << ");\n";
    _spx.setBoolParam(param, value);
  }

  Rational objValueRational()
  {
    return _spx.objValueRational();
  }

  int numRowsRational() const
  {
    return _spx.numRowsRational();
  }

  int numColsRational() const
  {
    return _spx.numColsRational();
  }

  int numNonzerosRational() const
  {
    return _spx.numNonzerosRational();
  }

  void changeObjRational(const VectorRational& objective)
  {
    _stream << "{\n";
    _stream << "  DVectorRational objectiveVector;\n";
    _stream << "  objectiveVector.reDim(" << objective.dim() << ");\n";
    for (std::size_t v = 0; v < objective.dim(); ++v)
    {
      _stream << "  objectiveVector[" << v << "] = ";
      printRational(objective[v]);
      _stream << ";\n";
    }
    _stream << "  spx.changeObjRational(objectiveVector);\n";
    _stream << "}\n";
    _spx.changeObjRational(objective);
  }

  void addRowRational(const LPRowRational& lprow)
  {
    _stream << "{\n";
    _stream << "  DSVectorRational vector;\n";
    for (int p = 0; p < lprow.rowVector().size(); ++p)
    {
      _stream << "  vector.add(" << lprow.rowVector().index(p) << ", ";
      printRational(lprow.rowVector().value(p));
      _stream << ");\n";
    }
    _stream << "  spx.addRowRational(LPRowRational(";
    printRational(lprow.lhs());
    _stream << ", vector, ";
    printRational(lprow.rhs());
    _stream << "));\n";
    _stream << "}\n";
    _spx.addRowRational(lprow);
  }

  void changeRowRational(int row, const LPRowRational& lprow)
  {
    _stream << "{\n";
    _stream << "  DSVectorRational vector;\n";
    for (int p = 0; p < lprow.rowVector().size(); ++p)
    {
      _stream << "  vector.add(" << lprow.rowVector().index(p) << ", ";
      printRational(lprow.rowVector().value(p));
      _stream << ");\n";
    }
    _stream << "  spx.changeRowRational(" << row << ", LPRowRational(";
    printRational(lprow.lhs());
    _stream << ", vector, ";
    printRational(lprow.rhs());
    _stream << ");\n";
    _stream << "}\n";
    _spx.changeRowRational(row, lprow);
  }

  void reproduceSolve(const std::string& fileName)
  {
    std::ofstream file(fileName.c_str(), std::ios::out);
    file << "#include <soplex.h>\n";
    file << "\n";
    file << "using namespace soplex;\n";
    file << "\n";
    file << "int main(int argc, char** argv)\n";
    file << "{\n";
    file << "  SoPlex spx;\n";
    file << "  spx.printVersion();\n";
    file << _stream.str();
    file << "  spx.solve();\n";
    file << "}\n";
    file.close();
  }

  SPxSolver::Status solve()
  {
    _stream << "spx.solve();\n";
    return _spx.solve();
  }

  void getPrimalRational(VectorRational& solution)
  {
    _spx.getPrimalRational(solution);
  }

  void getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[])
  {
    _spx.getBasis(rows, cols);
  }

  void writeFileRational(const char* fileName, const NameSet* rowNames = 0, const NameSet* colNames = 0,
    const DIdxSet* intVars = 0)
  {
    _spx.writeFileRational(fileName, rowNames, colNames, intVars);
  }

protected:
  void printRational(const Rational& x)
  {
    if (x == -infinity)
      _stream << "-infinity";
    else if (x == infinity)
      _stream << "infinity";
    else
    {
      _stream << "Rational(" << double(x) << ")";
    }
  }

protected:
  soplex::SoPlex _spx;
  std::stringstream _stream;
};


static void DUMP_printRational(const soplex::Rational& x)
{
  if (x <= -soplex::infinity)
    std::cout << "-soplex::infinity";
  else if (x >= soplex::infinity)
    std::cout << "soplex::infinity";
  else
  {
    std::cout << "soplex::Rational(" << double(x) << ")";
  }
}


static void DUMP_addRowRationalLPRowRational(const soplex::Rational& lhs, const soplex::SVectorRational& vector, const
soplex::Rational& rhs, const std::string& spxAccess = "spx.", const std::string& prefix = "")
{
  std::cout << prefix << "{\n";
  std::cout << prefix << "  soplex::DSVectorRational rowVector;\n";
  for (int p = 0; p < vector.size(); ++p)
  {
    std::cout << prefix << "  rowVector.add(" << vector.index(p) << ", ";
    DUMP_printRational(vector.value(p));
    std::cout << ");\n";
  }
  std::cout << prefix << "  " << spxAccess << "addRowRational(soplex::LPRowRational(";
  DUMP_printRational(lhs);
  std::cout << ", rowVector, ";
  DUMP_printRational(rhs);
  std::cout << "));\n";
  std::cout << prefix << "}\n";
}

static void DUMP_changeObjRational(const soplex::VectorRational& objective, const std::string& spxAccess = "spx.",
  const std::string& prefix = "")
{
  std::cout << prefix << "{\n";
  std::cout << prefix << "  soplex::DVectorRational vector;\n";
  std::cout << prefix << "  vector.reDim(" << objective.dim() << ");\n";
  for (std::size_t v = 0; v < objective.dim(); ++v)
    std::cout << prefix << "  vector[" << v << "] = " << double(objective[v]) << ";\n";
  std::cout << prefix << "  " << spxAccess << "changeObjRational(vector);\n";
  std::cout << prefix << "}\n";
}

}

#endif /* _SOPLEX_REPRODUCE_H */
