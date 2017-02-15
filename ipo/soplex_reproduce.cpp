#include "soplex_reproduce.h"

namespace soplex
{

  ReproSoPlex::ReproSoPlex()
    : _counter(0)
  {
    
  }

  ReproSoPlex::~ReproSoPlex()
  {

  }

  double ReproSoPlex::intParam(const SoPlex::IntParam param)
  {
    return _spx.intParam(param);
  }

  Real ReproSoPlex::realParam(const SoPlex::RealParam param)
  {
    return _spx.realParam(param);
  }

  void ReproSoPlex::setIntParam(const SoPlex::IntParam param, const int value)
  {
    _stream << "spx.setIntParam(SoPlex::IntParam(" << int(param) << "), " << value << ");\n";
    _spx.setIntParam(param, value);
  }

  void ReproSoPlex::setRealParam(const SoPlex::RealParam param, const Real value)
  {
    _stream << "spx.setRealParam(SoPlex::RealParam(" << int(param) << "), " << value << ");\n";
    _spx.setRealParam(param, value);
  }

  void ReproSoPlex::setBoolParam(const SoPlex::BoolParam param, const bool value)
  {
    _stream << "spx.setBoolParam(SoPlex::BoolParam(" << int(param) << "), " << (value ? "true" : "false") << ");\n";
    _spx.setBoolParam(param, value);
  }

  Rational ReproSoPlex::objValueRational()
  {
    return _spx.objValueRational();
  }

  int ReproSoPlex::numRowsRational() const
  {
    return _spx.numRowsRational();
  }

  int ReproSoPlex::numColsRational() const
  {
    return _spx.numColsRational();
  }

  int ReproSoPlex::numNonzerosRational() const
  {
    return _spx.numNonzerosRational();
  }

  void ReproSoPlex::changeObjRational(int v, const Rational& x)
  {
    _stream << "spx.changeObjRational(" << v << ",";
    printRational(x);
    _stream << ");\n";
    _spx.changeObjRational(v, x);
  }

  void ReproSoPlex::changeObjRational(const VectorRational& objective)
  {
    _stream << "{\n";
    _stream << "  DVectorRational objectiveVector;\n";
    _stream << "  objectiveVector.reDim(" << objective.dim() << ");\n";
    for (int v = 0; v < objective.dim(); ++v)
    {
      _stream << "  objectiveVector[" << v << "] = ";
      printRational(objective[v]);
      _stream << ";\n";
    }
    _stream << "  spx.changeObjRational(objectiveVector);\n";
    _stream << "}\n";
    _spx.changeObjRational(objective);
  }

  void ReproSoPlex::addRowRational(const LPRowRational& lprow)
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

  void ReproSoPlex::getRowsRational(int start, int end, LPRowSetRational& lprowset) const
  {
    _spx.getRowsRational(start, end, lprowset);
  }

  void ReproSoPlex::addRowsRational(const LPRowSetRational& rows)
  {
    _stream << "{\n";
    _stream << "  LPRowSetRational rows;\n";
    for (int r = 0; r < rows.num(); ++r)
    {
      _stream << "  {\n";
      _stream << "    DSVectorRational vector;\n";
      for (int p = 0; p < rows.rowVector(r).size(); ++p)
      {
        _stream << "    vector.add(" << rows.rowVector(r).index(p) << ", ";
        printRational(rows.rowVector(r).value(p));
        _stream << ");\n";
      }
      _stream << "    rows.add(LPRowRational(";
      printRational(rows.lhs(r));
      _stream << ", vector, ";
      printRational(rows.rhs(r));
      _stream << ", ";
      printRational(rows.obj(r));
      _stream << "));\n";
      _stream << "  }\n";
    }
    _stream << "  spx.addRowsRational(rows);\n";
    _stream << "}\n";
    _spx.addRowsRational(rows);
  }

  void ReproSoPlex::getColsRational(int start, int end, LPColSetRational& lpcolset) const
  {
    _spx.getColsRational(start, end, lpcolset);
  }

  void ReproSoPlex::addColsRational(const LPColSetRational& lpcols)
  {
    _stream << "{\n";
    _stream << "  LPColSetRational cols;\n";
    for (int c = 0; c < lpcols.num(); ++c)
    {
      _stream << "  {\n";
      _stream << "    DSVectorRational vector;\n";
      for (int p = 0; p < lpcols.colVector(c).size(); ++p)
      {
        _stream << "    vector.add(" << lpcols.colVector(c).index(p) << ", ";
        printRational(lpcols.colVector(c).value(p));
        _stream << ");\n";
      }
      _stream << "    cols.add(";
      printRational(lpcols.maxObj(c));
      _stream << ", ";
      printRational(lpcols.lower(c));
      _stream << ", vector, ";
      printRational(lpcols.upper(c));
      _stream << ");\n";
      _stream << "  }\n";
    }
    _stream << "  spx.addColsRational(cols);\n";
    _stream << "}\n";
    _spx.addColsRational(lpcols);
  }

  const SVectorRational& ReproSoPlex::rowVectorRational(int i) const
  {
    return _spx.rowVectorRational(i);
  }

  const Rational& ReproSoPlex::lhsRational(int i) const
  {
    return _spx.lhsRational(i);
  }

  const Rational& ReproSoPlex::rhsRational(int i) const
  {
    return _spx.rhsRational(i);
  }

  LPRowBase<Rational>::Type ReproSoPlex::rowTypeRational(int i) const
  {
    return _spx.rowTypeRational(i);
  }

  void ReproSoPlex::changeRowRational(int row, const LPRowRational& lprow)
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

  void ReproSoPlex::removeRowRational(int i)
  {
    _spx.removeRowRational(i);
  }

  void ReproSoPlex::removeRowRangeRational(int start, int end, int* perm)
  {
    _spx.removeRowRangeRational(start, end, perm);
  }

  const Rational& ReproSoPlex::upperRational(int variable) const
  {
    return _spx.upperRational(variable);
  }

  const Rational& ReproSoPlex::lowerRational(int variable) const
  {
    return _spx.lowerRational(variable);
  }

  void ReproSoPlex::changeBoundsRational(int i, const Rational& lower, const Rational& upper)
  {
    _stream << "{\n";
    _stream << "  spx.changeBoundsRational(" << i << ", ";
    printRational(lower);
    _stream << ", ";
    printRational(upper);
    _stream << ");\n";
    _stream << "}\n";
    _spx.changeBoundsRational(i, lower, upper);
  }
  
  void ReproSoPlex::changeUpperRational(int i, const Rational& upper)
  {
    _stream << "{\n";
    _stream << "  spx.changeUpperRational(" << i << ", ";
    printRational(upper);
    _stream << ");\n";
    _stream << "}\n";
    _spx.changeUpperRational(i, upper);
  }

  void ReproSoPlex::changeLowerRational(int i, const Rational& lower)
  {
    _stream << "{\n";
    _stream << "  spx.changeLowerRational(" << i << ", ";
    printRational(lower);
    _stream << ");\n";
    _stream << "}\n";
    _spx.changeUpperRational(i, lower);
  }

  void ReproSoPlex::reproduceSolve(const std::string& fileNamePrefix)
  {
    std::stringstream fileName;
    fileName << fileNamePrefix << "-" << std::setfill('0') << std::setw(6) << _counter << ".cpp";
    _counter++;

    std::ofstream file(fileName.str().c_str(), std::ios::out);
    file << "#include <soplex.h>\n";
    file << "#include <rational.h>\n";
    file << "\n";
    file << "using namespace soplex;\n";
    file << "\n";
    file << "int main(int argc, char** argv)\n";
    file << "{\n";
    file << "  SoPlex spx;\n";
    file << "  spx.printVersion();\n";
    file << _stream.str();
    file << "std::cerr << \"!!! solve() !!!\" << std::endl;\n";
    file << "spx.solve();\n";
    file << "}\n";
    file.close();
  }

  SPxSolver::Status ReproSoPlex::solve()
  {
    _stream << "std::cerr << \"!!! solve() !!!\" << std::endl;\n";
    _stream << "spx.solve();\n";
    return _spx.solve();
  }

  void ReproSoPlex::getPrimalRational(VectorRational& vector)
  {
    _spx.getPrimalRational(vector);
  }

  void ReproSoPlex::getPrimalRayRational(VectorRational& vector)
  {
    _spx.getPrimalRayRational(vector);
  }

  void ReproSoPlex::getObjRational(VectorRational& obj) const
  {
    _spx.getObjRational(obj);
  }

  void ReproSoPlex::getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[])
  {
    _spx.getBasis(rows, cols);
  }

  void ReproSoPlex::writeFileRational(const char* fileName, const NameSet* rowNames, const NameSet* colNames, const DIdxSet* intVars)
  {
    _spx.writeFileRational(fileName, rowNames, colNames, intVars);
  }

  void ReproSoPlex::printRational(const Rational& x)
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


}
