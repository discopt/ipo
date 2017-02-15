#ifndef _SOPLEX_REPRODUCE_H
#define _SOPLEX_REPRODUCE_H

#include <string>

#include <rational.h>
#include <dsvector.h>
#include <soplex.h>

namespace soplex
{

  class ReproSoPlex
  {
  public:
    ReproSoPlex();
    ~ReproSoPlex();

    double intParam(const SoPlex::IntParam param);
    Real realParam(const SoPlex::RealParam param);
    void setIntParam(const SoPlex::IntParam param, const int value);
    void setRealParam(const SoPlex::RealParam param, const Real value);
    void setBoolParam(const SoPlex::BoolParam param, const bool value);
    Rational objValueRational();
    int numRowsRational() const;
    int numColsRational() const;
    int numNonzerosRational() const;
    void changeObjRational(int v, const Rational& x);
    void changeObjRational(const VectorRational& objective);
    void addRowRational(const LPRowRational& lprow);
    void getRowsRational(int start, int end, LPRowSetRational& lprowset) const;
    void addRowsRational(const LPRowSetRational& rows);
    void getColsRational(int start, int end, LPColSetRational& lpcolset) const;
    void addColsRational(const LPColSetRational& lpcols);
    const SVectorRational& rowVectorRational(int i) const;
    const Rational& lhsRational(int i) const;
    const Rational& rhsRational(int i) const;
    LPRowBase<Rational>::Type rowTypeRational(int i) const;
    void changeRowRational(int row, const LPRowRational& lprow);
    void removeRowRational(int i);
    void removeRowRangeRational(int start, int end, int* perm = NULL);
    const Rational& upperRational(int variable) const;
    const Rational& lowerRational(int variable) const;
    void changeBoundsRational(int i, const Rational& lower, const Rational& upper);
    void changeUpperRational(int i, const Rational& upper);
    void changeLowerRational(int i, const Rational& lower);
    void reproduceSolve(const std::string& fileNamePrefix);
    SPxSolver::Status solve();
    void getPrimalRational(VectorRational& vector);
    void getPrimalRayRational(VectorRational& vector);
    void getObjRational(VectorRational& obj) const;
    void getBasis(SPxSolver::VarStatus rows[], SPxSolver::VarStatus cols[]);
    void writeFileRational(const char* fileName, const NameSet* rowNames = 0, const NameSet* colNames = 0, const DIdxSet* intVars = 0);
    void printRational(const Rational& x);

  protected:
    soplex::SoPlex _spx;
    std::stringstream _stream;
    std::size_t _counter;
  };

}

#endif /* _SOPLEX_REPRODUCE_H */
