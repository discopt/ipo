#include <ipo/lp.hpp>

#include <cmath>
#include <string>


#if defined(IPO_WITH_SOPLEX)
#define SOPLEX_WITH_GMP
#include <soplex.h>
#endif /* IPO_WITH_SOPLEX */


namespace ipo
{
  std::ostream& operator<<(std::ostream& stream, LPStatus status)
  {
    switch (status)
    {
    case OPTIMAL:
      return stream << "optimal";
    case UNBOUNDED:
      return stream << "unbounded";
    case INFEASIBLE:
      return stream << "infeasible";
    case TIME_LIMIT:
      return stream << "time limit reached";
    case ITERATION_LIMIT:
      return stream << "iteration limit reached";
    case CYCLING:
      return stream << "cycling detected";
    case NUMERICS:
      return stream << "encountered numerical difficulties";
    default:
      return stream << "Unknown status <" << (int)status << ">";
    }
  }

  template <bool Removable>
  class ElementMapping
  {
  };

  template <>
  class ElementMapping<true>
  {
  public:
    ElementMapping()
    {

    }

    std::size_t externalToInternal(std::size_t external) const
    {
      assert(external < _externalToInternal.size());
      assert(_externalToInternal[external] < std::numeric_limits<std::size_t>::max());
      return _externalToInternal[external];
    }

    std::size_t internalToExternal(std::size_t internal) const
    {
      assert(internal < _internalToExternal.size()); 
      return _internalToExternal[internal];
    }

    std::size_t add(std::size_t internal, std::size_t external)
    {
      assert(internal == _internalToExternal.size());
      if (external == std::numeric_limits<std::size_t>::max())
        external = _externalToInternal.size();
      assert(external >= _externalToInternal.size()
        || _externalToInternal[external] == std::numeric_limits<std::size_t>::max());
      if (external >= _externalToInternal.size())
        _externalToInternal.resize(external + 1, std::numeric_limits<std::size_t>::max());
      else
        assert(_externalToInternal[external] == std::numeric_limits<std::size_t>::max());
      _externalToInternal[external] = internal;
      _internalToExternal.push_back(external);
      return external;
    }

  protected:
    std::vector<std::size_t> _externalToInternal;
    std::vector<std::size_t> _internalToExternal;
  };

  template <>
  class ElementMapping<false>
  {
  public:
    ElementMapping()
    {

    }

    std::size_t externalToInternal(std::size_t external) const
    {
      return external;
    }

    std::size_t internalToExternal(std::size_t internal) const
    {
      return internal;
    }

    std::size_t add(std::size_t internal, std::size_t external)
    { 
      return internal;
    }
  };

#if defined(IPO_WITH_SOPLEX)
  static
  bool createNameSet(soplex::NameSet& nameSet, const std::vector<std::string>& internalNames)
  {
    for (std::size_t i = 0; i < internalNames.size(); ++i)
    {
      soplex::DataKey key;
      nameSet.add(key, internalNames[i].c_str());
      if (!key.isValid())
        return false;
    }
    return true;
  }

#endif /* IPO_WITH_SOPLEX */

#if !defined(IPO_WITH_GUROBI) && defined(IPO_WITH_SOPLEX)

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  class LPImplementation;

  template <bool RowsRemovable, bool ColumnsRemovable>
  class LPImplementation<double, RowsRemovable, ColumnsRemovable>
  {
  public:
    typedef double Number;
    typedef ElementMapping<RowsRemovable> RowMapping;
    typedef ElementMapping<ColumnsRemovable> ColumnMapping;

    LPImplementation()
    {
      _spx.setIntParam(soplex::SoPlex::SIMPLIFIER, soplex::SoPlex::SIMPLIFIER_AUTO);
      _spx.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);
    }

    ~LPImplementation()
    {

    }

    std::size_t numRows() const
    {
      return _spx.numRows();

    }

    std::size_t numColumns() const
    {
      return _spx.numCols();
    }

    LPStatus status() const
    {
      switch (_spx.status())
      {
        case soplex::SPxSolver::SINGULAR:
          return NUMERICS;
        case soplex::SPxSolver::ABORT_CYCLING:
          return CYCLING;
        case soplex::SPxSolver::ABORT_TIME:
          return TIME_LIMIT;
        case soplex::SPxSolver::ABORT_ITER:
          return ITERATION_LIMIT;
        case soplex::SPxSolver::OPTIMAL:
          return OPTIMAL;
        case soplex::SPxSolver::UNBOUNDED:
          return UNBOUNDED;
        case soplex::SPxSolver::INFEASIBLE:
          return INFEASIBLE;
      default:
        std::ostringstream stream;
        stream << "Unhandled SoPlex status " << _spx.status();
        throw std::runtime_error(stream.str());
      }
    }

    double getObjectiveValue() const
    {
      return const_cast<LPImplementation*>(this)->_spx.objValueReal();
    }

    bool hasPrimalSolution() const
    {
      return !_primalSolution.empty();
    }

    double getPrimalValue(LPColumn column) const
    {
      size_t j = _columnMapping.externalToInternal(column);
      assert(j < _primalSolution.size());
      return _primalSolution[j];
    }

    std::vector<double> getPrimalSolution(const std::vector<LPColumn>& columns) const
    {
      std::vector<double> solution;
      if (columns.empty())
      {
        for (std::size_t j = 0; j < numColumns(); ++j)
        {
          LPColumn external = _columnMapping.internalToExternal(j);
          if (external >= solution.size())
            solution.resize(external + 1, 0);
          solution[external] = _primalSolution[j];
        }
      }
      else
      {
        solution.reserve(columns.size());
        for (const auto column : columns)
          solution.push_back(getPrimalValue(column));
      }
      return solution;
    }

    bool hasDualSolution() const
    {
      return !_dualSolution.empty();
    }

    double getDualValue(LPRow row) const
    {
      size_t i = _rowMapping.externalToInternal(row);
      assert(i < _dualSolution.size());
      return _dualSolution[i];
    }

    void setSense(LPSense newSense)
    {
      assert(newSense == LPSense::MAXIMIZE || newSense == LPSense::MINIMIZE);
      _spx.setIntParam(soplex::SoPlex::OBJSENSE,
        newSense == LPSense::MAXIMIZE ? soplex::SoPlex::OBJSENSE_MAXIMIZE : soplex::SoPlex::OBJSENSE_MINIMIZE);
    }

    LPColumn addColumn(
      double lowerBound,
      double upperBound,
      double objectiveCoefficient,
      const std::string& name,
      LPColumn unusedColumn)
    {
      std::size_t internalColumn = numColumns();
      _vector.clear();
      _spx.addColReal(soplex::LPColBase<double>(objectiveCoefficient, _vector,
        std::isnormal(upperBound) ? upperBound : soplex::infinity,
        std::isnormal(lowerBound) ? lowerBound : -soplex::infinity));
      _columnNames.push_back(name);
      return _columnMapping.add(internalColumn, unusedColumn);
    }

    LPRow addRow(
      double lhs,
      std::size_t numNonzeros,
      const LPColumn* nonzeroVariables,
      const double* nonzeroCoefficients,
      double rhs,
      const std::string& name,
      LPRow unusedRow)
    {
      std::size_t internalRow = numRows();
      _vector.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
      {
        std::size_t j = _columnMapping.externalToInternal(nonzeroVariables[i]);
        assert(j < numColumns());
        _vector.add(j, nonzeroCoefficients[i]);
      }
      _spx.addRowReal(soplex::LPRowBase<double>(lhs, _vector, rhs));
      _rowNames.push_back(name);
      return _rowMapping.add(internalRow, unusedRow);
    }

    void update()
    {

    }

    void changeUpper(LPColumn column, double newUpperBound)
    {
      _spx.changeUpperReal(_columnMapping.externalToInternal(column), newUpperBound);
    }

    void changeLower(LPColumn column, double newLowerBound)
    {
      _spx.changeLowerReal(_columnMapping.externalToInternal(column), newLowerBound);
    }

    void changeBounds(LPColumn column, double newLowerBound, double newUpperBound)
    {
      _spx.changeBoundsReal(_columnMapping.externalToInternal(column), newLowerBound, newUpperBound);
    }

    void changeObjective(LPColumn column, double newObjectiveCoefficient)
    {
      _spx.changeObjReal(_columnMapping.externalToInternal(column), newObjectiveCoefficient);
    }

    void write(const std::string& fileName) const
    {
      soplex::NameSet rowNames, columnNames;;
      bool validRowNames = createNameSet(rowNames, _rowNames);
      bool validColumnNames = createNameSet(columnNames, _columnNames);
      _spx.writeFileReal(fileName.c_str(), validRowNames ? &rowNames : NULL, validColumnNames ? &columnNames : NULL);
    }

    LPStatus solve()
    {
      _spx.optimize();
      if (_spx.isPrimalFeasible())
      {
        _primalSolution.resize(numColumns());
        _spx.getPrimalReal(&_primalSolution[0], numColumns());
      }
      else
        _primalSolution.clear();
      if (_spx.isDualFeasible())
      {
        _dualSolution.resize(numRows());
        _spx.getDualReal(&_dualSolution[0], numRows());
      }
      else
        _dualSolution.clear();

      return status();
    }

  private:
    RowMapping _rowMapping;
    ColumnMapping _columnMapping;
    std::vector<double> _primalSolution;
    std::vector<double> _primalRay;
    std::vector<double> _dualSolution;
    std::vector<double> _dualRay;

    soplex::SoPlex _spx;
    soplex::DSVectorReal _vector;
    std::vector<std::string> _rowNames;
    std::vector<std::string> _columnNames;
  };

#if defined(IPO_WITH_SOPLEX) && defined(IPO_WITH_GMP)
  
  template <bool RowsRemovable, bool ColumnsRemovable>
  class LPImplementation<mpq_class, RowsRemovable, ColumnsRemovable>
  {
  public:
    typedef mpq_class Number;
    typedef ElementMapping<RowsRemovable> RowMapping;
    typedef ElementMapping<ColumnsRemovable> ColumnMapping;

  private:
    soplex::Rational mpq2rational(const mpq_class& number)
    {
      mpq_t x;
      mpq_init(x);
      mpq_set(x, number.get_mpq_t());
      soplex::Rational result(x);
      mpq_clear(x);
      return result;
    }

    mpq_class rational2mpq(const soplex::Rational& number)
    {
      return mpq_class(number.getMpqRef());
    }

  public:
    LPImplementation()
    {
      _spx.setIntParam(soplex::SoPlex::SOLVEMODE, soplex::SoPlex::SOLVEMODE_RATIONAL);
      _spx.setIntParam(soplex::SoPlex::SYNCMODE, soplex::SoPlex::SYNCMODE_AUTO);
      _spx.setIntParam(soplex::SoPlex::SIMPLIFIER, soplex::SoPlex::SIMPLIFIER_AUTO);
      _spx.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_ERROR);
      _spx.setRealParam(soplex::SoPlex::FEASTOL, 0.0);
      _spx.setBoolParam(soplex::SoPlex::RATREC, true);
      _spx.setBoolParam(soplex::SoPlex::RATFAC, true);
    }

    ~LPImplementation()
    {

    }

    std::size_t numRows() const
    {
      return _spx.numRowsRational();

    }

    std::size_t numColumns() const
    {
      return _spx.numColsRational();
    }

    LPStatus status() const
    {
      switch (_spx.status())
      {
        case soplex::SPxSolver::SINGULAR:
          return NUMERICS;
        case soplex::SPxSolver::ABORT_CYCLING:
          return CYCLING;
        case soplex::SPxSolver::ABORT_TIME:
          return TIME_LIMIT;
        case soplex::SPxSolver::ABORT_ITER:
          return ITERATION_LIMIT;
        case soplex::SPxSolver::OPTIMAL:
          return OPTIMAL;
        case soplex::SPxSolver::UNBOUNDED:
          return UNBOUNDED;
        case soplex::SPxSolver::INFEASIBLE:
          return INFEASIBLE;
      default:
        std::ostringstream stream;
        stream << "Unhandled SoPlex status " << _spx.status();
        throw std::runtime_error(stream.str());
      }
    }

    double getObjectiveValue() const
    {
      return const_cast<LPImplementation*>(this)->_spx.objValueRational();
    }

    bool hasPrimalSolution() const
    {
      return !_primalSolution.empty();
    }

    const mpq_class& getPrimalValue(LPColumn column) const
    {
      size_t j = _columnMapping.externalToInternal(column);
      assert(j < _primalSolution.size());
      return _primalSolution[j];
    }

    std::vector<mpq_class> getPrimalSolution(const std::vector<LPColumn>& columns) const
    {
      std::vector<mpq_class> solution;
      if (columns.empty())
      {
        for (std::size_t j = 0; j < numColumns(); ++j)
        {
          LPColumn external = _columnMapping.internalToExternal(j);
          if (external >= solution.size())
            solution.resize(external + 1, 0);
          solution[external] = _primalSolution[j];
        }
      }
      else
      {
        solution.reserve(columns.size());
        for (const auto column : columns)
          solution.push_back(getPrimalValue(column));
      }
      return solution;
    }

    bool hasDualSolution() const
    {
      return !_dualSolution.empty();
    }

    const mpq_class& getDualValue(LPRow row) const
    {
      size_t i = _rowMapping.externalToInternal(row);
      assert(i < _dualSolution.size());
      return _dualSolution[i];
    }

    void setSense(LPSense newSense)
    {
      assert(newSense == LPSense::MAXIMIZE || newSense == LPSense::MINIMIZE);
      _spx.setIntParam(soplex::SoPlex::OBJSENSE,
        newSense == LPSense::MAXIMIZE ? soplex::SoPlex::OBJSENSE_MAXIMIZE : soplex::SoPlex::OBJSENSE_MINIMIZE);
    }

    LPColumn addColumn(
      const mpq_class& lowerBound,
      const mpq_class& upperBound,
      const mpq_class& objectiveCoefficient,
      const std::string& name,
      LPColumn unusedColumn)
    {
      std::size_t internalColumn = numColumns();
      _sparse.clear();
      _spx.addColRational(soplex::LPColBase<soplex::Rational>(mpq2rational(objectiveCoefficient), _sparse,
        mpq2rational(upperBound), mpq2rational(lowerBound)));
      _columnNames.push_back(name);
      return _columnMapping.add(internalColumn, unusedColumn);
    }

    LPRow addRow(
      const mpq_class& lhs,
      std::size_t numNonzeros,
      const LPColumn* nonzeroVariables,
      const mpq_class* nonzeroCoefficients,
      const mpq_class& rhs,
      const std::string& name,
      LPRow unusedRow)
    {
      std::size_t internalRow = numRows();
      _sparse.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
      {
        std::size_t j = _columnMapping.externalToInternal(nonzeroVariables[i]);
        assert(j < numColumns());
        _sparse.add(j, mpq2rational(nonzeroCoefficients[i]));
      }
      _spx.addRowRational(soplex::LPRowBase<soplex::Rational>(mpq2rational(lhs), _sparse, mpq2rational(rhs)));
      _rowNames.push_back(name);
      return _rowMapping.add(internalRow, unusedRow);
    }

    void update()
    {

    }

    void changeUpper(LPColumn column, const mpq_class& newUpperBound)
    {
      _spx.changeUpperRational(_columnMapping.externalToInternal(column), mpq2rational(newUpperBound));
    }

    void changeLower(LPColumn column, const mpq_class& newLowerBound)
    {
      _spx.changeLowerRational(_columnMapping.externalToInternal(column), mpq2rational(newLowerBound));
    }

    void changeBounds(LPColumn column, const mpq_class& newLowerBound, const mpq_class& newUpperBound)
    {
      _spx.changeBoundsRational(_columnMapping.externalToInternal(column), mpq2rational(newLowerBound),
        mpq2rational(newUpperBound));
    }

    void changeObjective(LPColumn column, const mpq_class& newObjectiveCoefficient)
    {
      _spx.changeObjRational(_columnMapping.externalToInternal(column), mpq2rational(newObjectiveCoefficient));
    }

    void write(const std::string& fileName) const
    {
      soplex::NameSet rowNames, columnNames;;
      bool validRowNames = createNameSet(rowNames, _rowNames);
      bool validColumnNames = createNameSet(columnNames, _columnNames);
      _spx.writeFileRational(fileName.c_str(), validRowNames ? &rowNames : NULL, validColumnNames ? &columnNames : NULL);
    }

    LPStatus solve()
    {
      _spx.optimize();
      if (_spx.isPrimalFeasible())
      {
        _dense.reDim(numColumns());
        _spx.getPrimalRational(_dense);
        _primalSolution.resize(numColumns());
        for (std::size_t c = 0; c < _primalSolution.size(); ++c)
          _primalSolution[c] = rational2mpq(_dense[c]);
      }
      else
        _primalSolution.clear();
      if (_spx.isDualFeasible())
      {
        _dense.reDim(numRows());
        _spx.getPrimalRational(_dense);
        _dualSolution.resize(numRows());
        for (std::size_t r = 0; r < _dualSolution.size(); ++r)
          _dualSolution[r] = rational2mpq(_dense[r]);
      }
      else
        _dualSolution.clear();

      return status();
    }

  private:
    RowMapping _rowMapping;
    ColumnMapping _columnMapping;
    std::vector<mpq_class> _primalSolution;
    std::vector<mpq_class> _primalRay;
    std::vector<mpq_class> _dualSolution;
    std::vector<mpq_class> _dualRay;

    soplex::SoPlex _spx;
    soplex::DVectorRational _dense;
    soplex::DSVectorRational _sparse;
    std::vector<std::string> _rowNames;
    std::vector<std::string> _columnNames;
  };

#endif /* IPO_WITH_SOPLEX && IPO_WITH_GMP */

  static double plusInfinityDouble = soplex::infinity;
  static double minusInfinityDouble = -soplex::infinity;
  
  template <>
  const double& LP<double>::plusInfinity()
  {
    return plusInfinityDouble;
  };

  template <>
  const double& LP<double>::minusInfinity()
  {
    return minusInfinityDouble;
  };

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  LP<NumberType, RowsRemovable, ColumnsRemovable>::LP()
  {
    _implementation = new LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  LP<NumberType, RowsRemovable, ColumnsRemovable>::~LP()
  {
    delete static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  std::size_t LP<NumberType, RowsRemovable, ColumnsRemovable>::numRows() const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->numRows();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  std::size_t LP<NumberType, RowsRemovable, ColumnsRemovable>::numColumns() const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->numColumns();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  LPStatus LP<NumberType, RowsRemovable, ColumnsRemovable>::status() const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->status();
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  Number LP<Number, RowsRemovable, ColumnsRemovable>::getObjectiveValue() const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->getObjectiveValue();
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  bool LP<Number, RowsRemovable, ColumnsRemovable>::hasPrimalSolution() const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->hasPrimalSolution();
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  Number LP<Number, RowsRemovable, ColumnsRemovable>::getPrimalValue(LPColumn column) const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->getPrimalValue(
      column);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  std::vector<Number> LP<Number, RowsRemovable, ColumnsRemovable>::getPrimalSolution(
    const std::vector<LPColumn>& columns) const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->getPrimalSolution(
      columns);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  bool LP<Number, RowsRemovable, ColumnsRemovable>::hasDualSolution() const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->hasDualSolution();
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  Number LP<Number, RowsRemovable, ColumnsRemovable>::getDualValue(LPRow row) const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->getDualValue(row);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  void LP<Number, RowsRemovable, ColumnsRemovable>::setSense(LPSense newSense)
  {
    static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->setSense(newSense);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  LPColumn LP<Number, RowsRemovable, ColumnsRemovable>::addColumn(
    const Number& lowerBound,
    const Number& upperBound,
    const Number& objectiveCoefficient,
    const std::string& name,
    LPColumn unusedColumn)
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->addColumn(
      lowerBound, upperBound, objectiveCoefficient, name, unusedColumn);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  LPRow LP<Number, RowsRemovable, ColumnsRemovable>::addRow(
    const Number& lhs,
    std::size_t numNonzeros,
    const LPColumn* nonzeroVariables,
    const Number* nonzeroCoefficients,
    const Number& rhs,
    const std::string& name,
    LPRow unusedRow)
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->addRow(
      lhs, numNonzeros, nonzeroVariables, nonzeroCoefficients, rhs, name, unusedRow);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  void LP<NumberType, RowsRemovable, ColumnsRemovable>::update()
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->update();
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  void LP<Number, RowsRemovable, ColumnsRemovable>::changeUpper(LPColumn column, const Number& newUpperBound)
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->changeUpper(
      column, newUpperBound);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  void LP<Number, RowsRemovable, ColumnsRemovable>::changeLower(LPColumn column, const Number& newLowerBound)
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->changeLower(
      column, newLowerBound);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  void LP<Number, RowsRemovable, ColumnsRemovable>::changeBounds(LPColumn column, const Number& newLowerBound,
    const Number& newUpperBound)
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->changeBounds(
      column, newLowerBound, newUpperBound);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  void LP<Number, RowsRemovable, ColumnsRemovable>::changeObjective(LPColumn column, const Number& newObjectiveCoefficient)
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->changeObjective(
      column, newObjectiveCoefficient);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  void LP<Number, RowsRemovable, ColumnsRemovable>::write(const std::string& fileName) const
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->write(fileName);
  }

  template <typename Number, bool RowsRemovable, bool ColumnsRemovable>
  LPStatus LP<Number, RowsRemovable, ColumnsRemovable>::solve()
  {
    return static_cast<LPImplementation<Number, RowsRemovable, ColumnsRemovable>*>(_implementation)->solve();
  }

  template class LP<double, true, true>;
  template class LP<double, true, false>;
  template class LP<double, false, true>;
  template class LP<double, false, false>;

#endif /* !IPO_WITH_GUROBI && IPO_WITH_SOPLEX */

#if defined(IPO_WITH_SOPLEX) && defined(IPO_WITH_GMP)

  static mpq_class plusInfinityGMP = soplex::infinity;
  static mpq_class minusInfinityGMP = -soplex::infinity;

  template <>
  const mpq_class& LP<mpq_class>::plusInfinity()
  {
    return plusInfinityGMP;
  };

  template <>
  const mpq_class& LP<mpq_class>::minusInfinity()
  {
    return minusInfinityGMP;
  };

  template class LP<mpq_class, true, true>;
  template class LP<mpq_class, true, false>;
  template class LP<mpq_class, false, true>;
  template class LP<mpq_class, false, false>;

#endif /* IPO_WITH_SOPLEX && IPO_WITH_GMP */
  
} /* namespace lp */
