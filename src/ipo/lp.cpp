#include <ipo/lp.hpp>

#include <cmath>
#include <string>

#if defined(IPO_WITH_SOPLEX)
#include <soplex.h>
#endif /* IPO_WITH_SOPLEX */


namespace ipo
{

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
      default:
        std::ostringstream stream;
        stream << "Unknown SoPlex status " << _spx.status();
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
      _vector.clear();
      _spx.addColReal(soplex::LPColBase<double>(objectiveCoefficient, _vector,
        std::isnormal(upperBound) ? upperBound : soplex::infinity,
        std::isnormal(lowerBound) ? lowerBound : -soplex::infinity));
      _columnNames.push_back(name);
      return _columnMapping.add(numColumns(), unusedColumn);
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
      _vector.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
      {
        std::size_t j = _columnMapping.externalToInternal(nonzeroVariables[i]);
        assert(j < numColumns());
        _vector.add(j, nonzeroCoefficients[i]);
      }
      _spx.addRowReal(soplex::LPRowBase<double>(lhs, _vector, rhs));
      _rowNames.push_back(name);
      return _rowMapping.add(numRows(), unusedRow);
    }

    void update()
    {

    }

    void changeObjective(LPColumn column, double newObjectiveCoefficient)
    {
      _spx.changeObjReal(_columnMapping.externalToInternal(column), newObjectiveCoefficient);
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

    LPImplementation()
    {

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
      default:
        std::ostringstream stream;
        stream << "Unknown SoPlex status " << _spx.status();
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

    double getPrimalValue(LPColumn column) const
    {
      size_t j = _columnMapping.externalToInternal(column);
      assert(j < _primalSolution.size());
      return _primalSolution[j];
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
      _vector.clear();
      _spx.addColReal(soplex::LPColBase<double>(objectiveCoefficient, _vector,
        std::isnormal(upperBound) ? upperBound : soplex::infinity,
        std::isnormal(lowerBound) ? lowerBound : -soplex::infinity));
      _columnNames.push_back(name);
      return _columnMapping.add(numColumns(), unusedColumn);
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
      _vector.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
      {
        std::size_t j = _columnMapping.externalToInternal(nonzeroVariables[i]);
        assert(j < numColumns());
        _vector.add(j, nonzeroCoefficients[i]);
      }
      _spx.addRowReal(soplex::LPRowBase<double>(lhs, _vector, rhs));
      _rowNames.push_back(name);
      return _rowMapping.add(numRows(), unusedRow);
    }

    void update()
    {

    }

    void changeObjective(LPColumn column, double newObjectiveCoefficient)
    {
      _spx.changeObjReal(_columnMapping.externalToInternal(column), newObjectiveCoefficient);
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

#endif /* IPO_WITH_SOPLEX && IPO_WITH_GMP */

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

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  double LP<NumberType, RowsRemovable, ColumnsRemovable>::getObjectiveValue() const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->getObjectiveValue();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  bool LP<NumberType, RowsRemovable, ColumnsRemovable>::hasPrimalSolution() const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->hasPrimalSolution();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  double LP<NumberType, RowsRemovable, ColumnsRemovable>::getPrimalValue(LPColumn column) const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->getPrimalValue(
      column);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  bool LP<NumberType, RowsRemovable, ColumnsRemovable>::hasDualSolution() const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->hasDualSolution();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  double LP<NumberType, RowsRemovable, ColumnsRemovable>::getDualValue(LPRow row) const
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->getDualValue(row);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  void LP<NumberType, RowsRemovable, ColumnsRemovable>::setSense(LPSense newSense)
  {
    static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->setSense(newSense);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  LPColumn LP<NumberType, RowsRemovable, ColumnsRemovable>::addColumn(
    double lowerBound,
    double upperBound,
    double objectiveCoefficient,
    const std::string& name,
    LPColumn unusedColumn)
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->addColumn(
      lowerBound, upperBound, objectiveCoefficient, name, unusedColumn);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  LPRow LP<NumberType, RowsRemovable, ColumnsRemovable>::addRow(
    double lhs,
    std::size_t numNonzeros,
    const LPColumn* nonzeroVariables,
    const double* nonzeroCoefficients,
    double rhs,
    const std::string& name,
    LPRow unusedRow)
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->addRow(
      lhs, numNonzeros, nonzeroVariables, nonzeroCoefficients, rhs, name, unusedRow);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  void LP<NumberType, RowsRemovable, ColumnsRemovable>::update()
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->update();
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  void LP<NumberType, RowsRemovable, ColumnsRemovable>::changeObjective(LPColumn column, double newObjectiveCoefficient)
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->changeObjective(
      column, newObjectiveCoefficient);
  }

  template <typename NumberType, bool RowsRemovable, bool ColumnsRemovable>
  LPStatus LP<NumberType, RowsRemovable, ColumnsRemovable>::solve()
  {
    return static_cast<LPImplementation<NumberType, RowsRemovable, ColumnsRemovable>*>(_implementation)->solve();
  }

  template class LP<double, true, true>;
  template class LP<double, true, false>;
  template class LP<double, false, true>;
  template class LP<double, false, false>;

#endif /* !IPO_WITH_GUROBI && IPO_WITH_SOPLEX */

#if defined(IPO_WITH_SOPLEX) && defined(IPO_WITH_GMP)

  template class LP<mpq_class, true, true>;
  template class LP<mpq_class, true, false>;
  template class LP<mpq_class, false, true>;
  template class LP<mpq_class, false, false>;

#endif /* IPO_WITH_SOPLEX && IPO_WITH_GMP */
  
} /* namespace lp */
