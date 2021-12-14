#include <ipo/lp.hpp>

#include <ipo/arithmetic.hpp>

#include <cmath>
#include <string>
#include <unordered_map>
#include <ostream>

#if defined(IPO_DOUBLE_LP_SOPLEX) || defined(IPO_RATIONAL_LP_SOPLEX)
#include <soplex.h>
#endif /* IPO_DOUBLE_LP_SOPLEX || IPO_RATIONAL_LP_SOPLEX */

namespace ipo
{

#if defined(IPO_DOUBLE_LP) || defined(IPO_RATIONAL_LP)

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

#endif /* IPO_DOUBLE_LP */

#if defined(IPO_DOUBLE_LP_SOPLEX) || defined(IPO_RATIONAL_LP_SOPLEX)

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

#endif /* IPO_DOUBLE_LP_SOPLEX || IPO_RATIONAL_LP_SOPLEX */

#if defined(IPO_DOUBLE_LP) || defined(IPO_RATIONAL_LP)

  template <typename Number>
  class LPImplementation;

#endif /* IPO_DOUBLE_LP */
  
#if defined(IPO_DOUBLE_LP_SOPLEX)
  
  template <>
  class LPImplementation<double>
  {
  public:
    typedef double Number;

    LPImplementation()
      : _nextRowKey(0), _nextColumnKey(0)
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

    double getPrimalValue(int column) const
    {
      return _primalSolution[column];
    }

    std::vector<double> getPrimalSolution(const std::vector<int>& columns) const
    {
      std::vector<double> solution;
      if (columns.empty())
        return _primalSolution;
      else
      {
        solution.reserve(columns.size());
        for (std::size_t c = 0; c < columns.size(); ++c)
          solution.push_back(_primalSolution[columns[c]]);
      }
      return solution;
    }

    bool hasDualSolution() const
    {
      return !_dualSolution.empty();
    }

    double getDualValue(int row) const
    {
      return _dualSolution[row];
    }

    void setSense(LPSense newSense)
    {
      assert(newSense == LPSense::MAXIMIZE || newSense == LPSense::MINIMIZE);
      _spx.setIntParam(soplex::SoPlex::OBJSENSE,
        newSense == LPSense::MAXIMIZE ? soplex::SoPlex::OBJSENSE_MAXIMIZE : soplex::SoPlex::OBJSENSE_MINIMIZE);
    }

    LPKey addColumn(double lowerBound, double upperBound, double objectiveCoefficient, const std::string& name)
    {
      _sparse.clear();
      _spx.addColReal(soplex::LPColBase<double>(objectiveCoefficient, _sparse,
        std::isnormal(upperBound) ? upperBound : soplex::infinity,
        std::isnormal(lowerBound) ? lowerBound : -soplex::infinity));
      _columnNames.push_back(name);
      _columnKeys.push_back(_nextColumnKey);
      _columnMap[_nextColumnKey] = _spx.numCols() - 1;
      return LPKey(_nextColumnKey++);
    }

    LPKey addRow(double lhs, std::size_t numNonzeros, const int* nonzeroColumns, const double* nonzeroCoefficients,
      double rhs, const std::string& name)
    {
      _sparse.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
        _sparse.add(nonzeroColumns[i], nonzeroCoefficients[i]);
      _spx.addRowReal(soplex::LPRowBase<double>(lhs, _sparse, rhs));
      _rowNames.push_back(name);
      _rowKeys.push_back(_nextRowKey);
      _rowMap[_nextRowKey] = _spx.numRows() - 1;
      return LPKey(_nextRowKey++);
    }

    void update()
    {

    }

    void changeUpper(int column, double newUpperBound)
    {
      _spx.changeUpperReal(column, newUpperBound);
    }

    void changeLower(int column, double newLowerBound)
    {
      _spx.changeLowerReal(column, newLowerBound);
    }

    void changeBounds(int column, double newLowerBound, double newUpperBound)
    {
      _spx.changeBoundsReal(column, newLowerBound, newUpperBound);
    }

    void changeObjective(int column, double newObjectiveCoefficient)
    {
      _spx.changeObjReal(column, newObjectiveCoefficient);
    }

    void changeRow(LPKey rowKey, const double lhs, std::size_t numNonzeros, const int* nonzeroColumns,
      const double* nonzeroCoefficients, const double rhs)
    {
      std::size_t row = _rowMap[rowKey.id];
      std::cout << "changeRow affects row " << row << std::endl;
      for (size_t i = 0; i < numNonzeros; ++i)
        _sparse.add(nonzeroColumns[i], nonzeroCoefficients[i]);
      _spx.changeRowReal(row, soplex::LPRowReal(lhs, _sparse, rhs));
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
    std::vector<LPKey> _rowKeys;
    std::vector<LPKey> _columnKeys;
    int _nextRowKey;
    int _nextColumnKey;
    std::unordered_map<int, std::size_t> _rowMap;
    std::unordered_map<int, std::size_t> _columnMap;
    std::vector<double> _primalSolution;
    std::vector<double> _primalRay;
    std::vector<double> _dualSolution;
    std::vector<double> _dualRay;

    soplex::SoPlex _spx;
    soplex::DSVectorReal _sparse;
    std::vector<std::string> _rowNames;
    std::vector<std::string> _columnNames;
  };

#endif /* IPO_DOUBLE_LP_SOPLEX */

#if defined(IPO_RATIONAL_LP_SOPLEX)
  
  template <>
  class LPImplementation<rational>
  {
  public:
    typedef rational Number;

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

    rational getObjectiveValue() const
    {
      return const_cast<LPImplementation*>(this)->_spx.objValueRational();
    }

    bool hasPrimalSolution() const
    {
      return !_primalSolution.empty();
    }

    const rational& getPrimalValue(int column) const
    {
      return _primalSolution[column];
    }

    std::vector<rational> getPrimalSolution(const std::vector<int>& columns) const
    {
      std::vector<rational> solution;
      if (columns.empty())
        return _primalSolution;
      else
      {
        solution.reserve(columns.size());
        for (std::size_t c = 0; c < columns.size(); ++c)
          solution.push_back(_primalSolution[columns[c]]);
      }
      return solution;
    }

    bool hasDualSolution() const
    {
      return !_dualSolution.empty();
    }

    const rational& getDualValue(int row) const
    {
      return _dualSolution[row];
    }

    void setSense(LPSense newSense)
    {
      assert(newSense == LPSense::MAXIMIZE || newSense == LPSense::MINIMIZE);
      _spx.setIntParam(soplex::SoPlex::OBJSENSE,
        newSense == LPSense::MAXIMIZE ? soplex::SoPlex::OBJSENSE_MAXIMIZE : soplex::SoPlex::OBJSENSE_MINIMIZE);
    }

    LPKey addColumn(const rational& lowerBound, const rational& upperBound,
      const rational& objectiveCoefficient, const std::string& name)
    {
      _sparse.clear();
      _spx.addColRational(soplex::LPColBase<soplex::Rational>(objectiveCoefficient, _sparse, upperBound, lowerBound));
      _columnNames.push_back(name);
      _columnKeys.push_back(_nextColumnKey);
      _columnMap[_nextColumnKey] = _spx.numCols() - 1;
      return _nextColumnKey++;
    }

    LPKey addRow(const rational& lhs, std::size_t numNonzeros, const int* nonzeroColumns,
      const rational* nonzeroCoefficients, const rational& rhs, const std::string& name)
    {
      _sparse.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
        _sparse.add(nonzeroColumns[i], nonzeroCoefficients[i]);
      _spx.addRowRational(soplex::LPRowBase<soplex::Rational>(lhs, _sparse, rhs));
      _rowNames.push_back(name);
      _rowKeys.push_back(_nextRowKey);
      _rowMap[_nextRowKey] = _spx.numRows() - 1;
      return _nextRowKey++;
    }

    void update()
    {

    }

    void changeUpper(int column, const rational& newUpperBound)
    {
      _spx.changeUpperRational(column, newUpperBound);
    }

    void changeLower(int column, const rational& newLowerBound)
    {
      _spx.changeLowerRational(column, newLowerBound);
    }

    void changeBounds(int column, const rational& newLowerBound, const rational& newUpperBound)
    {
      _spx.changeBoundsRational(column, newLowerBound, newUpperBound);
    }

    void changeObjective(int column, const rational& newObjectiveCoefficient)
    {
      _spx.changeObjRational(column, newObjectiveCoefficient);
    }

    void changeRow(LPKey rowKey, const rational& lhs, std::size_t numNonzeros, const int* nonzeroColumns,
      const rational* nonzeroCoefficients, const rational& rhs)
    {
      std::size_t row = _rowMap[rowKey.id];
      _sparse.clear();
      for (size_t i = 0; i < numNonzeros; ++i)
        _sparse.add(nonzeroColumns[i], nonzeroCoefficients[i]);
      _spx.changeRowRational(row, soplex::LPRowRational(lhs, _sparse, rhs));
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
          _primalSolution[c] = _dense[c];
      }
      else
        _primalSolution.clear();
      if (_spx.isDualFeasible())
      {
        _dense.reDim(numRows());
        _spx.getPrimalRational(_dense);
        _dualSolution.resize(numRows());
        for (std::size_t r = 0; r < _dualSolution.size(); ++r)
          _dualSolution[r] = _dense[r];
      }
      else
        _dualSolution.clear();

      return status();
    }

  private:
    std::vector<LPKey> _rowKeys;
    std::vector<LPKey> _columnKeys;
    int _nextRowKey;
    int _nextColumnKey;
    std::unordered_map<int, std::size_t> _rowMap;
    std::unordered_map<int, std::size_t> _columnMap;
    std::vector<rational> _primalSolution;
    std::vector<rational> _primalRay;
    std::vector<rational> _dualSolution;
    std::vector<rational> _dualRay;

    soplex::SoPlex _spx;
    soplex::DVectorRational _dense;
    soplex::DSVectorRational _sparse;
    std::vector<std::string> _rowNames;
    std::vector<std::string> _columnNames;
  };

#endif /* IPO_RATIONAL_LP_SOPLEX */

#if defined(IPO_DOUBLE_LP) || defined(IPO_RATIONAL_LP)

  template <typename Number>
  LP<Number>::LP()
  {
    _implementation = new LPImplementation<Number>();
  }

  template <typename Number>
  LP<Number>::~LP()
  {
    delete static_cast<LPImplementation<Number>*>(_implementation);
  }

  template <typename Number>
  std::size_t LP<Number>::numRows() const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->numRows();
  }

  template <typename Number>
  std::size_t LP<Number>::numColumns() const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->numColumns();
  }

  template <typename Number>
  LPStatus LP<Number>::status() const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->status();
  }

  template <typename Number>
  Number LP<Number>::getObjectiveValue() const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->getObjectiveValue();
  }

  template <typename Number>
  bool LP<Number>::hasPrimalSolution() const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->hasPrimalSolution();
  }

  template <typename Number>
  Number LP<Number>::getPrimalValue(int column) const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->getPrimalValue(column);
  }

  template <typename Number>
  std::vector<Number> LP<Number>::getPrimalSolution(const std::vector<int>& columns) const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->getPrimalSolution(columns);
  }

  template <typename Number>
  bool LP<Number>::hasDualSolution() const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->hasDualSolution();
  }

  template <typename Number>
  Number LP<Number>::getDualValue(int row) const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->getDualValue(row);
  }

  template <typename Number>
  void LP<Number>::setSense(LPSense newSense)
  {
    static_cast<LPImplementation<Number>*>(_implementation)->setSense(newSense);
  }

  template <typename Number>
  LPKey LP<Number>::addColumn(const Number& lowerBound, const Number& upperBound, const Number& objectiveCoefficient,
    const std::string& name)
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->addColumn(lowerBound, upperBound,
      objectiveCoefficient, name);
  }

  template <typename Number>
  LPKey LP<Number>::addRow(const Number& lhs, std::size_t numNonzeros, const int* nonzeroColumns,
    const Number* nonzeroCoefficients, const Number& rhs, const std::string& name)
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->addRow(lhs, numNonzeros, nonzeroColumns,
      nonzeroCoefficients, rhs, name);
  }

  template <typename Number>
  void LP<Number>::update()
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->update();
  }

  template <typename Number>
  void LP<Number>::changeUpper(int column, const Number& newUpperBound)
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->changeUpper(column, newUpperBound);
  }

  template <typename Number>
  void LP<Number>::changeLower(int column, const Number& newLowerBound)
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->changeLower(column, newLowerBound);
  }

  template <typename Number>
  void LP<Number>::changeBounds(int column, const Number& newLowerBound, const Number& newUpperBound)
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->changeBounds(column, newLowerBound, newUpperBound);
  }

  template <typename Number>
  void LP<Number>::changeObjective(int column, const Number& newObjectiveCoefficient)
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->changeObjective(column, newObjectiveCoefficient);
  }

  template <typename Number>
  void LP<Number>::changeRow(LPKey rowKey, const Number& lhs, std::size_t numNonzeros, const int* nonzeroColumns,
    const Number* nonzeroCoefficients, const Number& rhs)
  {
    static_cast<LPImplementation<Number>*>(_implementation)->changeRow(rowKey, lhs, numNonzeros, nonzeroColumns,
      nonzeroCoefficients, rhs);
  }

  template <typename Number>
  void LP<Number>::write(const std::string& fileName) const
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->write(fileName);
  }

  template <typename Number>
  LPStatus LP<Number>::solve()
  {
    return static_cast<LPImplementation<Number>*>(_implementation)->solve();
  }

  template class LP<double>;

#endif /* IPO_DOUBLE_LP || IPO_RATIONAL_LP */

#if defined(IPO_DOUBLE_LP)

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
  
#endif /* IPO_DOUBLE_LP */

#if defined(IPO_RATIONAL_LP_SOPLEX)
  
  static rational plusInfinityRational = soplex::infinity;
  static rational minusInfinityRational = -soplex::infinity;

  template <>
  const rational& LP<rational>::plusInfinity()
  {
    return plusInfinityRational;
  };

  template <>
  const rational& LP<rational>::minusInfinity()
  {
    return minusInfinityRational;
  };

  template class LP<rational>;

#endif /* IPO_RATIONAL_LP_SOPLEX */
  
} /* namespace ipo */
