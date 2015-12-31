#ifndef IPO_ORACLES_H_
#define IPO_ORACLES_H_

#include <vector>
#include <limits>

#include "spx_gmp.h"
#include "statistics.h"
#include "rows.h"

namespace ipo {

  /**
   * \brief Defines a face of a polyhedron by a set of inequalities.
   * 
   * Defines a face \f$F\f$ of a polyhedron \f$P\f$ by a set of inequalities.
   * It is used to create an optimization oracle for \f$F\f$.
   * If the optimization oracle class for \f$P\f$ inherits from
   * \ref FaceOptimizationOracleBase, the face can be controlled directly.
   * In any case one can construct an instance of \ref \FaceOptimizationOracle
   * that calls the optimization oracle for \f$P\f$ (maybe multiple times per call).
   **/

  class Face
  {
  public:
    /**
     * Creates a face
     **/

    Face(std::size_t numVariables);
    Face(std::size_t numVariables, const soplex::LPRowRational& inequality);
    Face(std::size_t numVariables, const soplex::LPRowSetRational& inequalities);
    virtual ~Face();

    void add(const soplex::LPRowRational& inequality);
    void add(const soplex::LPRowSetRational& inequalities);

    inline const soplex::LPRowSetRational& inequalities() const
    {
      return _inequalities;
    }

    inline const soplex::SVectorRational& normal()
    {
      ensureSync();
      return _normal;
    }

    inline const soplex::Rational& rhs()
    {
      ensureSync();
      return _rhs;
    }

    inline const soplex::Rational& largestAbsCoefficient()
    {
      ensureSync();
      return _largestAbsCoefficient;
    }

  protected:
    void ensureSync();

    soplex::LPRowSetRational _inequalities;
    bool _synced;
    soplex::DVectorRational _worker;
    soplex::DSVectorRational _normal;
    soplex::Rational _largestAbsCoefficient;
    soplex::Rational _rhs;
  };

  struct OptimizationResult
  {
    bool optimal;
    std::size_t bestIndex;
    soplex::Rational bestValue;
    std::vector<soplex::DSVectorRational*> points;
    std::vector<soplex::Rational> objectives;
    std::vector<soplex::DSVectorRational*> rays;

    inline bool isInfeasible() const
    {
      return objectives.empty() && rays.empty();
    }

    inline bool isUnbounded() const
    {
      return !rays.empty();
    }

    inline bool isFeasible() const
    {
      return !objectives.empty();
    }

    void reset(std::size_t numVariables);
    soplex::DSVectorRational& newPoint();
    void setFeasible(const soplex::VectorRational& objective);
    void setBest(const soplex::Rational& value, std::size_t index = std::numeric_limits<std::size_t>::max());
    void setInfeasible();
    void setUnbounded();

    void filterDuplicates();

#ifndef NDEBUG
  public:
    void checkConsistent() const;
    bool hasDuplicates() const;
#endif

  protected:
    void filterDuplicates(std::vector<soplex::DSVectorRational*>& vectors);

  protected:
    std::size_t _numVariables;

  };

  class OptimizationOracleBase
  {
  public:
    virtual ~OptimizationOracleBase();

    inline std::size_t numVariables() const
    {
      assert(_initialized);
      return _variableNames.size();
    }

    inline const std::string& variableName(std::size_t var) const
    {
      assert(_initialized);
      return _variableNames[var];
    }

    inline const std::string& name() const
    {
      return _name;
    }

    void printRow(std::ostream& stream, const soplex::LPRowRational& row) const;
    void printRow(std::ostream& stream, const soplex::LPRowSetRational& rows, std::size_t index) const;
    void printRows(std::ostream& stream, const soplex::LPRowSetRational& rows) const;
    void printVector(std::ostream& stream, const soplex::SVectorRational* vector) const;

    void maximize(OptimizationResult& result, const soplex::VectorRational& objective, bool forceOptimal = true);
    void maximize(OptimizationResult& result, const soplex::VectorReal& objective, bool forceOptimal = true);
    void maximize(OptimizationResult& result, const soplex::SVectorRational& objective, bool forceOptimal = true);
    void maximize(OptimizationResult& result, const soplex::SVectorReal& objective, bool forceOptimal = true);
    void improve(OptimizationResult& result, const soplex::VectorRational& objective, const soplex::Rational& value,
        bool forceOptimal = true);
    void improve(OptimizationResult& result, const soplex::VectorReal& objective, const soplex::Rational& value,
        bool forceOptimal = true);
    void improve(OptimizationResult& result, const soplex::SVectorRational& objective, const soplex::Rational& value,
        bool forceOptimal = true);
    void improve(OptimizationResult& result, const soplex::SVectorReal& objective, const soplex::Rational& value,
        bool forceOptimal = true);

  protected:
    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    OptimizationOracleBase(const std::string& name);

    void printRow(std::ostream& stream, const soplex::Rational* lhs, const soplex::Rational* rhs,
        const soplex::SVectorRational& vector) const;

    void initialize(const std::vector<std::string>& variableNames);
    void initialize(const OptimizationOracleBase* oracle);

  private:
    OptimizationOracleBase();

    std::string _name;
    std::vector<std::string> _variableNames;
    soplex::DVectorRational _objective;
#ifndef NDEBUG
    bool _initialized;
#endif
  };

  class ChainedOptimizationOracle: public OptimizationOracleBase
  {
  public:
    ChainedOptimizationOracle(OptimizationOracleBase* first, OptimizationOracleBase* second);
    virtual ~ChainedOptimizationOracle();

  protected:
    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

  protected:
    OptimizationOracleBase* _first;
    OptimizationOracleBase* _second;
  };

  class FaceOptimizationOracleBase: public OptimizationOracleBase
  {
  public:
    virtual ~FaceOptimizationOracleBase();

    virtual Face* setFace(Face* face = NULL);

  protected:
    FaceOptimizationOracleBase(const std::string& name);

    virtual void faceEnabled(Face* face) = 0;
    virtual void faceDisabled(Face* face) = 0;

  protected:
    Face* _face;
  };

  class Projection
  {
  public:
    Projection();
    Projection(const OptimizationOracleBase* oracle, const std::vector<std::size_t>& variableSubset);
    virtual ~Projection();

    void addVariable(const std::string& name, const soplex::SVectorRational& variableMap,
        const soplex::Rational& shift = soplex::Rational(0));
    void addVariable(const OptimizationOracleBase* oracle, std::size_t originalVariable, const soplex::Rational& shift =
        soplex::Rational(0));
    void projectPoint(const soplex::DVectorRational& point, soplex::DSVectorRational& image) const;

    inline std::size_t numVariables() const
    {
      return _map.size();
    }

    inline const std::vector<std::string>& names() const
    {
      return _names;
    }

    inline const soplex::SVectorRational& map(std::size_t var) const
    {
      return _map[var];
    }

    inline const soplex::Rational& shift(std::size_t var) const
    {
      return _shift[var];
    }

  protected:
    std::vector<std::string> _names;
    std::vector<soplex::DSVectorRational> _map;
    std::vector<soplex::Rational> _shift;
  };

  class ProjectedOptimizationOracle: public OptimizationOracleBase
  {
  public:
    ProjectedOptimizationOracle(const std::string& name, const Projection& projection, OptimizationOracleBase* oracle);
    virtual ~ProjectedOptimizationOracle();

  protected:
    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

  protected:
    const Projection& _projection;
    OptimizationOracleBase* _oracle;
    soplex::DVectorRational _liftedObjective;
    OptimizationResult _result;
  };

} /* namespace polycomb */

#endif /* IPO_ORACLE_H_ */
