#ifndef IPO_SCIP_ORACLES_H_
#define IPO_SCIP_ORACLES_H_

#include "oracles.h"
#include "mixed_integer_program.h"

#include <set>
#include <map>
#include <limits>

#include <scip/scip.h>

namespace ipo {

  typedef std::map<SCIP_VAR*, std::size_t> SCIPvarToIndexMap;

  void getSCIPvarToIndexMap(SCIP* originalSCIP, SCIPvarToIndexMap& map);

  void getSCIPObjective(SCIP* originalSCIP, soplex::DVectorRational& objective,
      bool makeMaximization = true);

  class SCIPOptimizationOracle: public FaceOptimizationOracleBase
  {
  public:
    SCIPOptimizationOracle(const std::string& name, SCIP* originalSCIP, bool isHeuristic = true);
    SCIPOptimizationOracle(const std::string& name, const MixedIntegerProgram& mip, bool isHeuristic = true);
    virtual ~SCIPOptimizationOracle();

  protected:
    virtual void run(OptimizationResult& result,
        const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    virtual void faceEnabled(Face* face);
    virtual void faceDisabled(Face* face);

    soplex::Rational computeVectorScaling(const soplex::VectorRational& vector);

  protected:
    SCIP* _scip;
    bool _isHeuristic;
    soplex::Rational _maxLargestCoefficient;
    soplex::Rational _minLargestCoefficient;
    soplex::Rational _bestLargestCoefficient;
    std::vector<SCIP_VAR*> _variables;
    SCIP_CONS* _faceConstraint;
  };

  class ExactSCIPOptimizationOracle: public FaceOptimizationOracleBase
  {
  public:
    ExactSCIPOptimizationOracle(const std::string& name,
        const std::string& exactBinary, MixedIntegerProgram& mip,
        FaceOptimizationOracleBase* heuristic, double timeLimit = std::numeric_limits<double>::max());
    virtual ~ExactSCIPOptimizationOracle();

    bool scaleObjective;

  protected:

    virtual void run(OptimizationResult& result,
        const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    virtual void faceEnabled(Face* face);
    virtual void faceDisabled(Face* face);

    void createTempDirectory();
    void deleteTempDirectory();
    void writeModel(const soplex::VectorRational& objective);
    void callSolver();
    soplex::DSVectorRational* parseOutput();

  protected:

    const std::string _binary;
    MixedIntegerProgram& _mip;
    FaceOptimizationOracleBase* _heuristic;
    std::string _path;
    double _timeLimit;
  };

} /* namespace polycomb */

#endif /* IPO_SCIP_ORACLES_H_ */
