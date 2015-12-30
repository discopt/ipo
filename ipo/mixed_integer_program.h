#ifndef IPO_MIXED_INTEGER_PROGRAM_H_
#define IPO_MIXED_INTEGER_PROGRAM_H_

#include <string>
#include <vector>

#include <soplex.h>
#include <scip/scip.h>
#include "oracles.h"

namespace ipo {

  class MixedIntegerProgram
  {
  public:
    MixedIntegerProgram(SCIP* scip);

    virtual ~MixedIntegerProgram();

    inline std::size_t numVariables() const
    {
      return _variableNames.size();
    }

    inline std::size_t numConstraints() const
    {
      return _constraintNames.size();
    }

    inline const soplex::LPColSetRational& columns() const
    {
      return _columns;
    }

    inline const soplex::LPRowSetRational& rows() const
    {
      return _rows;
    }

    inline bool isIntegral(std::size_t var) const
    {
      return _integrality[var];
    }

    inline const std::string& variableName(std::size_t var) const
    {
      return _variableNames[var];
    }

    inline const std::string& constraintName(std::size_t cons) const
    {
      return _constraintNames[cons];
    }

    bool checkPointBounds(const soplex::SVectorRational* point) const;
    bool checkPointConstraints(const soplex::SVectorRational* point);
    bool checkPointIntegral(const soplex::SVectorRational* point) const;
    bool checkPoint(const soplex::SVectorRational* point);
    bool checkRayBounds(const soplex::SVectorRational* ray) const;
    bool checkRayConstraints(const soplex::SVectorRational* ray);
    bool checkRay(const soplex::SVectorRational* ray);

    void faceEnabled(Face* face);
    void faceDisabled(Face* face);

    void getConstraints(soplex::LPRowSetRational& rows, bool inequalities, bool equations,
        std::vector<std::string>* names = NULL);
    void getFixedVariableEquations(soplex::LPRowSetRational& rows, std::vector<std::string>* names = NULL);

  protected:
    soplex::LPColSetRational _columns;
    soplex::LPRowSetRational _rows;
    std::vector<std::string> _variableNames;
    std::vector<std::string> _constraintNames;
    std::vector<bool> _integrality;
    soplex::DVectorRational _worker;
    Face* _face;
  };

  class MixedIntegerProgramCorrectorOracle: public FaceOptimizationOracleBase
  {
  public:
    MixedIntegerProgramCorrectorOracle(const std::string& name, MixedIntegerProgram& mip,
        FaceOptimizationOracleBase* inexact, bool correctAlways = true);
    virtual ~MixedIntegerProgramCorrectorOracle();

  protected:

    soplex::DSVectorRational* correctPoint(const soplex::SVectorRational* point,
        const soplex::VectorRational& objective);
    soplex::DSVectorRational* correctRay(const soplex::SVectorRational* ray, const soplex::VectorRational& objective);

    virtual void run(OptimizationResult& result, const soplex::VectorRational& objective,
        const soplex::Rational* improveValue, bool forceOptimal);

    virtual void faceEnabled(Face* face);
    virtual void faceDisabled(Face* face);

  protected:

    MixedIntegerProgram& _mip;
    FaceOptimizationOracleBase* _inexact;
    soplex::SoPlex _spx;
    soplex::DVectorRational _worker;
    bool _correctAlways;
  };

} /* namespace polycomb */

#endif /* IPO_MIXED_INTEGER_PROGRAM_H_ */
