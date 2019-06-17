#include <ipo/oracles.hpp>

#include <algorithm>
#include <cassert>

namespace ipo
{

   Oracle::Oracle(const std::string& name)
      : _name(name), _space(nullptr)
   {

   }

   OptimizationOracle::Query::Query()
   {
      reset();
   }
   
   void OptimizationOracle::Query::reset()
   {
#ifdef IPO_WITH_GMP
      rational = false;
#endif
      minNumSolutions = 1;
      maxNumSolutions = 10;
      minObjectiveValue = std::numeric_limits<double>::infinity();
      timeLimit = std::numeric_limits<double>::infinity();
   }

   OptimizationOracle::Result::Result()
   {
      reset();
   }

   void OptimizationOracle::Result::reset()
   {
      objectiveValues.clear();
      firstIndices.clear();
      nonzeroCoordinates.clear();
      nonzeroValues.clear();
#ifdef IPO_WITH_GMP
      rationalNonzeroValues.clear();
#endif /* IPO_WITH_GMP */
   }

   OptimizationOracle::OptimizationOracle(const std::string& name)
      : Oracle(name)
   {

   }

#ifdef IPO_WITH_GMP

   void OptimizationOracle::maximize(const mpq_class* objectiveVector,
      const OptimizationOracle::Query& query, OptimizationOracle::Result& result)
   {
      // Create floating-point approximation of objective vector.
      double* approximateObjectiveVector = new double[space()->dimension()];
      for (std::size_t i = 0; i < space()->dimension(); ++i)
         approximateObjectiveVector[i] = objectiveVector[i].get_d();

      this->maximize(approximateObjectiveVector, query, result);

      delete[] approximateObjectiveVector;
   }

#endif /* IPO_WITH_GMP */

   SeparationOracle::Query::Query()
   {
      reset();
   }

   void SeparationOracle::Query::reset()
   {
#ifdef IPO_WITH_GMP
      rational = false;
#endif /* IPO_WITH_GMP */
      maxNumInequalities = 50;
   }

   SeparationOracle::Result::Result()
   {
      reset();
   }

   void SeparationOracle::Result::reset()
   {
      rightHandSides.clear();
      firstIndices.clear();
      nonzeroCoefficients.clear();
      nonzeroCoordinates.clear();
#ifdef IPO_WITH_GMP
      rationalNonzeroCoefficients.clear();
      rationalRightHandSides.clear();
#endif /* IPO_WITH_GMP */
   }

   SeparationOracle::SeparationOracle(const std::string& name)
      : Oracle(name)
   {

   }

   void SeparationOracle::separate(const mpq_class* vector, bool isPoint,
      const SeparationOracle::Query& query, SeparationOracle::Result& result)
   {
      // Create floating-point approximation of point.
      double* approximateVector = new double[space()->dimension()];
      for (std::size_t i = 0; i < space()->dimension(); ++i)
         approximateVector[i] = vector[i].get_d();

      this->separate(approximateVector, isPoint, query, result);

      delete[] approximateVector;
   }

} /* namespace ipo */
