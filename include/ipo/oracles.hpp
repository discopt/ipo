#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>
#include <ipo/space.hpp>

#include <memory>

namespace ipo
{
   /**
    * \brief Base class for all IPO oracles.
    */

   class Oracle
   {
   public:
      /**
       * \brief Constructs an oracle with given \p name.
       */

      IPO_EXPORT
      Oracle(const std::string& name);

      /**
       * \brief Returns the oracle's name.
       */

      IPO_EXPORT
      inline const std::string& name() const
      {
         return _name;
      }

      /**
       * \brief Returns the oracle's ambient space.
       */

      IPO_EXPORT
      inline std::shared_ptr<Space> space() const
      {
         return _space;
      }

   protected:
      std::string _name;
      std::shared_ptr<Space> _space;
   };

   

   class OptimizationOracle: public Oracle
   {
   public:
      /**
       * \brief Structure for storing the query.
       **/

      struct Query
      {
#ifdef IPO_WITH_GMP
         /// Are rational values requested? If yes, floating-point and rational values are required.
         bool rational;
#endif /* IPO_WITH_GMP */
         /// The caller will only use points having at least this value.
         double minObjectiveValue;
         /// Minimum number of solutions for early termination.
         int minNumSolutions;
         /// Maximum number of solutions to return.
         int maxNumSolutions;
         /// Time limit
         double timeLimit;

         /**
          * \brief Constructs the query structure.
          */

         IPO_EXPORT
         Query();

         /**
          * \brief Clears the query data.
          */

         IPO_EXPORT
         void reset();
      };
      
      /**
       * \brief Structure for storing the query result.
       **/

      struct Result
      {
         /**
          * \brief Array of objective values of all points (\c std::numeric_limits<double>::nan 
          * otherwise).
          */
         std::vector<double> objectiveValues;
         /// Array that maps a point/ray to the index of its first entry.
         std::vector<std::size_t> firstIndices;
         /// Array that maps indices to coordinates.
         std::vector<std::size_t> nonzeroCoordinates;
         /// Array that maps indices to floating-point values.
         std::vector<double> nonzeroValues;
         /**
          * \brief Upper bound on the optimal solution value (may be \c
          * std::numeric_limits<double>::infinity).
          */
         double dualBound;

#ifdef IPO_WITH_GMP
         /// Array that maps indices to rational values. 
         std::vector<mpq_class> rationalNonzeroValues;
         ///Rational value of \ref upperBound (if finite).
         double rationalDualBound;
#endif /* IPO_WITH_GMP */

         /**
          * \brief Constructs the result structure.
          */

         IPO_EXPORT
         Result();

         /**
          * \brief Clears the result data.
          */

         IPO_EXPORT
         void reset();
      };

      /**
       * \brief Constructs the oracle.
       */

      IPO_EXPORT
      OptimizationOracle(const std::string& name);

      /**
       * \brief Returns true iff the oracle is exact.
       *
       * Returns true iff the oracle is exact, i.e., upon request it can return solutions as exact
       * rational vectors.
       */

      IPO_EXPORT
      bool isExact() const;

      /**
       * \brief Maximize a floating-point objective vector.
       * 
       * \param objectiveVector Array that maps coordinates to objective value coefficients.
       * \param query Structure for query.
       * \param result Structure for returning the result.
       **/

      virtual void maximize(const double* objectiveVector, const Query& query, Result& result) = 0;

#ifdef IPO_WITH_GMP

      /**
       * \brief Maximize a rational objective vector.
       *
       * Maximize a rational objective vector. The default implementation just converts the
       * objective vector to a floating-point vector and calls \ref maximize.
       *
       * \param objectiveVector Array that maps coordinates to objective value coefficients.
       * \param query Structure for query.
       * \param result Structure for returning the result.
       **/

      IPO_EXPORT
      virtual void maximize(const mpq_class* objectiveVector, const Query& query, Result& result);

#endif /* IPO_WITH_GMP */

   };

   /**
    * \brief Base class for a separation oracle.
    * 
    * Base class for a separation oracle. When queried with a point, the oracle returns any number
    * (including none) of (less-than-or-equal) inequalities.
    */

   class SeparationOracle: public Oracle
   {
   public:
      /**
       * \brief Structure for storing the query.
       **/

      struct Query
      {
#ifdef IPO_WITH_GMP
         /// Are rational values requested? If yes, floating-point and rational values are required.
         bool rational;
#endif /* IPO_WITH_GMP */
         /// Maximum number of solutions to return.
         int maxNumInequalities;
         /// Time limit
         double timeLimit;

         /**
          * \brief Constructs the query structure.
          */

         IPO_EXPORT
         Query();

         /**
          * \brief Clears the query data.
          */

         IPO_EXPORT
         void reset();
      };

      /**
       * \brief Structure for storing the query result.
       **/

      struct Result
      {
         /// Array that maps an inequality index to its right-hand side.
         std::vector<double> rightHandSides;
         /// Array that maps an inequality to the index of its first entry.
         std::vector<std::size_t> firstIndices;
         /// Array that maps indices to coordinates.
         std::vector<std::size_t> nonzeroCoordinates;
         /// Array that maps indices to floating-point coefficients.
         std::vector<double> nonzeroCoefficients;

#ifdef   IPO_WITH_GMP
         /// Array that maps an inequality index to its rational right-hand side.
         std::vector<mpq_class> rationalRightHandSides;
         /// Array that maps indices to rational coefficients. 
         std::vector<mpq_class> rationalNonzeroCoefficients;
#endif /* IPO_WITH_GMP */

         /**
          * \brief Constructs the result structure.
          */

         IPO_EXPORT
         Result();

         /**
          * \brief Clears the result data.
          */

         IPO_EXPORT
         void reset();
      };

      
      /**
       * \brief Constructs the oracle.
       */

      IPO_EXPORT
      SeparationOracle(const std::string& name);

      /**
       * \brief Separates a point/ray with floating-point coordinates.
       * 
       * \param vector Array that maps coordinates to point/ray coordinates.
       * \param query Structure for query.
       * \param isPoint Whether a point shall be separated.
       * \param result Structure for returning the result.
       * 
       * \returns \c true if and only if the point/ray was separated.
       **/

      virtual bool separate(const double* vector, bool isPoint, const Query& query,
         Result& result) = 0;

#ifdef IPO_WITH_GMP

      /**
       * \brief Separates a point/ray with rational coordinates.
       * 
       * \param vector Array that maps coordinates to point/ray coordinates.
       * \param query Structure for query.
       * \param isPoint Whether a point shall be separated.
       * \param result Structure for returning the result.
       *
       * \returns \c true if and only if the point/ray was separated.
       **/

      IPO_EXPORT
      virtual void separate(const mpq_class* vector, bool isPoint, const Query& query,
         Result& result);

#endif /* IPO_WITH_GMP */
   };

} /* namespace ipo */
