#pragma once

#include <ipo/config.hpp>
#include <ipo/export.hpp>

#ifdef IPO_WITH_SCIP

#include <ipo/oracles.hpp>
#include <map>

// This is necessary due to a bug in SCIP. Whether some functionality is in a macro or not depends 
// on how you include it.
#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
#endif

namespace ipo
{
   class SCIPOptimizationOracle;
   class SCIPSeparationOracle;
   
   class SCIPSolver
   {
   public:
      typedef std::size_t Face;

      IPO_EXPORT
      SCIPSolver(SCIP* scip);

      IPO_EXPORT
      SCIPSolver(const std::string& fileName);

      IPO_EXPORT
      ~SCIPSolver();

      inline const std::string& name() const
      {
         return _name;
      }

      inline std::shared_ptr<Space> space() const
      {
         return _space;
      }
   
      inline const double* instanceObjective() const
      {
         return _instanceObjective;
      }

      IPO_EXPORT
      Face addFace(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
         double* nonzeroCoefficients, double rhs);

#ifdef IPO_WITH_GMP
      IPO_EXPORT
      Face addFace(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
         mpq_class* nonzeroCoefficients, const mpq_class& rhs);
#endif /* IPO_WITH_GMP */

      inline
      std::shared_ptr<SCIPOptimizationOracle> getOptimizationOracle(Face face = 0)
      {
         std::shared_ptr<SCIPSolver> solver(this);
         return std::make_shared<SCIPOptimizationOracle>(solver, face);
      }

      inline
      std::shared_ptr<SCIPSeparationOracle> getSeparationOracle(Face face = 0)
      {
         std::shared_ptr<SCIPSolver> solver(this);
         return std::make_shared<SCIPSeparationOracle>(solver, face);
      }

   protected:
      
      friend SCIPOptimizationOracle;
      friend SCIPSeparationOracle;
      
      void initialize();

      IPO_EXPORT
      void setFace(SCIPSolver::Face face);

   protected:
      SCIP* _scip;
      std::vector<SCIP_VAR*> _variables;
      std::map<SCIP_VAR*, std::size_t> _variablesToCoordinates;
      double* _instanceObjective;
      std::string _name;
      std::shared_ptr<Space> _space;
      Face _currentFace;
      std::vector<SCIP_CONS*> _faceConstraints;
   };

   /**
    * \brief OptimizationOracle based on the SCIP solver.
    */

   class SCIPOptimizationOracle: public OptimizationOracle
   {
   public:

      /**
       * \brief Constructs oracle using the \p solver.
       * 
       * \param solver The solver instance that is used to answer the queries.
       * \param faceIndex The face we are optimizing over.
       */

      IPO_EXPORT
      SCIPOptimizationOracle(std::shared_ptr<SCIPSolver> solver, SCIPSolver::Face face);

      /**
       * \brief Destructor.
       */

      IPO_EXPORT
      virtual ~SCIPOptimizationOracle();

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

      IPO_EXPORT
      virtual void maximize(const double* objectiveVector, const Query& query, Result& result);

   protected:
      friend SCIPSolver;

      /// The solver instance
      std::shared_ptr<SCIPSolver> _solver;
      /// The index of the face we are optimizing over.
      SCIPSolver::Face _face;
   };

   /**
    * \brief SeparationOracle for the LP relaxation based on the SCIP solver.
    */
   
   class SCIPSeparationOracle: public SeparationOracle
   {
   public:
      /**
       * \brief Constructs oracle using the \p solver.
       * 
       * \param solver The solver instance that is used to answer the queries.
       * \param faceIndex Indexes the face we are separating for.
       */

      IPO_EXPORT
      SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver, SCIPSolver::Face face);

      /**
       * \brief Destructor.
       */

      IPO_EXPORT
      virtual ~SCIPSeparationOracle();

      /**
       * \brief Separates a point/ray with floating-point coordinates.
       * 
       * \param vector Array that maps coordinates to point/ray coordinates.
       * \param result Structure for returning the result.
       * \param timeLimit Time limit for this call (in seconds).
       * 
       * \returns \c true if and only if the point/ray was separated.
       */

      IPO_EXPORT
      virtual bool separate(const double* vector, bool isPoint, 
         const SeparationOracle::Query& query, SeparationOracle::Result& result);

   protected:
      friend SCIPSolver;

      /// The solver instance
      std::shared_ptr<SCIPSolver> _solver;
      /// The index of the face we are separating for.
      SCIPSolver::Face _face;
   };


//   class MIPOracleBase;
// 
//   /**
//    * A map to associate SCIP variables with variable indices.
//    */
// 
//   typedef std::map<SCIP_VAR*, std::size_t> SCIPvarToIndexMap;
// 
//   /**
//    * \brief Returns a canonical mapping from SCIP variables to variable indices.
//    *
//    * Returns a canonical mapping from SCIP variables to variable indices.
//    * It maps original variables in the order SCIP returns them.
//    * If the SCIP instance was transformed already, then
//    * it also maps the corresponding transformed variables to the same indices.
//    *
//    * \param originalSCIP
//    *   SCIP instance
//    * \param map
//    *   Variable map the method writes into.
//    */
// 
//   void getSCIPvarToIndexMap(SCIP* originalSCIP, SCIPvarToIndexMap& map);
// 
//   /**
//    * \brief Extracts the objective from a SCIP instance.
//    *
//    * Extracts the objective from the given SCIP instance. The variables are ordered canonically (see \ref getSCIPvarToIndexMap()).
//    * By default, the objective is scaled such that it corresponds to a maximization problem.
//    *
//    * \param scip
//    *   SCIP instance
//    * \param makeMaximization
//    *   If \c true, the objective is scaled such that it corresponds to a maximization problem.
//    */
// 
//   Vector getSCIPObjective(SCIP* scip, bool makeMaximization = true);
// 
//   /**
//    * \brief Extracts the objective from a SCIP instance that is read from \p fileName.
//    *
//    * Extracts the objective from the given SCIP instance that is read from \p fileName. The variables are ordered canonically (see 
//    * \ref getSCIPvarToIndexMap()). By default, the objective is scaled such that it corresponds to a maximization problem.
//    *
//    * \param fileName
//    *   File name for SCIP instance.
//    * \param makeMaximization
//    *   If \c true, the objective is scaled such that it corresponds to a maximization problem.
//    */
// 
//   Vector getSCIPObjective(const std::string& fileName, bool makeMaximization = true);
// 
//   /**
//    * \brief An oracle based on the SCIP solver.
//    *
//    * An oracle for the convex hull of the solutions returned by the SCIP instance. The computed floating-point solutions are
//    * turned into rational ones by the underlying \ref MIPOracleBase.
//    *
//    * An instance is either constructed from a \c SCIP instance or from a \ref MixedIntegerLinearSet.
//    near*/
// 
//   class SCIPOracle: public MIPOracleBase
//   {
//   public:
//     /**
//      * \brief Constructs a SCIP oracle with given \p name, optionally associated to \p nextOracle.
//      *
//      * Constructs a SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is
//      * defined via the \p originalSCIP instance (and must be equal to that of \p nextOracle). The oracle is implemented by
//      * calling SCIP on a copy of the given \p originalSCIP instance.
//      */
// 
//     SCIPOracle(const std::string& name, SCIP* originalSCIP, const std::shared_ptr<OracleBase>& nextOracle = NULL);
// 
//     /**
//      * \brief Constructs a SCIP oracle with given \p name in given \p space.
//      *
//      * Constructs a SCIP oracle with given \p name that is optionally associated to \p nextOracle. The ambient space is equal to
//      * that of \p nextOracle and of the space of the given \p mixedIntegerSet. The oracle is implemented by calling SCIP in order
//      * to solve mixed-integer programs over the \p mixedIntegerSet.
//      */
// 
//     SCIPOracle(const std::string& name, const std::shared_ptr<MixedIntegerLinearSet>& mixedIntegerLinearSet,
//       const std::shared_ptr<OracleBase>& nextOracle = NULL);
// 
//     /**
//      * \brief Constructs a SCIP oracle from a file specified by \p fileName, optionally associated to \p nextOracle.
//      *
//      * Constructs a SCIP oracle from a file specified by \p fileName, optionally associated to \p nextOracle. The ambient space is
//      * defined by the model read from the file (and must be equal to that of \p nextOracle). The oracle is implemented by calling 
//      * SCIP in order to solve mixed-integer programs.
//      */
// 
//     SCIPOracle(const std::string& fileName, const std::shared_ptr<OracleBase>& nextOracle = NULL);
// 
//     /**
//      * \brief Destructor.
//      */
// 
//     virtual ~SCIPOracle();
// 
//     /**
//      * \brief Restricts the oracle to the face defined by \p newFace.
//      *
//      * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
//      * For \p newFace equal to \c NULL we define \f$ F := P \f$.
//      *
//      * This implementation adds an equation to the underlying SCIP instance.
//      */
// 
//     virtual void setFace(const LinearConstraint& newFace = completeFaceConstraint());
// 
//     /**
//      * \brief Sets a time limit for each oracle call.
//      *
//      * Sets a time limit (in seconds) for each oracle call. If \c heuristicLevel is 0, this raises an exception.
//      */
// 
//     double setTimeLimit(double timeLimit);
// 
//     /**
//      * \brief Returns the time limit set for each oracle call.
//      *
//      * Returns the time limit set via setTimeLimit().
//      */
// 
//     double getTimeLimit();
// 
//   protected:
// 
//     std::shared_ptr<MixedIntegerLinearSet> constructFromSCIP(SCIP* originalSCIP);
// 
//     void constructFromMixedIntegerLinearSet(const std::shared_ptr<MixedIntegerLinearSet>& mixedIntegerLinearSet);
// 
//     std::shared_ptr<MixedIntegerLinearSet> constructFromFile(const std::string& fileName);
// 
//     virtual void solverMaximize(double* objective, double objectiveBound, std::vector<double*>& points,
//       std::vector<double*>& rays, bool& hitLimit);
// 
//   protected:
//     SCIP* _scip; // SCIP instance
//     std::vector<SCIP_VAR*> _variables; // SCIP variables.
//     SCIP_CONS* _faceConstraint; // Special equation constraint for optimizing over a face.
//   };

} /* namespace ipo */

#endif /* IPO_WITH_SCIP */
