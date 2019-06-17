#include <ipo/oracles_scip.hpp>

#include <cassert>

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <scip/pub_misc.h>
  #define NDEBUG
#else
  #include <scip/scipdefplugins.h>
  #include <scip/cons_linear.h>
  #include <scip/pub_misc.h>
#endif


/**
 * \brief Macro to raise a SCIPException in case of a SCIP error.
 **/

#define SCIP_CALL_EXC(x) \
{ \
   SCIP_RETCODE _retcode; \
   if ((_retcode = (x)) != SCIP_OKAY) \
      throw ipo::SCIPException(_retcode); \
}

namespace ipo
{
   /**
    * \brief Exception handling class for SCIP.
    *
    * Represents a SCIP error in C++.
    **/

   class SCIPException: public std::exception
   {
   public:

      /**
       * \brief Constructs a SCIPException from an error code.
       **/

      SCIPException(SCIP_RETCODE retcode)
         : _retcode(retcode)
      {
         switch (retcode)
         {
         case SCIP_OKAY:
            SCIPsnprintf(_message, 256, "normal termination");
         break;
         case SCIP_ERROR:
            SCIPsnprintf(_message, 256, "unspecified error");
            break;
         case SCIP_NOMEMORY:
            SCIPsnprintf(_message, 256, "insufficient memory error");
            break;
         case SCIP_READERROR:
            SCIPsnprintf(_message, 256, "read error");
            break;
         case SCIP_WRITEERROR:
            SCIPsnprintf(_message, 256, "write error");
            break;
         case SCIP_NOFILE:
            SCIPsnprintf(_message, 256, "file not found error");
            break;
         case SCIP_FILECREATEERROR:
            SCIPsnprintf(_message, 256, "cannot create file");
            break;
         case SCIP_LPERROR:
            SCIPsnprintf(_message, 256, "error in LP solver");
            break;
         case SCIP_NOPROBLEM:
            SCIPsnprintf(_message, 256, "no problem exists");
            break;
         case SCIP_INVALIDCALL:
            SCIPsnprintf(_message, 256, "method cannot be called at this time in solution process");
            break;
         case SCIP_INVALIDDATA:
            SCIPsnprintf(_message, 256, "method cannot be called with this type of data");
            break;
         case SCIP_INVALIDRESULT:
            SCIPsnprintf(_message, 256, "method returned an invalid result code");
            break;
         case SCIP_PLUGINNOTFOUND:
            SCIPsnprintf(_message, 256, "a required plugin was not found");
            break;
         case SCIP_PARAMETERUNKNOWN:
            SCIPsnprintf(_message, 256, "the parameter with the given name was not found");
            break;
         case SCIP_PARAMETERWRONGTYPE:
            SCIPsnprintf(_message, 256, "the parameter is not of the expected type");
            break;
         case SCIP_PARAMETERWRONGVAL:
            SCIPsnprintf(_message, 256, "the value is invalid for the given parameter");
            break;
         case SCIP_KEYALREADYEXISTING:
            SCIPsnprintf(_message, 256, "the given key is already existing in table");
            break;
         case SCIP_MAXDEPTHLEVEL:
            SCIPsnprintf(_message, 256, "maximal branching depth level exceeded");
            break;
         case SCIP_BRANCHERROR:
            SCIPsnprintf(_message, 256, "branching could not be performed (e.g. too large values in variable domain)");
            break;
         default:
            SCIPsnprintf(_message, 256, "unknown error code %d", retcode);
            break;
         }
      }

      /**
       * \brief Destructor.
       **/

      ~SCIPException(void) throw ()
      {

      }

      /**
       * \brief Returns the error message
       **/
      
      const char* what(void) const throw ()
      {
         return _message;
      }

   private:
      /// Buffer for the error message.
      char _message[256];
      /// SCIP error code.
      SCIP_RETCODE _retcode;
   };

   SCIPSolver::SCIPSolver(SCIP* scip)
      : _scip(scip)
   {
      initialize();
   }
   
   SCIPSolver::SCIPSolver(const std::string& fileName)
   {
      SCIP_CALL_EXC( SCIPcreate(&_scip) );
      SCIP_CALL_EXC( SCIPincludeDefaultPlugins(_scip) );
//       SCIP_CALL_EXC(SCIPsetIntParam(_scip, "display/verblevel", 0) );
      SCIP_CALL_EXC( SCIPreadProb(_scip, fileName.c_str(), NULL) );

      initialize();
   }

   void SCIPSolver::initialize()
   {
      _faceConstraints.push_back(nullptr);
      _currentFace = 0;

      std::size_t n = SCIPgetNOrigVars(_scip);
      _variables.resize(n);
      std::vector<std::string> variableNames;
      variableNames.resize(n);
      _instanceObjective = new double[n];
      SCIP_VAR** vars = SCIPgetOrigVars(_scip);
      double scale = SCIPgetObjsense(_scip) == SCIP_OBJSENSE_MAXIMIZE ? 1.0 : -1.0;
      for (std::size_t i = 0; i < n; ++i)
      {
         _variables[i] = vars[i];
         variableNames[i] = SCIPvarGetName(vars[i]);
         _variablesToCoordinates[vars[i]] = i;
         _instanceObjective[i] = scale * SCIPvarGetObj(vars[i]);
      }
      SCIP_CALL_EXC( SCIPsetObjsense(_scip, SCIP_OBJSENSE_MAXIMIZE) );
      
      _space = std::make_shared<Space>(variableNames);
   }

   SCIPSolver::~SCIPSolver()
   {
      delete[] _instanceObjective;
      for (std::size_t i = 1; i < _faceConstraints.size(); ++i)
      {
         SCIPreleaseCons(_scip, &_faceConstraints[i]);
      }
      SCIPfree(&_scip);
      
      assert(_faceConstraints.size() == 1);
   }

   SCIPSolver::Face SCIPSolver::addFace(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
         double* nonzeroCoefficients, double rhs)
   {
      _faceConstraints.push_back(nullptr);
      char consName[16];
      SCIPsnprintf(consName, 16, "face%d", _faceConstraints.size() - 1);
      std::vector<SCIP_VAR*> vars;
      vars.resize(numNonzeros);
      for (std::size_t i = 0; i < numNonzeros; ++i)
         vars[i] = _variables[i];
      SCIP_CALL_EXC( SCIPcreateConsBasicLinear(_scip, &_faceConstraints.back(), consName, numNonzeros,
         &vars[0], &nonzeroCoefficients[0], -SCIPinfinity(_scip), rhs));

      return _faceConstraints.size() - 1;
   }

#ifdef IPO_WITH_GMP

   SCIPSolver::Face SCIPSolver::addFace(std::size_t numNonzeros, std::size_t* nonzeroCoordinates,
         mpq_class* nonzeroCoefficients, const mpq_class& rhs)
   {
      double* approxVector = new double[numNonzeros];
      double approxRhs = rhs.get_d();
   
      Face face = addFace(numNonzeros, nonzeroCoordinates, approxVector, approxRhs);

      delete[] approxVector;
      
      return face;
   }

#endif


   void SCIPSolver::setFace(Face face)
   {
      assert(face >= 0);
      assert(face < _faceConstraints.size());
      
      if (face == _currentFace)
         return;

      if (_currentFace > 0)
      {
         SCIP_CALL_EXC( SCIPdelCons(_scip, _faceConstraints[_currentFace]) );
      }
      if (face > 0)
      {
         SCIP_CALL_EXC( SCIPaddCons(_scip, _faceConstraints[face]) );
      }
      _currentFace = face;
   }

   SCIPOptimizationOracle::SCIPOptimizationOracle(std::shared_ptr<SCIPSolver> solver,
      SCIPSolver::Face face)
      : OptimizationOracle(solver->name()), _solver(solver), _face(face)
   {
      _space = solver->space();
   }
   
   SCIPOptimizationOracle::~SCIPOptimizationOracle()
   {

   }

   bool SCIPOptimizationOracle::isExact() const
   {
#ifdef IPO_WITH_SOPLEX
      return true;
#else
      return false;
#endif /* IPO_WITH_SOPLEX */
   }

   void SCIPOptimizationOracle::maximize(const double* objectiveVector, const Query& query,
      Result& result)
   {
      _solver->setFace(_face);

      SCIP_CALL_EXC( SCIPsetRealParam(_solver->_scip, "limits/time",
         query.timeLimit == std::numeric_limits<double>::infinity() 
         ? SCIPinfinity(_solver->_scip) : query.timeLimit) );

      double oldObjectiveLimit = SCIPgetObjlimit(_solver->_scip);
      int oldSolutionLimit;
      SCIP_CALL_EXC( SCIPgetIntParam(_solver->_scip, "limits/solutions", &oldSolutionLimit) );
      if (query.minObjectiveValue < std::numeric_limits<double>::infinity())
      {
         SCIP_CALL_EXC( SCIPsetObjlimit(_solver->_scip, query.minObjectiveValue) );
         SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "limits/solutions",
            std::max(query.minNumSolutions, 2)) );
      }

      std::size_t n = space()->dimension();

      for (std::size_t i = 0; i < n; ++i)
      {
         SCIP_CALL_EXC( SCIPchgVarObj(_solver->_scip, _solver->_variables[i], objectiveVector[i]) );
      }

      int oldMaxRounds;
      SCIP_CALL_EXC( SCIPgetIntParam(_solver->_scip, "presolving/maxrounds", &oldMaxRounds) );

      for (int attempt = 1; attempt <= 2; ++attempt)
      {
         SCIP_CALL_EXC( SCIPsolve(_solver->_scip) );

         bool hasRay = SCIPhasPrimalRay(_solver->_scip);
         if (hasRay)
         {
            result.firstIndices.push_back(result.nonzeroCoordinates.size());
            result.objectiveValues.push_back(std::numeric_limits<double>::signaling_NaN());
            result.dualBound = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < n; ++i)
            {
               double y = SCIPgetPrimalRayVal(_solver->_scip, _solver->_variables[i]);
               if (!SCIPisZero(_solver->_scip, y))
               {
                  result.nonzeroCoordinates.push_back(i);
                  result.nonzeroValues.push_back(y);
               }
            }
            break;
         }

         SCIP_STATUS status = SCIPgetStatus(_solver->_scip);
         if (status == SCIP_STATUS_UNBOUNDED && attempt == 2 && !hasRay)
         {
            throw std::runtime_error("SCIP reports unboundedness without ray!");
         }
      
         std::size_t numSolutions = SCIPgetNSols(_solver->_scip);
         if (status != SCIP_STATUS_INFEASIBLE && status != SCIP_STATUS_INFORUNBD
            && status != SCIP_STATUS_UNBOUNDED && numSolutions > 0)
         {
            result.dualBound = SCIPgetDualbound(_solver->_scip);
            SCIP_SOL** solutions = SCIPgetSols(_solver->_scip);
            for (std::size_t solIndex = 0; solIndex < numSolutions; ++solIndex)
            {
               SCIP_SOL* sol = solutions[solIndex];
               result.objectiveValues.push_back(SCIPgetSolOrigObj(_solver->_scip, sol));
               result.firstIndices.push_back(result.nonzeroCoordinates.size());
               for (std::size_t i = 0; i < n; ++i)
               {
                  double x = SCIPgetSolVal(_solver->_scip, sol, _solver->_variables[i]);
                  if (!SCIPisZero(_solver->_scip, x))
                  {
                     result.nonzeroCoordinates.push_back(i);
                     result.nonzeroValues.push_back(x);
                  }
               }
            }

            break;
         }
         else if (status == SCIP_STATUS_INFEASIBLE)
         {
            assert(numSolutions == 0);
            result.dualBound = -std::numeric_limits<double>::infinity();
         }

         // Disable presolving for the second round.

         SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", 0) );
         SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
         SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
      }

      SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "presolving/maxrounds", oldMaxRounds) );
      SCIP_CALL_EXC( SCIPfreeSolve(_solver->_scip, true) );
      SCIP_CALL_EXC( SCIPfreeTransform(_solver->_scip) );
      SCIP_CALL_EXC( SCIPsetObjlimit(_solver->_scip, oldObjectiveLimit) );
      SCIP_CALL_EXC( SCIPsetIntParam(_solver->_scip, "limits/solutions", oldSolutionLimit) );

      _solver->setFace(0);
   }

   SCIPSeparationOracle::SCIPSeparationOracle(std::shared_ptr<SCIPSolver> solver,
      SCIPSolver::Face face)
      : SeparationOracle(solver->name()), _solver(solver), _face(face)
   {
      _space = solver->space();
   }
   
   SCIPSeparationOracle::~SCIPSeparationOracle()
   {

   }

   bool SCIPSeparationOracle::separate(const double* vector, bool isPoint, const Query& query,
      Result& result)
   {
      _solver->setFace(_face);

      SCIP_CONS** conss = SCIPgetConss(_solver->_scip);
      std::size_t m = SCIPgetNConss(_solver->_scip);
      for (std::size_t c = 0; c < m; ++c)
      {
         SCIP_CONS* cons = conss[c];
         SCIP_CONSHDLR* conshdlr = SCIPconsGetHdlr(cons);
         std::size_t k = 0;
         SCIP_VAR** vars = nullptr;
         double* vals = nullptr;
         double lhs = -std::numeric_limits<double>::infinity();
         double rhs = std::numeric_limits<double>::infinity();
         const std::string name = SCIPconshdlrGetName(conshdlr);
         if (name == "linear")
         {
            k = SCIPgetNVarsLinear(_solver->_scip, cons);
            vars = SCIPgetVarsLinear(_solver->_scip, cons);
            vals = SCIPgetValsLinear(_solver->_scip, cons);
            lhs = SCIPgetLhsLinear(_solver->_scip, cons);
            rhs = SCIPgetRhsLinear(_solver->_scip, cons);
         }
         else if (name == "setppc")
         {
            k = SCIPgetNVarsSetppc(_solver->_scip, cons);
            vars = SCIPgetVarsSetppc(_solver->_scip, cons);
            lhs = SCIPgetTypeSetppc(_solver->_scip, cons) == SCIP_SETPPCTYPE_PACKING ? 0.0 : 1.0;
            rhs = SCIPgetTypeSetppc(_solver->_scip, cons) == SCIP_SETPPCTYPE_COVERING
               ? SCIPinfinity(_solver->_scip) : 1.0;
         }
         else
         {
            throw std::runtime_error("SCIPSeparationOracle does not implement constraint type <"
               + name + ">.");
         }

         // Set lhs/rhs to 0 if we are separating a ray.
         if (!isPoint && !SCIPisInfinity(_solver->_scip, rhs))
            rhs = 0.0;
         if (!isPoint && !SCIPisInfinity(_solver->_scip, -lhs))
            lhs = 0.0;

         double activity = 0.0;
         for (std::size_t i = 0; i < k; ++i)
         {
            assert(_solver->_variablesToCoordinates.find(vars[i])
               != _solver->_variablesToCoordinates.end());
            activity += vector[_solver->_variablesToCoordinates[vars[i]]]
               * (vals == nullptr ? 1.0 : vals[i]);
         }

         if (!SCIPisFeasLE(_solver->_scip, activity, rhs))
         {
            result.firstIndices.push_back(result.nonzeroCoefficients.size());
            result.rightHandSides.push_back(rhs);
            for (std::size_t i = 0; i < k; ++i)
            {
               result.nonzeroCoordinates.push_back(_solver->_variablesToCoordinates[vars[i]]);
               result.nonzeroCoefficients.push_back(vals == nullptr ? 1.0 : vals[i]);
            }
         }
         else if (!SCIPisFeasGE(_solver->_scip, activity, lhs))
         {
            result.firstIndices.push_back(result.nonzeroCoefficients.size());
            result.rightHandSides.push_back(-lhs);
            for (std::size_t i = 0; i < k; ++i)
            {
               result.nonzeroCoordinates.push_back(_solver->_variablesToCoordinates[vars[i]]);
               result.nonzeroCoefficients.push_back(vals == nullptr ? -1.0 : -vals[i]);
            }
         }
      }

      _solver->setFace(0);

      return !result.nonzeroCoordinates.empty();
   }

} /* namespace ipo */
