#include <iostream>
#include <cmath>
#include <sstream>
#include <map>

#include <ipo/ipo.h>

#ifdef WITH_SCIP
#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include <scip/cons_linear.h>
#include <ipo/scip_exception.hpp>
#include <ipo/scip_oracles.h>
#endif

#include "ipo/affine_hull.h"
#include "ipo/facets.h"
#include "ipo/smallest_face.h"
#include "ipo/min_norm_2d.h"
#include "ipo/console_app.h"

using namespace ipo;

#ifdef WITH_SCIP
#define ORACLE_DEFAULT_SCIP
#define ORACLE_DEFAULT "scip"
#elif WITH_EXACT_SCIP
#define ORACLE_DEFAULT_EXACT_SCIP
#define ORACLE_DEFAULT "exactscip"
#else
#define ORACLE_DEFAULT_NONE
#define ORACLE_DEFAULT ""
#endif

class IPOConsoleApplication: public ConsoleApplicationBase
{
public:
  IPOConsoleApplication(int numArguments, char** arguments) :
      ConsoleApplicationBase(numArguments, arguments)
  {
    _oracleArgument = ORACLE_DEFAULT;
    _heuristicArgument = "";
    _relaxationArgument = false;
    _taskPrintInstanceObjective = false;
    _useInstanceObjective = true;
    _useInstanceBounds = true;
    _useInstanceInequalities = true;

#ifdef WITH_SCIP
    _scipSpace = NULL;
    _scipMip = NULL;
    _scipOracle = NULL;
    _scipCorrectorOracle = NULL;
#endif
  }

  virtual ~IPOConsoleApplication()
  {
    /// ~ConsoleApplicationBase will free the oracle passed to setBasicOracle.

#ifdef WITH_SCIP
    if (_scipMip)
      delete _scipMip;
    if (_scipOracle)
      delete _scipOracle;
    if (_scipCorrectorOracle)
      delete _scipCorrectorOracle;
    if (_scipSpace)
      delete _scipSpace;
#endif
  }

  virtual bool run()
  {
    if (!ConsoleApplicationBase::run())
      return false;
    if (_taskPrintInstanceObjective)
    {
      if (!printInstanceObjective())
        return false;
    }
    return true;
  }

  virtual void printAdditionalOptionsPolyhedron(std::ostream& stream)
  {
    std::cerr << "  --oracle ORACLE\n";
    std::cerr << "      Use oracle ORACLE to define the polyhedron P.\n";
    std::cerr << "  --heuristic ORACLE\n";
    std::cerr << "      Use oracle ORACLE as a heuristic in addition to the oracle. See section on oracles below.\n";
    std::cerr << "  --relaxation\n";
    std::cerr << "      If the instance is a MIP, set all variables to continuous.\n";
  }

  virtual void printAdditionalOptionsTasks(std::ostream& stream)
  {
    std::cerr << "  --print-instance-objective\n";
    std::cerr << "      If the instance is a MIP, print its maximization objective.\n";
  }

  virtual void printAdditionalOptionsFurther(std::ostream& stream)
  {
    std::cerr << "  --instance-objective on|off\n";
    std::cerr << "      If on (default) and the instance is a MIP, it considers the maximization objective of the\n";
    std::cerr << "      instance for facet-generation.\n";
    std::cerr << "  --use-bounds on|off\n";
    std::cerr << "      If on (default) and the instance is a MIP, it uses its bounds for facet-generation.\n";
    std::cerr << "  --use-inequalities on|off\n";
    std::cerr << "      If on (default) and the instance is a MIP, it uses its inequalities for facet-generation.\n";
  }

  virtual void printAdditionalOptionsSpecific(std::ostream& stream)
  {
    std::cerr << "\n";
    std::cerr << "Oracles that can be used for --oracle and --heuristic (see IPO's build options):\n";
#ifdef WITH_SCIP
    std::cerr << "  scip      Use the MIP solver SCIP from the SCIP Optimization Suite (default oracle).\n";
#endif
#ifdef WITH_EXACT_SCIP
    std::cerr << "  exactscip Use the exact MIP solver SCIP-ex from the SCIP Optimization Suite";
#ifdef ORACLE_DEFAULt_EXACT_SCIP
    std::cerr << " (default oracle)";
#endif
    std::cerr << ".\n";
#endif
    std::cerr << "  PROGRAM   For each oracle call, execute PROGRAM and communicate via stdin and stdout.\n";
    std::cerr << "            See IPO's python directory for example implementations.\n";
#ifdef ORACLE_DEFAULT_NONE
    std::cerr << "IPO is built without a default oracle (e.g., SCIP). See build options!\n";
#endif
  }

  virtual int parseArgument(const std::string& firstArgument)
  {
    int parsed = ConsoleApplicationBase::parseArgument(firstArgument);
    if (parsed > 0)
      return parsed;

    if (firstArgument == "--oracle")
    {
      if (numArguments() == 1)
        throw std::runtime_error("Invalid option: --oracle must be followed by an oracle name.");
      _oracleArgument = argument(1);

#ifdef WITH_SCIP
      if (_oracleArgument == "scip")
      {
        return 2;
      }
#endif
#ifdef WITH_EXACT_SCIP
      if (_oracleArgument == "exactscip")
      {
        return 2;
      }
#endif
      std::ifstream file(_oracleArgument.c_str());
      if (!file.is_open())
      {
        std::stringstream ss;
        ss << "Invalid option: " << _oracleArgument << " is not a valid oracle.";
        throw std::runtime_error(ss.str());
      }
    }
    if (firstArgument == "--heuristic")
    {
      if (numArguments() == 1)
        throw std::runtime_error("Invalid option: --heuristic must be followed by an oracle name.");
      _heuristicArgument = argument(1);

#ifdef WITH_SCIP
      if (_heuristicArgument == "scip")
      {
        return 2;
      }
#endif
#ifdef WITH_EXACT_SCIP
      if (_heuristicArgument == "exactscip")
      {
        return 2;
      }
#endif
      std::ifstream file(_heuristicArgument.c_str());
      if (!file.is_open())
      {
        std::stringstream ss;
        ss << "Invalid option: " << _heuristicArgument << " is not a valid oracle.";
        throw std::runtime_error(ss.str());
      }
    }
    if (firstArgument == "--relaxation")
    {
      _relaxationArgument = true;
      return 1;
    }
    if (firstArgument == "--print-instance-objective")
    {
      _taskPrintInstanceObjective = true;
      return 1;
    }
    if (firstArgument == "--instance-objective")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _useInstanceObjective = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _useInstanceObjective = false;
          return 2;
        }
      }
      throw std::runtime_error("Invalid option: --instance-objective must be followed by either \"on\" or \"off\".");
    }
    if (firstArgument == "--use-bounds")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _useInstanceBounds = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _useInstanceBounds = false;
          return 2;
        }
      }
      throw std::runtime_error("Invalid option: --use-bounds must be followed by either \"on\" or \"off\".");
    }
    if (firstArgument == "--use-inequalities")
    {
      if (numArguments() > 1)
      {
        if (argument(1) == "on")
        {
          _useInstanceInequalities = true;
          return 2;
        }
        if (argument(1) == "off")
        {
          _useInstanceInequalities = false;
          return 2;
        }
      }
      throw std::runtime_error("Invalid option: --use-inequalities must be followed by either \"on\" or \"off\".");
    }

    /// Now the remaining arguments must be an instance.

    if (numArguments() > 1
#ifdef WITH_SCIP
        && _oracleArgument == "scip"
#endif
            )
    {
      std::stringstream ss;
      ss << "Expected one instance argument, but " << numArguments() << " non-options remain:\n";
      for (std::size_t i = 0; i < numArguments(); ++i)
        ss << "\n  " << argument(i);
      throw std::runtime_error(ss.str());
    }

    OracleBase* oracle = NULL;
#ifdef WITH_SCIP
    if (_oracleArgument == "scip")
    {
      try
      {
        SCIP* scip = NULL;
        try
        {
          SCIP_CALL_EXC(SCIPcreate(&scip));
          SCIP_CALL_EXC(SCIPincludeDefaultPlugins(scip));
          SCIP_CALL_EXC(SCIPsetIntParam(scip, "display/verblevel", 0));
          SCIP_CALL_EXC(SCIPreadProb(scip, firstArgument.c_str(), NULL));
          SCIP_CALL_EXC(SCIPtransformProb(scip));
        }
        catch (SCIPException& ex)
        {
          SCIP_CALL_EXC(SCIPfree(&scip));
          throw ex;
        }

        _scipSpace = new Space;
        _scipMip = new MixedIntegerProgram(*_scipSpace, scip);
        SCIP_CALL_EXC(SCIPfree(&scip));
      }
      catch (SCIPException& exc)
      {
        std::stringstream ss;
        ss << "Invalid SCIP instance \"" << firstArgument << "\".";
        throw std::runtime_error(ss.str());
      }

      /// Extract instance objective.

      assert(_instanceObjective.size() == 0);
      for (std::size_t i = 0; i < _scipSpace->dimension(); ++i)
      {
        const soplex::Rational& x = _scipMip->columns().maxObj(i);
        if (x != 0)
          _instanceObjective.add(i, x);
      }

      _scipOracle = new SCIPOracle(firstArgument, *_scipMip);
      _scipCorrectorOracle = new MixedIntegerProgramCorrectorOracle(firstArgument + "-corrected", *_scipMip,
        _scipOracle);
      oracle = _scipCorrectorOracle;
    }
#endif

    if (_heuristicArgument != "")
    {
      throw std::runtime_error("Heuristic is not implemented!");
    }

    if (oracle == NULL)
    {
      std::stringstream ss;
      ss << "Unable to create oracle \"" << _oracleArgument << "\"." << std::endl;
      throw std::runtime_error(ss.str());
    }

    setBasicOracle(oracle);

    return 1;
  }

  virtual bool processArguments()
  {
    if (!ConsoleApplicationBase::processArguments())
      return false;

    /// TODO: Instance objective may be invalid if projection is enabled!
    if (_useInstanceObjective)
    {
      soplex::DSVectorRational* obj = new soplex::DSVectorRational;
      *obj = _instanceObjective;
      addObjective(obj, "instance");
    }

#ifdef WITH_SCIP
    if (_useInstanceBounds)
    {
      soplex::DVectorRational lower(space().dimension());
      soplex::DVectorRational upper(space().dimension());
      const soplex::LPColSetRational& cols = _scipMip->columns();
      setRelaxationBounds(cols.lower(), cols.upper());
    }
    if (_useInstanceInequalities)
    {
      addRelaxationRows(_scipMip->rows());
    }
#endif

    return true;
  }

protected:

  virtual bool printInstanceObjective()
  {
    std::cout << "Instance objective:\n ";
    space().printLinearForm(std::cout, &_instanceObjective);
    std::cout << "\n" << std::flush;
    return true;
  }

  std::string _oracleArgument;
  std::string _heuristicArgument;

  bool _relaxationArgument;
  bool _taskPrintInstanceObjective;
  bool _useInstanceObjective;
  bool _useInstanceBounds;
  bool _useInstanceInequalities;

#ifdef WITH_SCIP
  Space* _scipSpace;
  MixedIntegerProgram* _scipMip;
  SCIPOracle* _scipOracle;
  MixedIntegerProgramCorrectorOracle* _scipCorrectorOracle;
#endif
  soplex::DSVectorRational _instanceObjective;
};

int main(int argc, char** argv)
{
  IPOConsoleApplication app(argc, argv);
  app.run();
  soplex::Rational::freeListMem();

  return EXIT_SUCCESS;
}

