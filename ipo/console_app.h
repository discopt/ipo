#ifndef CONSOLE_APP_H_
#define CONSOLE_APP_H_

#include <ostream>
#include <deque>

#include "common.h"
#include "oracles.h"
#include "cache_oracle.h"
#include "projection.h"

namespace ipo {

  /**
   * \brief Base class for an IPO console application.
   *
   * Base class for an IPO console application similar to the ipo binary.
   */

  class ConsoleApplicationBase
  {
  public:
    ConsoleApplicationBase(int numArguments, char** arguments);
    virtual ~ConsoleApplicationBase();

    /**
     * \brief Returns the number of unprocessed arguments.
     *
     * Returns the number of unprocessed arguments.
     */

    inline std::size_t numArguments() const
    {
      return _arguments.size();
    }

    /**
     * \brief Returns the unprocessed argument \c index.
     *
     * Returns the unprocessed argument \c index.
     */

    inline const std::string& argument(std::size_t index) const
    {
      return _arguments[index];
    }

    /**
     * \brief Returns the \c std::deque containing the unprocessed arguments.
     *
     * Returns the \c std::deque containing the unprocessed arguments.
     */

    inline const std::deque<std::string>& arguments() const
    {
      return _arguments;
    }

    /**
     * \brief Sets the basic oracle.
     *
     * Sets the basic oracle. If a restriction to a face, a projection or caching of results are
     * requested, it gets wrapped. The implementation assumes that the latter options were already
     * parsed.
     */

    void setBasicOracle(OracleBase* oracle);

    /**
     * \brief Parses all arguments by repeatedly calling \ref parseArgument().
     *
     * Parses all arguments by repeatedly calling \ref parseArgument().
     * If no arguments are left, the oracle must have been created.
     *
     * \return \c false on error.
     */

    virtual bool parseArguments();

    /**
     * \brief Processes the arguments after oracle creation.
     *
     * Processes the arguments after oracle creation.
     * Most importantly, names of the variables are available,
     * allowing to check the validity of given input data.
     *
     * \return \c false on error.
     */

    virtual bool processArguments();

    /**
     * \brief Runs the requested default tasks.
     *
     * Runs the default tasks that were
     * requested according to the arguments.
     *
     * \return \c false on error.
     */

    virtual bool runDefaultTasks();

    /**
     * \brief Runs the typical methods.
     *
     * Runs the typical methods, i.e., \ref parseArguments(),
     * \ref processArguments() and \ref runDefaultTasks().
     *
     * \return \c false on error.
     */

    virtual bool run();

    /**
     * \brief Prints the usage of the tool.
     *
     * Prints the usage of the tool to stderr.
     */

    virtual void printUsage();



  protected:
    virtual void printAdditionalOptionsPolyhedron(std::ostream& stream);
    virtual void printAdditionalOptionsTasks(std::ostream& stream);
    virtual void printAdditionalOptionsInput(std::ostream& stream);
    virtual void printAdditionalOptionsFurther(std::ostream& stream);
    virtual void printAdditionalOptionsSpecific(std::ostream& stream);

    /**
     * \brief Parses remaining arguments, returning the number of processed arguments.
     *
     * Parses remaining arguments, returning the number of processed arguments.
     * To implement an own IPO console application, a user should implement
     * this method by calling the inherited version first. In case it does not succeed,
     * the implementation should try to identify one option. If an option is recognized,
     * but its potential additional arguments have a wrong syntax, the implementation
     * may raise an exception whose text is then printed before exiting.
     * If no option is recognized, it should check whether all remaining arguments
     * specify a problem instance suitable for the oracle(s).
     *
     * \return
     *   A positive number indicates the number of processed arguments,
     *   a nonpositive number indicates that it could not be parsed.
     */

    virtual int parseArgument(const std::string& firstArgument);

    virtual bool printAmbientDimension();
    virtual bool printVariables();
    virtual bool computeAffineHull(const LinearConstraint& face, const std::string& faceName);
    virtual bool optimizeObjective(const soplex::VectorRational* objective, bool maximize);
    virtual bool generateFacets(const soplex::VectorRational* objective, bool print);
    virtual bool separateRayFacet(const Vector& ray, bool& isFeasible);
    virtual bool separatePointFacet(const Vector& point, bool& isFeasible);
    virtual bool computeSmallestFace(const Vector& point);
    virtual bool printCached();

    void setRelaxationBounds(const soplex::VectorRational& lowerBounds, const soplex::VectorRational& upperBounds);
    void addRelaxationRows(const soplex::LPRowSetRational& rows);

  protected:
    inline bool taskPrintAmbientDimension() const
    {
      return _taskPrintAmbientDimension;
    }

    inline bool taskPrintVariables() const
    {
      return _taskPrintVariables;
    }

    inline bool taskMaximize() const
    {
      return _taskMaximize;
    }

    inline bool taskMinimize() const
    {
      return _taskMinimize;
    }

    inline bool taskDimension() const
    {
      return _taskDimension;
    }

    inline bool taskEquations() const
    {
      return _taskEquations;
    }

    inline bool taskFacet() const
    {
      return _taskSeparateFacet;
    }

    inline bool taskSmallestFace() const
    {
      return _taskSmallestFace;
    }

    inline bool taskFacets() const
    {
      return _taskGenerateFacets;
    }

    inline bool taskPrintCached() const
    {
      return _taskPrintCached;
    }

    inline std::size_t numFaces() const
    {
      return _faces.size();
    }

    inline const LinearConstraint& face(std::size_t index) const
    {
      return _faces[index];
    }

    inline std::size_t numObjectives() const
    {
      return _objectives.size();
    }

    inline const soplex::DVectorRational* objective(std::size_t index) const
    {
      return _objectives[index];
    }

    inline void addObjective(soplex::DVectorRational* objective, const std::string& name)
    {
      _objectives.push_back(objective);
      _objectiveNames.push_back(name);
    }

    inline std::size_t numRays() const
    {
      return _rays.size();
    }

    inline const Vector& ray(std::size_t index) const
    {
      return _rays[index];
    }

    inline std::size_t numPoints() const
    {
      return _points.size();
    }

    inline const Vector& point(std::size_t index) const
    {
      return _points[index];
    }

    inline const Space& space() const
    {
      return _space;
    }

    inline OracleBase* oracle()
    {
      return _oracle;
    }

    inline bool projectedSpace()
    {
      return _projectedOracle != NULL;
    }

    inline ProjectedOracle* projectedOracle()
    {
      return _projectedOracle;
    }

  private:
    std::string _program;
    std::deque<std::string> _arguments;

    std::string _projectionArgument;
    std::string _faceRestrictionArgument;
    std::vector<std::string> _faceArguments;
    std::vector<std::string> _faceFiles;
    std::vector<std::string> _objectiveArguments;
    std::vector<std::string> _objectiveFiles;
    std::vector<std::string> _rayArguments;
    std::vector<std::string> _rayFiles;
    std::vector<std::string> _pointArguments;
    std::vector<std::string> _pointFiles;
    std::size_t _numRandomObjectives;

    bool _taskPrintAmbientDimension;
    bool _taskPrintVariables;
    bool _taskMaximize;
    bool _taskMinimize;
    bool _taskDimension;
    bool _taskEquations;
    bool _taskSeparateFacet;
    bool _taskSmallestFace;
    bool _taskGenerateFacets;
    bool _taskPrintCached;
    bool _optionReadable;
    bool _optionCertificates;
    bool _optionReuseFacets;
    int _optionPrintRandom;
    bool _optionCache;

    Space _space;
    CacheOracle* _cacheOracle;
    Projection* _projection;
    ProjectedOracle* _projectedOracle;
    OracleBase* _oracle;

    std::vector<LinearConstraint> _faces;
    std::vector<std::string> _faceNames;
    std::vector<soplex::DVectorRational*> _objectives;
    std::vector<std::string> _objectiveNames;
    std::vector<Vector> _rays;
    std::vector<std::string> _rayNames;
    std::vector<Vector> _points;
    std::vector<std::string> _pointNames;
    soplex::LPColSetRational _relaxationColumns;
    soplex::LPRowSetRational _relaxationRows;

    soplex::LPRowSetRational* _equations;
    std::vector<Vector> _spanningPoints;
    std::vector<Vector> _spanningRays;
    std::vector<std::size_t> _basicColumns;
  };

} /* namespace ipo */

#endif /* CONSOLE_APP_H_ */

