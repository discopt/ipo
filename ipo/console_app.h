#ifndef CONSOLE_APP_H_
#define CONSOLE_APP_H_

#include <ostream>
#include <deque>

#include "ipo.h"
#include "unique_rational_vectors.h"
#include "oracles.h"
#include "spx_gmp.h"

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
     * Sets the basic oracle. If a restriction to a face or a projection
     * are requested, it gets wrapped.
     * The implementation assumes that the latter options were already parsed.
     */

    void setBasicOracle(FaceOptimizationOracleBase* oracle);

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
    virtual bool computeAffineHull(Face* face, const std::string& faceName);
    virtual bool optimizeObjective(const soplex::SVectorRational* objective, bool maximize);
    virtual bool generateFacets(const soplex::SVectorRational* objective, bool print);
    virtual bool separateDirectionFacet(const Direction* direction, bool& isFeasible);
    virtual bool separatePointFacet(const Point* point, bool& isFeasible);
    virtual bool computeSmallestFace(const Point* point);
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

    inline const Face* face(std::size_t index) const
    {
      return _faces[index];
    }

    inline std::size_t numObjectives() const
    {
      return _objectives.size();
    }

    inline const soplex::SVectorRational* objective(std::size_t index) const
    {
      return _objectives[index];
    }

    inline void addObjective(soplex::DSVectorRational* objective, const std::string& name)
    {
      _objectives.push_back(objective);
      _objectiveNames.push_back(name);
    }

    inline std::size_t numDirections() const
    {
      return _directions.size();
    }

    inline const Direction* directions(std::size_t index) const
    {
      return _directions[index];
    }

    inline std::size_t numPoints() const
    {
      return _points.size();
    }

    inline const Point* points(std::size_t index) const
    {
      return _points[index];
    }
    
    inline const Space& space() const
    {
      return _space;
    }

    inline FaceOptimizationOracleBase* oracle()
    {
      return _oracle;
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
    std::vector<std::string> _directionArguments;
    std::vector<std::string> _directionFiles;
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
    bool _optionPrintRandom;

    Space _space;
    FaceOptimizationOracleBase* _oracle;

    std::vector<Face*> _faces;
    std::vector<std::string> _faceNames;
    std::vector<soplex::DSVectorRational*> _objectives;
    std::vector<std::string> _objectiveNames;
    std::vector<Direction*> _directions;
    std::vector<std::string> _directionNames;
    std::vector<Point*> _points;
    std::vector<std::string> _pointNames;
    soplex::LPColSetRational _relaxationColumns;
    soplex::LPRowSetRational _relaxationRows;

    UniqueRationalVectors* _cachedPoints;
    UniqueRationalVectors* _cachedDirections;
    soplex::LPRowSetRational* _equations;
    VectorSubset _spanningCachedPoints;
    VectorSubset _spanningCachedDirections;
    VectorSubset _basicColumns;
  };

} /* namespace ipo */

#endif /* CONSOLE_APP_H_ */

