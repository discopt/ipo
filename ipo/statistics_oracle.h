#ifndef STATISTICS_ORACLE_H_
#define STATISTICS_ORACLE_H_

#include "common.h"
#include "timer.h"
#include "oracles.h"

namespace ipo {

  class StatisticsOracle : public OracleBase
  {
  public:
    /**
     * \brief Constructor.
     * 
     * Constructor. The \p targetOracle is the one for which this wrapper collects statistics.
     */

    StatisticsOracle(const std::shared_ptr<OracleBase>& targetOracle);

    /**
     * \brief Destructor.
     * 
     * Destructor.
     */

    virtual ~StatisticsOracle();

    /**
     * \brief Restricts the oracle to the face defined by \p newFace.
     *
     * Restricts the optimization oracle to the face \f$ F \f$ of \f$ P \f$ defined by \p newFace.
     * For \p newFace equal to \c NULL we define \f$ F := P \f$.
     *
     * This implementation forwards it to the target oracle.
     */    
    
    virtual void setFace(const LinearConstraint& newFace);
    
    /**
     * \brief Returns the wrapped oracle.
     * 
     * Returns the wrapped oracle.
     */

    inline const std::shared_ptr<OracleBase>& targetOracle() const
    {
      return _targetOracle;
    }

    inline std::size_t numCalls() const
    {
      return _numCalls;
    }

    inline std::size_t numForwards() const
    {
      return _numForwards;
    }

    inline std::size_t numSuccess() const
    {
      return _numCalls - _numForwards;
    }

    inline double time() const
    {
      return _timer.time();
    }

  protected:
    class ForwardWrapperOracle: public OracleBase
    {
    public:
      ForwardWrapperOracle(StatisticsOracle& host, const std::shared_ptr<OracleBase>& nextOracle);

      virtual ~ForwardWrapperOracle();

    protected:

      virtual std::size_t maximizeController(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);

      virtual std::size_t maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);

    protected:
      StatisticsOracle& _host;
    };

    virtual std::size_t maximizeController(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);

    void beforeForward();

    void afterForward();

    virtual std::size_t maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, std::size_t minHeuristic, std::size_t maxHeuristic, bool& sort, bool& checkDups);
  
  protected:

    std::shared_ptr<OracleBase> _targetOracle;
    std::size_t _numCalls;
    std::size_t _numForwards;
    Timer _timer;
  };
  
} /* namespace */

#endif /* STATISTICS_ORACLE_H_ */