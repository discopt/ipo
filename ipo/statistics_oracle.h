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

    /**
     * \brief Resets all values.
     *
     * Resets all values.
     */

    void reset();

    /**
     * \brief Returns how often the oracle was called.
     *
     * Returns how often the oracle was called.
     */

    inline std::size_t numCalls() const
    {
      return _numCalls;
    }

    /**
     * \brief Returns how often the oracle had to forward the call to the next oracle.
     *
     * Returns how often the oracle had to forward the call to the next oracle.
     */

    inline std::size_t numForwards() const
    {
      return _numForwards;
    }

    /**
     * \brief Returns how often the oracle's own answer was satisfactory.
     *
     * Returns how often the oracle's own answer was satisfactory.
     */

    inline std::size_t numSuccess() const
    {
      return _numCalls - _numForwards;
    }

    /**
     * \brief Returns the cumulative time used for this oracle, excluding the time used for forwarded calls.
     *
     * Returns the cumulative time used for this oracle, excluding the time used for forwarded calls.
     */

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

      virtual HeuristicLevel maximizeController(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
        bool& checkDups);

      virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
        const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
        bool& checkDups);

    protected:
      StatisticsOracle& _host;
    };

    virtual HeuristicLevel maximizeController(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
      bool& checkDups);

    void beforeForward();

    void afterForward();

    virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
      bool& checkDups);

  protected:

    std::shared_ptr<OracleBase> _targetOracle;
    std::size_t _numCalls;
    std::size_t _numForwards;
    Timer _timer;
  };

} /* namespace */

#endif /* STATISTICS_ORACLE_H_ */
