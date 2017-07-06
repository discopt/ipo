#ifndef IPO_DOMINANT_H_
#define IPO_DOMINANT_H_

#include "oracles.h"

namespace ipo {

  class DominantOracle: public OracleBase
  {
  public:
    /**
     * \brief Constructor for source oracle.
     *
     * Constructor for source \c sourceOracle.
     *
     * \param name         Name of the new oracle.
     * \param sourceOracle Source oracle.
     * \param nextOracle   Next oracle.
     */

    DominantOracle(const std::string& name, const std::shared_ptr<OracleBase>& sourceOracle, const std::shared_ptr<OracleBase>& nextOracle = NULL);

    /**
     * \brief Destructor.
     */

    virtual ~DominantOracle();

    /**
     * \brief Returns the ambient \c space.
     *
     * Returns a const-reference to the ambient \c space.
     */

    inline const Space& space() const
    {
      return _sourceOracle->space();
    }

  protected:

    /**
     * \brief Oracle's implementation to maximize the dense rational \p objective.
     *
     * This method is called by maximizeController() and contains the implementation of the oracle.
     *
     * \param result         After the call, contains the oracle's answer.
     * \param objective      Objective vector \f$ c \in \mathbb{Q}^n \f$ to be maximized.
     * \param objectiveBound Objective value \f$ \gamma \f$ that should be exceeded.
     * \param sort           Set this variable to true if points must be sorted.
     * \param checkDups      Set this variable to true if points or rays must be checked for duplicates.
     *
     * For requirements on the behavior, see Detailed Description of \ref OracleBase.
     */

    virtual HeuristicLevel maximizeImplementation(OracleResult& result, const soplex::VectorRational& objective,
      const ObjectiveBound& objectiveBound, HeuristicLevel minHeuristic, HeuristicLevel maxHeuristic, bool& sort,
      bool& checkDups);

  protected:

    std::shared_ptr<OracleBase> _sourceOracle; // Source oracle.
    bool _isFeasible;
  };

} /* namespace ipo */



#endif /* IPO_DOMINANT_H_ */
