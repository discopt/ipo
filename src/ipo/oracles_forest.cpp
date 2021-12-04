// #define IPO_DEBUG /* Uncomment to debug this file. */

#include <ipo/oracles_forest.hpp>

#include <sstream>

#include <iostream>

namespace ipo
{
  std::size_t getRepresentative(std::vector<std::size_t>& nodesComponentRoot, std::size_t node)
  {
    assert(node < nodesComponentRoot.size());

    std::size_t current = node;
    std::size_t next;
    while ((next = nodesComponentRoot[current]) != std::numeric_limits<std::size_t>::max())
      current = next;
    std::size_t root = current;
    current = node;
    while ((next = nodesComponentRoot[current]) != std::numeric_limits<std::size_t>::max())
    {
      if (next != root)
        nodesComponentRoot[current] = root;
      current = next;
    }
    return root;
  }

  template <typename NumberType>
  ForestOptimizationOracle<NumberType>::ForestOptimizationOracle(std::size_t numNodes,
    std::pair<std::size_t, std::size_t>* edgesFirst, std::pair<std::size_t, std::size_t>* edgesBeyond, bool spanning,
    const std::string& name)
    : OptimizationOracle<NumberType>(name), _spanning(spanning), _edges(edgesBeyond - edgesFirst)
  {
    _numNodes = numNodes;
    std::vector<std::string> variableNames(edgesBeyond - edgesFirst);
    for (size_t e = 0; e < _edges.size(); ++e)
    {
      assert(edgesFirst->first < numNodes);
      assert(edgesFirst->second < numNodes);
      _edges[e].first = edgesFirst->first;
      _edges[e].second = edgesFirst->second;
      std::ostringstream ss;
      ss << "x_" << edgesFirst->first << "#" << edgesFirst->second;
      variableNames[e] = ss.str();
      ++edgesFirst;
    }
    this->_space = std::make_shared<Space>(variableNames);
  }

#if defined(IPO_WITH_GMP)

  template <>
  IPO_EXPORT
  ForestOptimizationOracle<mpq_class>::ForestOptimizationOracle(std::size_t numNodes,
    std::pair<std::size_t, std::size_t>* edgesFirst, std::pair<std::size_t, std::size_t>* edgesBeyond, bool spanning,
    const std::string& name)
    : OptimizationOracle<mpq_class>(name), _spanning(spanning), _edges(edgesBeyond - edgesFirst)
  {
    _numNodes = numNodes;
    std::vector<std::string> variableNames(edgesBeyond - edgesFirst);
    for (size_t e = 0; e < _edges.size(); ++e)
    {
      assert(edgesFirst->first < numNodes);
      assert(edgesFirst->second < numNodes);
      _edges[e].first = edgesFirst->first;
      _edges[e].second = edgesFirst->second;
      std::ostringstream ss;
      ss << "x_" << edgesFirst->first << "#" << edgesFirst->second;
      variableNames[e] = ss.str();
      ++edgesFirst;
    }
    this->_space = std::make_shared<Space>(variableNames);
  }

#endif /* IPO_WITH_GMP */

  template <typename NumberType>
  ForestOptimizationOracle<NumberType>::~ForestOptimizationOracle()
  {

  }

#if defined(IPO_WITH_GMP)

  template <>
  IPO_EXPORT
  ForestOptimizationOracle<mpq_class>::~ForestOptimizationOracle()
  {

  }

#endif /* IPO_WITH_GMP */

  template <typename NumberType>
  OptimizationResponse<NumberType> ForestOptimizationOracle<NumberType>::maximize(const NumberType* objectiveVector,
    const OptimizationQuery<NumberType>& query)
  {
    std::vector<std::pair<NumberType, std::size_t>> weights(_edges.size());
    for (size_t e = 0; e < _edges.size(); ++e)
      weights[e] = std::make_pair(objectiveVector[e], e);

    // Initialize singleton components.
    std::vector<std::size_t> nodesComponentRoot(_numNodes, std::numeric_limits<std::size_t>::max());

    // Sort edges by weight.
    std::sort(weights.begin(), weights.end(),
      [](const std::pair<NumberType, std::size_t>& a, const std::pair<NumberType, std::size_t>& b) -> bool
      {
        return a.first > b.first;
      }
    );

    // Go through edges by descending weight and try to merge components.
    std::vector<std::pair<std::size_t, NumberType>> optimalSolution;
    NumberType optimum = 0;
    for (size_t i = 0; i < weights.size(); ++i)
    {
      std::size_t edge = weights[i].second;
      std::size_t u = getRepresentative(nodesComponentRoot, _edges[edge].first);
      std::size_t v = getRepresentative(nodesComponentRoot, _edges[edge].second);

      // If endnodes are the same component, we skip it.
      if (u == v)
        continue;

      auto weight = weights[i].first;
      if (!_spanning && weight < 0)
        continue;

      nodesComponentRoot[u] = v;
      optimalSolution.push_back(std::make_pair(edge, 1.0));
      optimum += weight;
    }

    OptimizationResponse<NumberType> response;
    response.outcome = OptimizationOutcome::FEASIBLE;
    response.primalBound = optimum;
    response.hasPrimalBound = true;
    response.dualBound = optimum;
    response.hasDualBound = true;
    response.points.push_back(typename OptimizationResponse<NumberType>::Point(
      std::make_shared<sparse_vector<NumberType>>(optimalSolution, true), optimum));

    return response;
  }

  template class ForestOptimizationOracle<double>;

#if defined(IPO_WITH_GMP)
  template class ForestOptimizationOracle<mpq_class>;
#endif /* IPO_WITH_GMP */  
  
} /* namespace ipo */
