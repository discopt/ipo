#pragma once

#include <iostream>
#include <sstream>
#include <random>
#include <unordered_map>
#include <chrono>
#include <ctime>
#include <ratio>

#include <rapidxml.hpp>
#include <rapidxml_utils.hpp>
#include <rapidxml_print.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <ipo/oracles.hpp>
#include <ipo/affine_hull.hpp>
#include <ipo/lp.hpp>
#include <ipo/oracles_radialcone.hpp>

namespace inverse
{
  template <typename Solver, typename Number>
  struct InverseMIP
  {
    std::string name;
    std::shared_ptr<ipo::Space> space;
    std::shared_ptr<sparse_vector<Number>> targetObjective;
    std::vector<std::shared_ptr<Solver>> solvers;
    std::vector<std::shared_ptr<sparse_vector<Number>>> targetSolutions;

    InverseMIP(const std::string& initialName)
      : name(initialName)
    {

    }

    InverseMIP(InverseMIP&& other)
      : space(other.space), targetObjective(std::move(other.targetObjective)),
      solvers(std::move(other.solvers)), targetSolutions(std::move(other.targetSolutions))
    {

    }
  };

  template <typename Number>
  std::shared_ptr<sparse_vector<Number>> parseVector(std::shared_ptr<ipo::Space> space, rapidxml::xml_node<>* node,
    const char* value_attribute)
  {
    std::unordered_map<std::string, std::size_t> namesToVariables;
    for (std::size_t v = 0; v < space->dimension(); ++v)
      namesToVariables[space->variable(v)] = v;

    std::vector<typename sparse_vector<Number>::value_type> values;

    for (auto child = node->first_node("variable"); child; child = child->next_sibling("variable"))
    {
      auto name_attr = child->first_attribute("name");
      if (name_attr == nullptr)
        throw std::runtime_error("A variable-node is missing a name-attribute.");

      auto iter = namesToVariables.find(name_attr->value());
      if (iter == namesToVariables.end())
        throw std::runtime_error(std::string("A variable-node has unknown name \"") + name_attr->value() + "\".");
      std::size_t coordinate = iter->second;

      auto value_attr = child->first_attribute(value_attribute);
      if (value_attr == nullptr)
      {
        throw std::runtime_error(std::string("A variable-node \"") + name_attr->value()
          + "\" is missing a " + value_attribute + "-attribute.");
      }
      std::istringstream str(value_attr->value());
      double value;
      str >> value;

      values.push_back(std::make_pair(coordinate, value));
    }

    return std::make_shared<sparse_vector<Number>>(std::move(values), true);
  }

  template <typename Solver, typename Number>
  InverseMIP<Solver, Number> readInverseProblem(const std::string& instanceFileName)
  {
    InverseMIP<Solver, Number> result(instanceFileName);
    std::string fileData;

    if (instanceFileName.substr(std::max(3UL, instanceFileName.length()) - 3) == ".gz")
    {
      std::ifstream file(instanceFileName, std::ios_base::in | std::ios_base::binary);
      boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
      in.push(boost::iostreams::gzip_decompressor());
      in.push(file);
      std::istream incoming(&in);
      fileData = std::string(std::istreambuf_iterator<char>(incoming), {});
    }
    else
    {
      rapidxml::file<> xmlFile(instanceFileName.c_str());
      fileData = xmlFile.data();
    }

    rapidxml::xml_document<> xml_doc;
    xml_doc.parse<0>(const_cast<char*>(fileData.c_str()));

    rapidxml::xml_node<>* xml_inverse_mip = xml_doc.first_node("inverse-mip");
    rapidxml::xml_attribute<>* xml_attr = nullptr;
    xml_attr = xml_inverse_mip->first_attribute("dimension");
    std::size_t dimension = std::numeric_limits<std::size_t>::max();
    if (xml_attr)
    {
      std::istringstream str(xml_attr->value());
      str >> dimension;
    }

    std::cout << "Inverse problem in dimension " << dimension << "." << std::endl;

    for (auto *xml_mip = xml_inverse_mip->first_node("mip"); xml_mip; xml_mip = xml_mip->next_sibling("mip"))
    {
      xml_attr = xml_mip->first_attribute("file");
      if (!xml_attr)
        throw std::runtime_error("Inverse MIP's mip node has no 'file' attribute.");

      std::string fileName = xml_attr->value();

      std::shared_ptr<Solver> solver = std::make_shared<Solver>(fileName);
      result.space = solver->space();

      if (result.space->dimension() != dimension)
      {
        throw std::runtime_error("Dimension of mip in file \"" + fileName
          + "\" does not match the one from the inverse MIP.");
      }

      result.solvers.push_back(solver);
      result.targetSolutions.push_back(parseVector<Number>(result.space,
        xml_mip->first_node("target-solution"), "value"));
    }

    auto *xml_target_objective = xml_inverse_mip->first_node("target-objective");
    xml_attr = xml_target_objective->first_attribute("norm");
    std::string target_objective_norm = xml_attr ? xml_attr->value() : "";

    result.targetObjective = parseVector<Number>(result.space, xml_target_objective, "coefficient");

    return result;
  }

  template<typename Solver, typename Number>
  void solve(const InverseMIP<Solver, Number>& instance)
  {
    std::vector<std::size_t> nonzeroColumns;
    std::vector<Number> nonzeroCoefficients;

    std::cout << "Solving inverse MIP <" << instance.name << ">." << std::endl;

    std::vector<std::shared_ptr<ipo::Polyhedron<Number>>> polyhedra;
    std::shared_ptr<ipo::Space> space;
    for (std::size_t i = 0; i < instance.solvers.size(); ++i)
    {
      auto oracle = instance.solvers[i]->template getOptimizationOracle<Number>();
      std::cout << "Creating RadialConeOptimizationOracle..." << std::endl;
      auto radialConeOracle = std::make_shared<ipo::RadialConeOptimizationOracle<Number>>(oracle,
        instance.targetSolutions[i]);

      if (radialConeOracle->isTrustRegionCapable())
      {
        std::cout << "Activating trust region for " << oracle->name() << "." << std::endl;
        radialConeOracle->enableManhattanTrustRegion(16, oracle->space()->dimension(), 2, 0);
      }

      polyhedra.push_back( std::make_shared<ipo::Polyhedron<Number>>(radialConeOracle));
      space = oracle->space();
    }

    std::size_t n = space->dimension();

    ipo::LP<Number> lp;
    lp.setSense(ipo::LPSense::MINIMIZE);
    for (std::size_t v = 0; v < n; ++v)
      lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), 0, space->variable(v));
    std::size_t firstTargetSolutionValueColumn = lp.numColumns();
    for (std::size_t p = 0; p < polyhedra.size(); ++p)
    {
      std::stringstream stream;
      stream << "targetsolval#" << (p+1);
      lp.addColumn(lp.minusInfinity(), lp.plusInfinity(), 0, stream.str());
    }
    std::size_t firstNormColumn = lp.numColumns();
    for (std::size_t v = 0; v < n; ++v)
      lp.addColumn(0, lp.plusInfinity(), 1, "obj#" + space->variable(v));

    std::vector<Number> targetObjectiveDense(n, 0);
    for (auto iter : *instance.targetObjective)
      targetObjectiveDense[iter.first] = iter.second;
    for (std::size_t v = 0; v < n; ++v)
    {
      // x_i - targetobj_i <= obj#x_i   <=> -x_i + obj#x_i >= -c_i
      Number coefs[2] = { -1, 1 };
      std::size_t columns[2] = { v, firstNormColumn + v };
      lp.addRow(-targetObjectiveDense[v], 2, columns, coefs, lp.plusInfinity(), "norm1#" + space->variable(v));

      // -x_i + targetobj_i <= obj#x_i  <=> x_i + obj#x_i >= c_i
      coefs[0] = 1;
      lp.addRow(targetObjectiveDense[v], 2, columns, coefs, lp.plusInfinity(), "norm2#" + space->variable(v));
    }

    // y_p = \sum tgtsol^p_i x_i
    for (std::size_t p = 0; p < polyhedra.size(); ++p)
    {
      nonzeroColumns.clear();
      nonzeroCoefficients.clear();

      for (auto iter : *instance.targetSolutions[p])
      {
        nonzeroColumns.push_back(iter.first);
        nonzeroCoefficients.push_back(iter.second);
      }
      nonzeroColumns.push_back(firstTargetSolutionValueColumn + p);
      nonzeroCoefficients.push_back(-1);
      std::stringstream stream;
      stream << "targetsolval#" << (p+1);
      lp.addRow(0, nonzeroCoefficients.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0, stream.str());
    }

    std::size_t iteration = 0;
    double timeTotalLP = 0.0;
    double timeTotalOracles = 0.0;
    double timeTotalIterations = 0.0;
    while (true)
    {
      const auto iteration_start = std::chrono::high_resolution_clock::now();

      ++iteration;

      auto status = lp.solve();
      timeTotalLP += lp.getSolveTime();

      if (status == ipo::LPStatus::OPTIMAL)
      {
        std::cout.setf( std::ios_base::fmtflags(), std::ios_base::floatfield );
        std::cout << "LP #" << iteration << " with " << lp.numColumns() << " variables and " << lp.numRows() << " rows solved."
          << " Optimum is " << ipo::formatNumberApprox(lp.getObjectiveValue()) << ". Time: "
          << lp.getSolveTime() << "s. Total LP time: " << timeTotalLP << "s." << std::endl;

        assert(lp.hasPrimalSolution());
        std::vector<Number> solutionObjective = lp.getPrimalSolution();
        // for (std::size_t v = 0; v < n; ++v)
        // {
        //   if (solutionObjective[v] || targetObjectiveDense[v])
        //   {
        //     std::cout << "  Variable #" << v << " " << space->variable(v) << ": candidate = " << solutionObjective[v]
        //       << ", target = " << targetObjectiveDense[v] << ", |diff| = "
        //       << fabs(ipo::convertNumber<double>(solutionObjective[v] - targetObjectiveDense[v])) << std::endl;
        //   }
        // }

        std::size_t numAddedCuts = 0;
        for (std::size_t p = 0; p < polyhedra.size(); ++p)
        {
          auto poly = polyhedra[p];
          auto targetSolution = instance.targetSolutions[p];

          Number targetSolutionValue = *targetSolution * solutionObjective;
          std::vector<Number> targetSolutionDense(n, 0);
          for (auto iter : *targetSolution)
            targetSolutionDense[iter.first] = iter.second;
          std::size_t targetSolutionSize = targetSolution->size();

          ipo::OptimizationQuery<Number> optQuery;
          // optQuery.setMinPrimalBound(targetSolutionValue);
          std::cout << "Target solution has objective value " << ipo::formatNumberApprox(targetSolutionValue)
            << std::endl;

          std::cout << optQuery << std::endl;

          const auto start = std::chrono::high_resolution_clock::now();
          ipo::OptimizationResponse<Number> optResponse = poly->maximize(&solutionObjective[0], optQuery);
          const auto end = std::chrono::high_resolution_clock::now();

          double timeOracle = 1.e-3 * std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
          timeTotalOracles += timeOracle;

          std::cout << optResponse << " in " << timeOracle << "s. Total oracle time: " << timeTotalOracles << "s."
            << std::endl;

          if (optResponse.hasPrimalBound() && optResponse.primalBound() > targetSolutionValue)
          {

          }

          for (const auto& point : optResponse.points)
          {
            if (point.objectiveValue <= targetSolutionValue)
              continue;

            if (point.objectiveValue > targetSolutionValue)
              throw std::runtime_error("RadialConeOptimizationOracle returned point better than apex.");
          }
          for (const auto& ray : optResponse.rays)
          {
            nonzeroColumns.clear();
            nonzeroCoefficients.clear();

            for (const auto& iter : *ray.vector)
            {
              nonzeroColumns.push_back(iter.first);
              nonzeroCoefficients.push_back(iter.second);
            }

            lp.addRow(lp.minusInfinity(), nonzeroColumns.size(), &nonzeroColumns[0], &nonzeroCoefficients[0], 0);
            ++numAddedCuts;
          }
        }

        if (numAddedCuts == 0)
        {
          std::cout << "All target solutions are optimal for proposed objective." << std::endl;
          break;
        }
      }
      else
      {
        std::cout << "LP status is " << status << std::endl;
        break;
      }

      const auto iteration_end = std::chrono::high_resolution_clock::now();
      double timeIteration = 1.e-3 * std::chrono::duration_cast<std::chrono::milliseconds>(iteration_end - iteration_start).count();
      timeTotalIterations += timeIteration;
      std::cout << "Iteration time: " << timeIteration << "s. Total iteration time: " << timeTotalIterations << "s.\n" << std::endl;
    }
  }

  int printUsage(const std::string& program)
  {
    std::cout << program << " [OPTIONS] FILE...\n";
    std::cout << "Solves inverse optimization problems on a polyhedron defined by FILE.\n";
    std::cout << "General options:\n";
    std::cout << " -h       Show this help and exit.\n";
    std::cout << "Oracle/polyhedron options:\n";
#if defined(IPO_WITH_RATIONAL)
    std::cout << " -x       Use exact arithmetic oracles instead of double precision.\n";
#endif /* IPO_WITH_RATIONAL */
    std::cout << std::flush;

    return EXIT_FAILURE;
  }

} /* namespace inverse */

