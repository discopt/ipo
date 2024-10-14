#pragma once

#include <iostream>
#include <sstream>
#include <random>
#include <unordered_map>

#include <rapidxml.hpp>
#include <rapidxml_utils.hpp>
#include <rapidxml_print.hpp>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <ipo/oracles.hpp>
#include <ipo/affine_hull.hpp>
#include <ipo/lp.hpp>

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
    std::cout << "Solving inverse MIP <" << instance.name << ">." << std::endl;

    std::vector<std::shared_ptr<ipo::Polyhedron<Number>>> polyhedra;
    for (std::size_t i = 0; i < instance.solvers.size(); ++i)
    {
      auto oracle = instance.solvers[i]->template getOptimizationOracle<Number>();
      polyhedra.push_back( std::make_shared<ipo::Polyhedron<Number>>(oracle));
    }

    for (std::size_t i = 0; i < instance.solvers.size(); ++i)
    {
      // TODO: For debugging purposes we start with the affine hull.

      ipo::AffineHull<Number> affineHull;
      std::cerr << "Starting affine hull computation.\n" << std::flush;
      ipo::AffineHullQuery affQuery;
      affineHull = ipo::affineHull(polyhedra[i], affQuery);
      std::cout << "Dimension: " << affineHull.dimension << " / " << polyhedra[i]->space()->dimension() << std::endl;
    }
  }

  // template <typename Number>
  // void run(std::shared_ptr<ipo::OptimizationOracle<Number>> oracle)
  // {
  //   auto poly = std::make_shared<ipo::Polyhedron<Number>>(oracle);
  //
  //   ipo::AffineHull<Number> affineHull;
  //   std::cerr << "Starting affine hull computation.\n" << std::flush;
  //   ipo::AffineHullQuery affQuery;
  //   // affQuery.timeLimit = timeLimit;
  //   affineHull = ipo::affineHull(poly, affQuery);
  //   std::cout << "Dimension: " << affineHull.dimension << " / " << poly->space()->dimension() << std::endl;
  //
  //   assert(false);
  // }

  int printUsage(const std::string& program)
  {
    std::cout << program << " [OPTIONS] FILE...\n";
    std::cout << "Solves inverse optimization problems on a polyhedron defined by FILE.\n";
    std::cout << "General options:\n";
    std::cout << " -h       Show this help and exit.\n";
    std::cout << "Oracle/polyhedron options:\n";
#if defined(IPO_RATIONAL)
    std::cout << " -x       Use exact arithmetic oracles instead of double precision.\n";
#endif /* IPO_RATIONAL */
    std::cout << std::flush;

    return EXIT_FAILURE;
  }

} /* namespace inverse */

