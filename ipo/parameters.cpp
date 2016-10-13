#include "parameters.h"

#include <sstream>

namespace ipo {

  Parameter::Parameter(ParameterGroup* group)
    : _group(group)
  {

  }

  Parameter::~Parameter()
  {

  }


  IntParameter::IntParameter(ParameterGroup* group, const std::string& name, const std::string& description,
    const int& defaultValue, int& valueStorage)
    : TypeParameter<int>(group, name, description, defaultValue, valueStorage)
  {

  }

  IntParameter::~IntParameter()
  {

  }

  bool IntParameter::set(const std::string& valueString)
  {
    std::stringstream stream(valueString);
    stream >> _valueStorage;
    return stream.good();
  }

  void IntParameter::get(std::ostream& stream)
  {
    stream << _valueStorage;
  }

  DoubleParameter::DoubleParameter(ParameterGroup* group, const std::string& name, const std::string& description,
    const double& defaultValue, double& valueStorage)
    : TypeParameter<double>(group, name, description, defaultValue, valueStorage)
  {

  }

  DoubleParameter::~DoubleParameter()
  {

  }

  bool DoubleParameter::set(const std::string& valueString)
  {
    std::stringstream stream(valueString);
    stream >> _valueStorage;
    return stream.good();
  }

  void DoubleParameter::get(std::ostream& stream)
  {
    stream << _valueStorage;
  }

  CharParameter::CharParameter(ParameterGroup* group, const std::string& name, const std::string& description,
    const char& defaultValue, char& valueStorage)
    : TypeParameter<char>(group, name, description, defaultValue, valueStorage)
  {

  }

  CharParameter::~CharParameter()
  {

  }

  bool CharParameter::set(const std::string& valueString)
  {
    std::stringstream stream(valueString);
    stream >> _valueStorage;
    return stream.good();
  }

  void CharParameter::get(std::ostream& stream)
  {

  }

  BoolParameter::BoolParameter(ParameterGroup* group, const std::string& name, const std::string& description,
    const bool& defaultValue, bool& valueStorage)
    : TypeParameter<bool>(group, name, description, defaultValue, valueStorage)
  {

  }

  BoolParameter::~BoolParameter()
  {

  }

  bool BoolParameter::set(const std::string& valueString)
  {
    if (valueString == "true" || valueString == "yes" || valueString == "y" || valueString == "1" || valueString == "on")
    {
      _valueStorage = true;
      return true;
    }
    else if (valueString == "false" || valueString == "no" || valueString == "n" || valueString == "0" || valueString == "off")
    {
      _valueStorage = false;
      return false;
    }
    else
      return false;
  }

  void BoolParameter::get(std::ostream& stream)
  {
    stream << (_valueStorage ? "true" : "false");
  }

  StringParameter::StringParameter(ParameterGroup* group, const std::string& name, const std::string& description,
    const std::string& defaultValue, std::string& valueStorage)
    : TypeParameter<std::string>(group, name, description, defaultValue, valueStorage)
  {

  }

  StringParameter::~StringParameter()
  {

  }

  bool StringParameter::set(const std::string& valueString)
  {
    _valueStorage = valueString;
  }

  void StringParameter::get(std::ostream& stream)
  {
    stream << _valueStorage;
  }



} /* namespace ipo */
