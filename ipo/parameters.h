#ifndef IPO_PARAMETERS_H_
#define IPO_PARAMETERS_H_

#include <string>

#include "common.h"

namespace ipo {

  class ParameterGroup;

  class Parameter
  {
  public:
    Parameter(ParameterGroup* group);
    virtual ~Parameter();

  protected:
    ParameterGroup* _group;
  };

  template <typename T>
  class TypeParameter : public Parameter
  {
  public:
    TypeParameter(ParameterGroup* group, const std::string& name, const std::string& description, const T& defaultValue,
      T& valueStorage)
      : Parameter(group), _name(name), _description(description), _defaultValue(defaultValue), _valueStorage(valueStorage)
    {

    }

    virtual ~TypeParameter()
    {

    }

    virtual bool set(const std::string& valueString) = 0;
    virtual void get(std::ostream& stream) = 0;

  protected:
    const std::string _name;
    const std::string _description;
    T _defaultValue;
    T& _valueStorage;
  };

  class IntParameter : public TypeParameter<int>
  {
  public:
    IntParameter(ParameterGroup* group, const std::string& name, const std::string& description, const int& defaultValue,
      int& valueStorage);
    virtual ~IntParameter();
    virtual bool set(const std::string& valueString);
    virtual void get(std::ostream& stream);
  };

  class DoubleParameter : public TypeParameter<double>
  {
  public:
    DoubleParameter(ParameterGroup* group, const std::string& name, const std::string& description, const double& defaultValue,
      double& valueStorage);
    virtual ~DoubleParameter();
    virtual bool set(const std::string& valueString);
    virtual void get(std::ostream& stream);
  };

  class CharParameter : public TypeParameter<char>
  {
  public:
    CharParameter(ParameterGroup* group, const std::string& name, const std::string& description, const char& defaultValue,
      char& valueStorage);
    virtual ~CharParameter();
    virtual bool set(const std::string& valueString);
    virtual void get(std::ostream& stream);
  };

  class BoolParameter : public TypeParameter<bool>
  {
  public:
    BoolParameter(ParameterGroup* group, const std::string& name, const std::string& description,
      const bool& defaultValue, bool& valueStorage);
    virtual ~BoolParameter();
    virtual bool set(const std::string& valueString);
    virtual void get(std::ostream& stream);
  };

  class StringParameter : public TypeParameter<std::string>
  {
  public:
    StringParameter(ParameterGroup* group, const std::string& name, const std::string& description,
      const std::string& defaultValue, std::string& valueStorage);
    virtual ~StringParameter();
    virtual bool set(const std::string& valueString);
    virtual void get(std::ostream& stream);
  };

} /* namespace ipo */

#endif /* IPO_PARAMETERS_H_ */
