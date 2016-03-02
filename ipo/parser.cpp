#include "parser.h"

#include <cassert>

#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>

using namespace soplex;

namespace ipo {

  LPToken::LPToken(char t) : type(t), nameIndex(0), number(0.0)
  {
    
  }

  LPToken::LPToken(std::size_t nameIdx) : type(NAME), nameIndex(nameIdx), number(0.0)
  {
    
  }

  LPToken::LPToken(const soplex::Rational& num) : type(NUMBER), nameIndex(0), number(num)
  {
    
  }

  LPToken::~LPToken()
  {
    
  }
  
  LPParser::LPParser(std::istream& stream) : _stream(stream), _lookAhead(EOF)
  {

  }
  
  LPParser::~LPParser()
  {
    
  }
  
  std::ostream& LPParser::printToken(std::ostream& stream, const LPToken& token)
  {
    if (token.type == LPToken::NAME)
      stream << "[" << getTokenName(token) << "]";
    else if (token.type == LPToken::NUMBER)
      stream << "[" << token.number << "]";
    else if (token.type == LPToken::ARROW)
      stream << "[->]";
    else if (token.type == LPToken::DOUBLE_COLON)
      stream << "[::]";
    else if (token.type == LPToken::END_OF_FILE)
      stream << "[EOF]";
    else
      stream << "[" << token.type << "]";
    return stream;
  }

  const std::string& LPParser::getTokenName(const LPToken& token)
  {
    assert(token.nameIndex < _nameList.size());
    return _nameList[token.nameIndex];
  }

  LPToken LPParser::extractName()
  {
    std::map<std::string, std::size_t>::const_iterator nameIter = _nameMap.find(_buffer);
    if (nameIter != _nameMap.end())
      return LPToken(nameIter->second);
    _nameMap[_buffer] = _nameList.size();
    _nameList.push_back(_buffer);
    return LPToken(_nameList.size() - 1);
  }
  
  LPToken LPParser::extractNumber()
  {
    soplex::Rational q;
    if (!soplex::readStringRational(_buffer.c_str(), q))
    {
      _buffer = "Failed to parse \"" + _buffer + "\" as a floating-point number.";
      throw std::runtime_error(_buffer);
    }
    else
      return LPToken(q);
  }

  LPToken LPParser::next()
  {
    static const char STATE_INITIAL = 0;
    static const char STATE_COLON = 1;
    static const char STATE_MINUS = 2;
    static const char STATE_LESS = 3;
    static const char STATE_GREATER = 4;
    static const char STATE_EQUAL = 5;
    static const char STATE_WHITESPACE = 6;
    static const char STATE_COMMENT = 7;
    static const char STATE_NAME = 8;
    static const char STATE_NUMBER_INTEGER = 9;
    static const char STATE_NUMBER_FRACTIONAL = 10;
    static const char STATE_NUMBER_EXPONENT_SIGN = 11;
    static const char STATE_NUMBER_EXPONENT_NUMBER = 12;

    char state = STATE_INITIAL;
    _buffer.clear();
    while (_lookAhead != EOF || !_stream.fail())
    {
      int c = _lookAhead != EOF ? _lookAhead : _stream.get();
      _lookAhead = EOF;

      switch (state)
      {
        case STATE_INITIAL:
          switch (c)
          {
            case ':':
              state = STATE_COLON;
              break;
            case '-':
              state = STATE_MINUS;
              break;
            case '<':
              state = STATE_LESS;
              break;
            case '>':
              state = STATE_GREATER;
              break;
            case '=':
              state = STATE_EQUAL;
              break;
            case ' ':
            case '\t':
              state = STATE_WHITESPACE;
              break;
            case '\\':
              state = STATE_COMMENT;
              break;
            case '\n':
            case '+': 
            case '^':
            case '/':
            case '[':
            case ']':
            case '*':
              return LPToken(char(c));
            default:
            {
              if (c == '.' || (c >= '0' && c <= '9'))
              {
                state = STATE_NUMBER_INTEGER;
                _buffer = c;
                break;
              }
              else if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '!' || c == '"' || c == '#' || c == '$'
                || c == '%' || c == '&' || c == ';' || c == '?' || c == '@' || c == '_' || c == '\'' || c == '{'
                || c == '}' || c == '~')
              {
                state = STATE_NAME;
                _buffer = c;
                break;
              }
              else if (c == EOF)
                return LPToken(LPToken::END_OF_FILE);
              std::stringstream stream;
              stream << "Parser error: Unexpected character " << c << " = '" << char(c) << "'.";
              throw std::runtime_error(stream.str());
            }
          }
        break;
        case STATE_COLON:
          if (c == ':')
            return LPToken(LPToken::DOUBLE_COLON);
          _lookAhead = c;
          return LPToken(':');
        case STATE_MINUS:
          if (c == '>')
            return LPToken(LPToken::ARROW);
          _lookAhead = c;
          return LPToken('-');
        case STATE_LESS:
          if (c != '=')
            _lookAhead = c;
          return LPToken('<');
        case STATE_GREATER:
          if (c != '=')
            _lookAhead = c;
          return LPToken('>');
        case STATE_EQUAL:
          if (c != '=')
            _lookAhead = c;
          return LPToken('=');
        case STATE_WHITESPACE:
          if (c != ' ' && c != '\t')
          {
            _lookAhead = c;
            return LPToken(' ');
          }
          break;
        case STATE_COMMENT:
          if (c == '\n')
          {
            state = STATE_INITIAL;
            break;
          }
          else if (c == EOF)
            return LPToken(LPToken::END_OF_FILE);
          break;
        case STATE_NAME:
        {
          if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '!' || c == '"' || c == '#' || c == '$'
            || c == '%' || c == '&' || c == ';' || c == '?' || c == '@' || c == '_' || c == '\'' || c == '{' 
            || c == '}' || c == '~' || (c >= '0' && c <= '9') || c == '.')
          {
            _buffer += c;
            break;
          }
          else
          {
            state = STATE_INITIAL;
            _lookAhead = c;
            return extractName();
          }
        }
        case STATE_NUMBER_INTEGER:
        {
          if (c >= '0' && c <= '9')
          {
            _buffer += c;
            break;
          }
          else if (c == '.')
          {
            _buffer += c;
            state = STATE_NUMBER_FRACTIONAL;
            break;
          }
          else if (c == 'e')
          {
            _buffer += c;
            state = STATE_NUMBER_EXPONENT_SIGN;
            break;
          }
          else
          {
            state = STATE_INITIAL;
            _lookAhead = c;
            return extractNumber();
          }
        }
        case STATE_NUMBER_FRACTIONAL:
        {
          if (c >= '0' && c <= '9')
          {
            _buffer += c;
            break;
          }
          else if (c == 'e')
          {
            _buffer += c;
            state = STATE_NUMBER_EXPONENT_SIGN;
            break;
          }
          else
          {
            state = STATE_INITIAL;
            _lookAhead = c;
            return extractNumber();
          }
        }
        case STATE_NUMBER_EXPONENT_SIGN:
        {
          if (c == '+' || c == '-' || (c >= '0' && c <= '9'))
          {
            _buffer += c;
            state = STATE_NUMBER_EXPONENT_NUMBER;
            break;
          }
          else
          {
            state = STATE_INITIAL;
            _lookAhead = c;
            return extractNumber();
          }
        }
        case STATE_NUMBER_EXPONENT_NUMBER:
        {
          if (c >= '0' && c <= '9')
          {
            _buffer += c;
            break;
          }
          else
          {
            state = STATE_INITIAL;
            _lookAhead = c;
            return extractNumber();
          }
        }
        default:
          throw std::runtime_error("Unknown parser state.");
      }
    }
    return LPToken(LPToken::END_OF_FILE);
  }
  
  LPToken LPParser::nextNonWhite()
  {
    LPToken token;
    do
    {
      token = next();
//       printToken(std::cout, token) << std::endl;
    }
    while(token.type == ' ');
    return token;
  }
  
  void LPParser::fetchNextNonWhite()
  {
    _token = nextNonWhite();
  }

  LPObjectiveParser::LPObjectiveParser(std::istream& stream): LPParser(stream)
  {
    
  }

  LPObjectiveParser::~LPObjectiveParser()
  {
    
  }
  
  
  void LPObjectiveParser::parseObjective(bool maximize)
  {
    std::string name = "";
    bool first = true;
    bool onlyName = true;
    std::map<std::string, soplex::Rational> coefficients;
    while (true)
    {
      if (token().type == '\n')
        fetchNextNonWhite();

      bool negate = false;
      if (token().type == '+' || token().type == '-')
      {
        negate = token().type == '-';
        fetchNextNonWhite();
        onlyName = false;
      }
      else
      {
        if (!first)
          return handleObjective(name, coefficients);
      }

      soplex::Rational value = maximize ? 1 : -1;
      if (token().type == LPToken::NUMBER)
      {
        value *= token().number;
        fetchNextNonWhite();
        if (token().type == '/')
        {
          fetchNextNonWhite();
          if (token().type == LPToken::NUMBER)
          {
            value /= token().number;
            fetchNextNonWhite();
          }
          else
            return;
        }
        else if (token().type == '*')
          fetchNextNonWhite();
        onlyName = false;
      }
      if (negate)
        value = -value;

      if (token().type == LPToken::NAME)
      {
        const std::string& var = getTokenName(token());
        std::map<std::string, soplex::Rational>::iterator iter = coefficients.find(var);
        if (iter != coefficients.end())
          iter->second += value;
        else
          coefficients.insert(std::make_pair(var, value));
        fetchNextNonWhite();
        if (first && onlyName)
        {
          assert(coefficients.size() == 1);
          if (token().type == ':')
          {
            coefficients.clear();
            name = var;
            fetchNextNonWhite();
          }
        }
        first = false;
      }
      else
      {
        if (!first)
          handleObjective(name, coefficients);
        return;
      }
    }
  }

  void LPObjectiveParser::parseGoal()
  {
    if (token().type == LPToken::NAME)
    {
      std::string lowered =  getTokenName(token());
      std::transform(lowered.begin(), lowered.end(), lowered.begin(), ::tolower);
      fetchNextNonWhite();
      if (lowered == "maximize" || lowered == "maximum" || lowered == "max")
        parseObjective(true);
      if (lowered == "minimize" || lowered == "minimum" || lowered == "min")
        parseObjective(false);
    }
    else
      fetchNextNonWhite();
  }

  void LPObjectiveParser::run()
  {
    fetchNextNonWhite();
    while (token().type != LPToken::END_OF_FILE)
    {
      parseGoal();
    }
  }
  
  LPInequalityParser::LPInequalityParser(std::istream& stream): LPParser(stream)
  {

  }

  LPInequalityParser::~LPInequalityParser()
  {

  }
  
  void LPInequalityParser::parseRhs(const std::string& name, const Rational& lhsValue, char lhsSign, const std::map< std::string, Rational>& coefficients)
  {
    char rhsSign = token().type;
    assert(rhsSign == '<' || rhsSign == '>' || rhsSign == '=');
    fetchNextNonWhite();
    
//     std::cout << "rhsSign = " << rhsSign << ", token = ";
//     printToken(std::cout, token()) << std::endl;
    
    bool negate = false;
    if (token().type == '+' || token().type == '-')
    {
      negate = token().type == '-';
      fetchNextNonWhite();
    }
    Rational rhsValue = negate ? -1 : 1;
    if (token().type == LPToken::NUMBER)
    {
      rhsValue *= token().number;
      fetchNextNonWhite();
      if (token().type == '/')
      {
        fetchNextNonWhite();
        if (token().type == LPToken::NUMBER)
        {
          rhsValue /= token().number;
          fetchNextNonWhite();
        }
        else
        {
          rhsSign = '<';
          rhsValue = infinity;
        }
      }
    }
    else
    {
      rhsSign = '<';
      rhsValue = infinity;
    }
    
    Rational left = -infinity;
    Rational right = infinity;
    if (lhsSign == '<' || lhsSign == '=')
      left = std::max(left, lhsValue);
    if (lhsSign == '>' || lhsSign == '=')
      right = std::min(right, lhsValue);
    if (rhsSign == '<' || rhsSign == '=')
      right = std::min(right, rhsValue);
    if (rhsSign == '>' || rhsSign == '=')
      left = std::max(left, rhsValue);
    return handleInequality(name, left, coefficients, right);
  }

  void LPInequalityParser::parseVector(const std::string& name, const Rational& lhsValue, char lhsSign, bool triedLhs, const Rational& coefficient, const std::string& triedVar)
  {
//     std::cout << "parseVector(" << name << ", " << lhsValue << " " << lhsSign << ", " << (triedLhs ? "true" : "false") << ", " << coefficient << ", " << triedVar << "), token = ";
//     printToken(std::cout, token()) << std::endl;
    
    std::string varName = triedVar;
    Rational value = coefficient;
    bool negate;
    std::map<std::string, soplex::Rational> values;
    bool first = true;
    while (token().type != LPToken::END_OF_FILE)
    {
      if (triedLhs)
      {
        assert(token().type == LPToken::NAME);
        varName = getTokenName(token());
        triedLhs = false;
        fetchNextNonWhite();
      }
      else if (varName == "")
      {
        if (token().type == '\n')
          fetchNextNonWhite();
        
        negate = false;
        if (token().type == '+' || token().type == '-')
        {
          negate = token().type == '-';
          fetchNextNonWhite();
        }
        else
        {
          if (!first)
          {
            if (token().type == '<' || token().type == '>' || token().type == '=')
            {
              return parseRhs(name, lhsValue, lhsSign, values);
            }
            else if (lhsSign != '<' || lhsValue > -infinity)
            {
              if (lhsSign == '<')
                return handleInequality(name, lhsValue, values, infinity);
              else if (lhsSign == '>')
                return handleInequality(name, -infinity, values, lhsValue);
              else
                return handleInequality(name, lhsValue, values, lhsValue);
            }
            else if (token().type == LPToken::NAME)
            {
//               std::cout << "Found a name ";
//               printToken(std::cout, token()) << " although +/- expected." << std::endl;
              return parseName();
            }
            else
              return;
          }
        }
        
        value = negate ? -1 : 1;
        if (token().type == LPToken::NUMBER)
        {
          value *= token().number;
          fetchNextNonWhite();
          if (token().type == '/')
          {
            fetchNextNonWhite();
            if (token().type == LPToken::NUMBER)
            {
              value /= token().number;
              fetchNextNonWhite();
            }
            else
              return;
          }
        }

        if (token().type == '*')
          fetchNextNonWhite();
        
        if (token().type == LPToken::NAME)
        {
          varName = getTokenName(token());
          fetchNextNonWhite();
        }
        else
          return;
      }

      std::map<std::string, Rational>::iterator iter = values.find(varName);
      if (iter != values.end())
        iter->second += value;
      else
        values.insert(std::make_pair(varName, value));
      
//       std::cout << "Added a coefficient with var = " << varName << ", next token = ";
//       printToken(std::cout, token()) << std::endl;

      varName = "";
      first = false;
    }
  }

  void LPInequalityParser::parseLhs(const std::string& name)
  {
//     std::cout << "parseLhs(" << name << ")" << std::endl;
    
    bool negate = false;
    if (token().type == '+' || token().type == '-')
    {
      negate = token().type == '-';
      fetchNextNonWhite();
    }
    Rational coefficient = negate ? -1 : 1;
    if (token().type == LPToken::NUMBER)
    {
      coefficient *= token().number;
      fetchNextNonWhite();
      if (token().type == '/')
      {
        fetchNextNonWhite();
        if (token().type == LPToken::NUMBER)
        {
          coefficient /= token().number;
          fetchNextNonWhite();
        }
        else
          return;
      }
    }
    else if (token().type == LPToken::NAME)
    {
      std::string lowered = getTokenName(token());
      std::transform(lowered.begin(), lowered.end(), lowered.begin(), ::tolower);
      if (lowered == "infinity" || lowered == "inf")
      {
        coefficient *= infinity;
        fetchNextNonWhite();
      }
      else
        return parseVector(name, -infinity, '<', true, coefficient, "");
    }
    else
      return;

    /// We now have read an optional sign and a number/infinity.
    
    if (token().type == '<' || token().type == '>' || token().type == '=')
    {
      char t = token().type;
      fetchNextNonWhite();
      return parseVector(name, coefficient, t, false, 0, "");
    }
    else if (token().type == '*')
      fetchNextNonWhite();
    
    if (token().type == LPToken::NAME)
      return parseVector(name, -infinity, '<', true, coefficient, "");
    else
      return;
  }

      
  void LPInequalityParser::parseName()
  {
    assert(token().type == LPToken::NAME);
    const std::string& name = getTokenName(token());
    fetchNextNonWhite();
    if (token().type == ':')
    {
      fetchNextNonWhite();
      if (token().type == '\n')
        fetchNextNonWhite();
      parseLhs(name);
    }
    else
    {
      parseVector("", -soplex::infinity, '<', false, 1, name);
    }
  }


  void LPInequalityParser::run()
  {
    fetchNextNonWhite();
    while (token().type != LPToken::END_OF_FILE)
    {
//       std::cout << "run: ";
//       printToken(std::cout, token()) << std::endl;
      
      if (token().type == '+' || token().type == '-' || token().type == LPToken::NUMBER)
        parseLhs("");
      else if (token().type == LPToken::NAME)
        parseName();
      else
        fetchNextNonWhite();
    }
  }


}