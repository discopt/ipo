#include "parser.h"

#include <cassert>

#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>

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
  
  LPToken LPParser::nextNoWhite()
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

  LPObjectiveParser::LPObjectiveParser(std::istream& stream): LPParser(stream)
  {
    
  }

  LPObjectiveParser::~LPObjectiveParser()
  {
    
  }
  
  
  void LPObjectiveParser::parseObjective(bool maximize)//, const std::string* name)
  {
//     std::cout << "parseObjective" << std::endl;
    std::string name = "";
    bool first = true;
    bool onlyName = true;
    while (true)
    {
//       printToken(std::cout, _token) << std::endl;
      
      if (_token.type == '\n')
        _token = nextNoWhite();
      bool negate = false;
      if (_token.type == '+' || _token.type == '-')
      {
        negate = _token.type == '-';
        _token = nextNoWhite();
        onlyName = false;
      }

      soplex::Rational value = maximize ? 1 : -1;
      if (_token.type == LPToken::NUMBER)
      {
        value *= _token.number;
        _token = nextNoWhite();
        if (_token.type == '/')
        {
          _token = nextNoWhite();
          if (_token.type == LPToken::NUMBER)
          {
            value /= _token.number;
            _token = nextNoWhite();
          }
          else
            return;
        }
        else if (_token.type == '*')
          _token = nextNoWhite();
        onlyName = false;
      }
      if (negate)
        value = -value;

      if (_token.type == LPToken::NAME)
      {
        const std::string& var = getTokenName(_token);
        std::map<std::string, soplex::Rational>::iterator iter = _values.find(var);
        if (iter != _values.end())
          iter->second += value;
        else
          _values.insert(std::make_pair(var, value));
        _token = nextNoWhite();
        if (first && onlyName)
        {
          assert(_values.size() == 1);
          if (_token.type == ':')
          {
            _values.clear();
            name = var;
            _token = nextNoWhite();
          }
        }
        first = false;
      }
      else
      {
        if (!first)
          handleObjective(name, _values);
        return;
      }
    }
  }

  void LPObjectiveParser::parseGoal()
  {
//     std::cout << "parseGoal: ";
//     printToken(std::cout, _token) << std::endl;

    if (_token.type == LPToken::NAME)
    {
      std::string lowered =  getTokenName(_token);
      std::transform(lowered.begin(), lowered.end(), lowered.begin(), ::tolower);
      _token = nextNoWhite();
      if (lowered == "maximize" || lowered == "maximum" || lowered == "max")
        parseObjective(true);
      if (lowered == "minimize" || lowered == "minimum" || lowered == "min")
        parseObjective(false);
    }
    else
      _token = nextNoWhite();
  }

  void LPObjectiveParser::run()
  {
    _token = nextNoWhite();
    while (true)
    {
      if (_token.type == LPToken::END_OF_FILE)
        break;
      parseGoal();
    }
  }


}