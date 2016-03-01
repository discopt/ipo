#ifndef IPO_PARSER_H_
#define IPO_PARSER_H_

#include <cstdlib>
#include <istream>
#include <vector>
#include <map>
#include <string>

#include "spx_gmp.h"

namespace ipo {

  struct LPToken
  {
    /// Special token types. All others are represented by themselves.
    
    static const char NAME = 'x';
    static const char ARROW = 'a';
    static const char DOUBLE_COLON = 'd';
    static const char NUMBER = 'n';
    static const char END_OF_FILE = (char)0;
    
    char type;
    std::size_t nameIndex;
    soplex::Rational number;
    
    LPToken(char type = 'Z');
    LPToken(std::size_t nameIndex);
    LPToken(const soplex::Rational& number);
    ~LPToken();
    
    inline
    bool operator()() const
    {
      return type != END_OF_FILE;
    }
  };
  
  class LPParser
  {
  public:
    LPParser(std::istream& stream);
    virtual ~LPParser();

  protected:
    LPToken next();
    LPToken nextNoWhite();

    std::ostream& printToken(std::ostream& stream, const LPToken& token);
    const std::string& getTokenName(const LPToken& token);

  private:
  public:
    LPToken extractName();
    LPToken extractNumber();
    
    std::string _buffer;
    int _lookAhead;
    std::vector<std::string> _nameList;
    std::map<std::string, std::size_t> _nameMap;
    std::istream& _stream;
  };

  class LPObjectiveParser : public LPParser
  {
  public:
    LPObjectiveParser(std::istream& stream);
    virtual ~LPObjectiveParser();
    
    void run();
    
  protected:
    virtual void handleObjective(const std::string& name, const std::map<std::string, soplex::Rational>& values) = 0;
    
    void parseObjective(bool maximize);
    void parseGoal();
    
  private:
    LPToken _token;
    std::map<std::string, soplex::Rational> _values;
  };
  
}

#endif