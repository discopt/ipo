#ifndef IPO_PARSER_H_
#define IPO_PARSER_H_

#include <cstdlib>
#include <istream>
#include <vector>
#include <map>
#include <string>

#include "ipo.h"

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
    LPToken nextNonWhite();
    void fetchNextNonWhite();

    inline const LPToken& token() const
    {
      return _token;
    }

    std::ostream& printToken(std::ostream& stream, const LPToken& token);
    const std::string& getTokenName(const LPToken& token);

  private:
    LPToken extractName();
    LPToken extractNumber();

    std::string _buffer;
    int _lookAhead;
    std::vector<std::string> _nameList;
    std::map<std::string, std::size_t> _nameMap;
    std::istream& _stream;
    LPToken _token;
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
  };

  class LPInequalityParser : public LPParser
  {
  public:
    LPInequalityParser(std::istream& stream);
    virtual ~LPInequalityParser();

    void run();

  protected:
    void parseRhs(const std::string& name, const soplex::Rational& lhsValue, char lhsSign, const std::map<std::string, soplex::Rational>& coefficients);
    void parseVector(const std::string& name, const soplex::Rational& lhsValue, char lhsSign, bool triedLhs, const soplex::Rational& coefficient, const std::string& triedVar);
    void parseLhs(const std::string& name);
    void parseName();

    virtual void handleInequality(const std::string& name, const soplex::Rational& lhs, const std::map<std::string, soplex::Rational>& values, const soplex::Rational& rhs) = 0;
  };

  class PointParser : public LPParser
  {
  public:
    PointParser(std::istream& stream);
    virtual ~PointParser();

    void run();

  protected:
    void parseVector(const std::string& name);
    void parseName();

    virtual void handlePoint(const std::string& name, const std::map<std::string, soplex::Rational>& values) = 0;
  };
}

#endif
