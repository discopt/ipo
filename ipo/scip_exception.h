#ifndef SCIP_EXCEPTION_H_
#define SCIP_EXCEPTION_H_

#include "common.h"

#include <exception>
#include <string>

#ifdef IPO_WITH_SCIP

#ifdef NDEBUG
  #undef NDEBUG
  #include <scip/scip.h>
  #include <scip/pub_misc.h>
  #define NDEBUG
#else
  #include <scip/scip.h>
  #include <scip/pub_misc.h>
#endif

#define SCIP_MSG_MAX 256

/**
 * @brief translates a SCIP_RETCODE into an error string
 *
 * @param[in] retcode SCIP_RETCODE you want to translate
 * @param[out] buffer_str buffer to character array to store translated message, this must be at least of size \ref SCIP_MSG_MAX
 * @return buffer_str or NULL, if retcode could not be translated
 */
inline char* SCIPgetErrorString(SCIP_RETCODE retcode, char* buffer_str, int buffersize)
{
  // the following was copied from SCIPprintError
  switch (retcode)
  {
  case SCIP_OKAY:
    SCIPsnprintf(buffer_str, buffersize, "normal termination");
    return buffer_str;
  case SCIP_ERROR:
    SCIPsnprintf(buffer_str, buffersize, "unspecified error");
    return buffer_str;
  case SCIP_NOMEMORY:
    SCIPsnprintf(buffer_str, buffersize, "insufficient memory error");
    return buffer_str;
  case SCIP_READERROR:
    SCIPsnprintf(buffer_str, buffersize, "file read error");
    return buffer_str;
  case SCIP_WRITEERROR:
    SCIPsnprintf(buffer_str, buffersize, "file write error");
    return buffer_str;
  case SCIP_NOFILE:
    SCIPsnprintf(buffer_str, buffersize, "file not found error");
    return buffer_str;
  case SCIP_FILECREATEERROR:
    SCIPsnprintf(buffer_str, buffersize, "cannot create file");
    return buffer_str;
  case SCIP_LPERROR:
    SCIPsnprintf(buffer_str, buffersize, "error in LP solver");
    return buffer_str;
  case SCIP_NOPROBLEM:
    SCIPsnprintf(buffer_str, buffersize, "no problem exists");
    return buffer_str;
  case SCIP_INVALIDCALL:
    SCIPsnprintf(buffer_str, buffersize, "method cannot be called at this time in solution process");
    return buffer_str;
  case SCIP_INVALIDDATA:
    SCIPsnprintf(buffer_str, buffersize, "method cannot be called with this type of data");
    return buffer_str;
  case SCIP_INVALIDRESULT:
    SCIPsnprintf(buffer_str, buffersize, "method returned an invalid result code");
    return buffer_str;
  case SCIP_PLUGINNOTFOUND:
    SCIPsnprintf(buffer_str, buffersize, "a required plugin was not found");
    return buffer_str;
  case SCIP_PARAMETERUNKNOWN:
    SCIPsnprintf(buffer_str, buffersize, "the parameter with the given name was not found");
    return buffer_str;
  case SCIP_PARAMETERWRONGTYPE:
    SCIPsnprintf(buffer_str, buffersize, "the parameter is not of the expected type");
    return buffer_str;
  case SCIP_PARAMETERWRONGVAL:
    SCIPsnprintf(buffer_str, buffersize, "the value is invalid for the given parameter");
    return buffer_str;
  case SCIP_KEYALREADYEXISTING:
    SCIPsnprintf(buffer_str, buffersize, "the given key is already existing in table");
    return buffer_str;
  case SCIP_MAXDEPTHLEVEL:
    SCIPsnprintf(buffer_str, buffersize, "maximal branching depth level exceeded");
    return buffer_str;
  default:
    return NULL;
  }
}

/** @brief exception handling class for SCIP
 *
 * this class enables you to handle the return code of SCIP functions in a C++ way
 */
class SCIPException: public std::exception
{
private:
  char _msg[SCIP_MSG_MAX]; ///< error message
  SCIP_RETCODE _retcode; ///< SCIP error code

public:

  /** @brief constructs a SCIPEexception from an error code
   *
   * this constructs a new SCIPException from given error code
   * @param[in] retcode SCIP error code
   */
  SCIPException(SCIP_RETCODE retcode) : _retcode(retcode)
  {
    if (SCIPgetErrorString(retcode, _msg, SCIP_MSG_MAX) == NULL)
      SCIPsnprintf(_msg, SCIP_MSG_MAX, "unknown SCIP retcode %d", retcode);
  }

  /** @brief returns the error message
   *
   * @return error message
   *
   * this overrides the corresponding method of std::exception in order to allow you to catch all of your exceptions as std::exception
   */
  const char * what(void) const throw ()
  {
    return _msg;
  }

  /** @brief get method for @ref _retcode
   *
   * @return stored SCIP_RETCODE
   */
  SCIP_RETCODE getRetcode(void) const
  {
    return _retcode;
  }

  /** destructor */
  ~SCIPException(void) throw ()
  {
  }
};

/**
 * @brief macro to call scip function with exception handling
 *
 * this macro calls a SCIP function and throws an instance of SCIPException if anything went wrong
 *
 **/

#define SCIP_CALL_EXC(x)                  \
   {                                      \
      SCIP_RETCODE _retcode;			        \
      if ((_retcode = (x)) != SCIP_OKAY)	\
      {						                        \
	      throw SCIPException(_retcode);		\
      }						                        \
   }

#endif

#endif
