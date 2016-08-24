#include "timer.h"

#include <cassert>
#include <set>
#include <sys/times.h>

using namespace soplex;

namespace ipo {

  std::set<Timer*> allTimers;

  /* determine TIMES_TICKS_PER_SEC for clock ticks delivered by times().
   * (don't use CLOCKS_PER_SEC since this is related to clock() only).
   */
#if defined(CLK_TCK)
#define TIMES_TICKS_PER_SEC CLK_TCK
#elif defined(_SC_CLK_TCK)
#define TIMES_TICKS_PER_SEC sysconf(_SC_CLK_TCK)
#elif defined(HZ)
#define TIMES_TICKS_PER_SEC HZ
#else // !CLK_TCK && !_SC_CLK_TCK && !HZ
#define TIMES_TICKS_PER_SEC 60
#endif // !CLK_TCK && !_SC_CLK_TCK && !HZ

  const long Timer::ticksPerSecond = long(TIMES_TICKS_PER_SEC);

  Timer::Timer() :
      soplex::Timer(), _userAccountedTicks(0), _userTicks(0), _lastTime(0.0)
  {
    assert(ticksPerSecond > 0);
  }

  Timer::Timer(const Timer& other) 
    : soplex::Timer(), _userAccountedTicks(other._userAccountedTicks), _userTicks(other._userTicks), _lastTime(other._lastTime)
  {
    assert(ticksPerSecond > 0);
  }

  Timer::~Timer()
  {

  }

  Timer& Timer::operator=(const Timer& other)
  {
    assert(ticksPerSecond > 0);
    _userAccountedTicks = other._userAccountedTicks;
    _userTicks = other._userTicks;
    _lastTime = other._lastTime;
    return *this;
  }

  void Timer::updateTicks()
  {
#if defined(_WIN32) || defined(_WIN64)

    _userTicks = clock();

#else   /* !(_WIN32 || _WIN64) */

    struct tms now;
    clock_t ret = times(&now);

    if (int(ret) == -1)
      now.tms_utime = now.tms_stime = ret = 0;

    _userTicks = long(now.tms_utime);

#endif  /* !(_WIN32 || _WIN64) */
  }

  void Timer::start()
  {
    if (status != RUNNING)
    {
      updateTicks();

      _userAccountedTicks -= _userTicks;
      status = RUNNING;
      allTimers.insert(this);
    }
    _lastTime = 0;
  }

  double Timer::stop()
  {
    if (status == RUNNING)
    {
      updateTicks();

      _userAccountedTicks += _userTicks;
      status = STOPPED;
      allTimers.erase(this);
    }
    return ticksToSeconds(_userAccountedTicks);
  }

  void Timer::addTime(double additionalTime)
  {
    _userAccountedTicks += double(ticksPerSecond) * additionalTime;
  }

  double Timer::time() const
  {
    if (status == RUNNING)
    {
      const_cast<Timer*>(this)->updateTicks();
      const_cast<Timer*>(this)->_lastTime = ticksToSeconds(_userTicks + _userAccountedTicks);
    }
    else
    {
      const_cast<Timer*>(this)->_lastTime = ticksToSeconds(_userAccountedTicks);
    }
    return _lastTime;
  }

  double Timer::lastTime() const
  {
    return _lastTime;
  }

  void addTimeToRunningTimers(double time)
  {
    for (std::set<Timer*>::iterator iter = allTimers.begin(); iter != allTimers.end(); ++iter)
    {
      (*iter)->addTime(time);
    }
  }

} // namespace soplex
