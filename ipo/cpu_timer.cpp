#include "cpu_timer.h"

#include <cassert>
#include <set>
#include <sys/times.h>

using namespace soplex;

namespace ipo {

  std::set<CPUTimer*> allTimers;

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

  const long CPUTimer::ticks_per_sec = long(TIMES_TICKS_PER_SEC);

  CPUTimer::CPUTimer() :
      Timer(), uAccount(0), uTicks(0), lasttime(0.0)
  {
    assert(ticks_per_sec > 0);
  }

  /// copy constructor
  CPUTimer::CPUTimer(const CPUTimer& old) :
      Timer(), uAccount(old.uAccount), uTicks(old.uTicks), lasttime(old.lasttime)
  {
    assert(ticks_per_sec > 0);
  }

  /// assignment operator
  CPUTimer& CPUTimer::operator=(const CPUTimer& old)
  {
    assert(ticks_per_sec > 0);
    uAccount = old.uAccount;
    uTicks = old.uTicks;
    lasttime = old.lasttime;
    return *this;
  }

  CPUTimer::~CPUTimer()
  {
    if (status == RUNNING)
      allTimers.erase(this);
  }

// get actual user, system and real time from system
  void CPUTimer::updateTicks() const
  {
#if defined(_WIN32) || defined(_WIN64)

    uTicks = clock();

#else   /* !(_WIN32 || _WIN64) */

    struct tms now;
    clock_t ret = times(&now);

    if (int(ret) == -1)
      now.tms_utime = now.tms_stime = ret = 0;

    uTicks = long(now.tms_utime);

#endif  /* !(_WIN32 || _WIN64) */
  }

// start timer, resume accounting user, system and real time.
  void CPUTimer::start()
  {
    // ignore start request if timer is running
    if (status != RUNNING)
    {
      updateTicks();

      uAccount -= uTicks;
      status = RUNNING;
      allTimers.insert(this);
    }
    lasttime = 0;
  }

// stop timer, return accounted user time.
  Real CPUTimer::stop()
  {
    // status remains unchanged if timer is not running
    if (status == RUNNING)
    {
      updateTicks();

      uAccount += uTicks;
      status = STOPPED;
      allTimers.erase(this);
    }
    return ticks2sec(uAccount);
  }

  void CPUTimer::addAccount(double additionalTime)
  {
    uAccount += Real(ticks_per_sec) * additionalTime;
  }

// get accounted user time.
  Real CPUTimer::time() const
  {
    if (status == RUNNING)
    {
      updateTicks();
      lasttime = ticks2sec(uTicks + uAccount);
    }
    else
    {
      lasttime = ticks2sec(uAccount);
    }
    return lasttime;
  }

  Real CPUTimer::lastTime() const
  {
    return lasttime;
  }

  void addTimeToActiveTimers(double time)
  {
    for (std::set<CPUTimer*>::iterator iter = allTimers.begin(); iter != allTimers.end(); ++iter)
    {
      (*iter)->addAccount(time);
    }
  }

} // namespace soplex
