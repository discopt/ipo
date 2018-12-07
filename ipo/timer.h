#ifndef IPO_TIMER_H_
#define IPO_TIMER_H_

#include "common.h"

#include <soplex/timer.h>

namespace ipo {

  class Timer: public soplex::Timer
  {
  public:

    /**
     * \brief Default constructor.
     *
     * Default constructor.
     */

    Timer();

    /**
     * \brief Copy constructor.
     *
     * Copy constructor.
     */

    Timer(const Timer& other);

    /**
     * \brief Destructor.
     *
     * Destructor.
     */

    virtual ~Timer();

    /**
     * \brief Assignment operator.
     */

    Timer& operator=(const Timer& other);

    /**
     * \brief Resets time to zero.
     *
     * Resets time to zero.
     */

    virtual void reset()
    {
      status = soplex::Timer::RESET;
      _userAccountedTicks = 0;
      _lastTime = 0.0;
    }

    /**
     * \brief Starts the timer.
     *
     * Starts the timer, resuming accounted user, system and real time.
     */

    virtual void start();

    /**
     * \brief Stops the timer and returns the accounted user time.
     *
     * Stops the timer and returns the accounted user time.
     */

    virtual soplex::Real stop();

    /**
     * \brief Returns the type of the timer.
     *
     * Returns the type of the timer.
     */

    virtual soplex::Timer::TYPE type()
    {
      return USER_TIME;
    }

    /**
     * \brief Adds \p additionalTime to the accounted user time.
     *
     * Adds \p additionalTime to the accounted user time.
     */

    void addTime(double additionalTime);

    /**
     * \brief Returns the accounted user time.
     *
     * Returns the accounted user time.
     */

    virtual double time() const;

    /**
     * \brief Returns the time from the last query.
     *
     * Returns the time from the last query.
     */

    virtual double lastTime() const;

  private:
    /**
     * \brief Converts ticks to seconds.
     *
     * Converts ticks to seconds.
     */

    double ticksToSeconds(long ticks) const
    {
      return (double(ticks) * 1000.0 / double(ticksPerSecond)) / 1000.0;
    }

    /**
     * \brief Updates the ticks from the system.
     *
     * Updates the ticks from the system.
     */

    void updateTicks();

  protected:
    static const long ticksPerSecond;

    long _userAccountedTicks;
    long _userTicks;
    soplex::Real _lastTime;
  };

  /**
   * \brief Adds the given number of seconds to all running timers.
   *
   * Adds the given number of seconds to all running timers.
   */

  void addTimeToRunningTimers(double time);

} // namespace soplex

#endif // IPO_TIMER_H_
