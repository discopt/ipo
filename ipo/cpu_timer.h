#ifndef IPO_CPU_TIMER_H_
#define IPO_CPU_TIMER_H_

#include <soplex.h>
#include <spxdefines.h>
#include <timer.h>

namespace ipo {

  class CPUTimer: public soplex::Timer
  {
  private:

    //------------------------------------
    /**@name number of ticks per second */
    //@{
    static const long ticks_per_sec; ///< ticks per secound, should be constant
    //@}

    //------------------------------------
    /**@name Data */
    //@{
    mutable long uAccount; ///< user time
    mutable long uTicks; ///< user ticks

    mutable soplex::Real lasttime;
    //@}

    //------------------------------------
    /**@name Internal helpers */
    //@{
    /// convert ticks to secounds.
    soplex::Real ticks2sec(long ticks) const
    {
      return (soplex::Real(ticks) * 1000.0 / soplex::Real(ticks_per_sec)) / 1000.0;
    }

    /// get actual user ticks from the system.
    void updateTicks() const;

    //@}

  public:

    //------------------------------------
    /**@name Construction / destruction */
    //@{
    /// default constructor
    CPUTimer();

    /// copy constructor
    CPUTimer(const CPUTimer& old);

    /// assignment operator
    CPUTimer& operator=(const CPUTimer& old);

    virtual ~CPUTimer();
    //@}

    //------------------------------------
    /**@name Control */
    //@{
    /// initialize timer, set timing accounts to zero.
    virtual void reset()
    {
      status = soplex::Timer::RESET;
      uAccount = 0;
      lasttime = 0.0;
    }

    /// start timer, resume accounting user, system and real time.
    virtual void start();

    /// stop timer, return accounted user time.
    virtual soplex::Real stop();

    /// return type of timer
    virtual soplex::Timer::TYPE type()
    {
      return USER_TIME;
    }

    void addAccount(double additionalTime);
    //@}

    //------------------------------------
    /**@name Access */
    //@{
    virtual soplex::Real time() const;

    virtual soplex::Real lastTime() const;
    //@}
  };

  void addTimeToActiveTimers(double time);

} // namespace soplex
#endif // CPU_TIMER_H_
