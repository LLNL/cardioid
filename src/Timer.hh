////////////////////////////////////////////////////////////////////////////////
//
//  Timer.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id$

#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <sys/time.h>

class Timer
{
  private:

  clock_t clk;
  double t,total_cpu,total_real;
  int running_;

  public:

  Timer() : total_cpu(0.0), total_real(0.0), running_(0) {};

  void reset() { total_cpu = 0.0; total_real = 0.0; running_ = 0; };

  void start()
  {
    clk = clock();
    t = gtod();
    running_ = 1;
  };

  int running() { return running_; };

  void stop()
  {
    if ( running_ ) 
    {
      total_cpu += ((double)(clock()-clk))/CLOCKS_PER_SEC;
      total_real += gtod()-t;
      running_ = 0;
    }
  };

  double cpu()
  {
    if ( running_ ) 
    {
      return total_cpu + ((double)(clock()-clk))/CLOCKS_PER_SEC;
    }
    else
    {
      return total_cpu;
    } 
  };

  double real()
  {
    if ( running_ ) 
    {
      return total_real + gtod()-t;
    }
    else
    {
      return total_real;
    } 
  };

  double gtod(void)
  {
    static struct timeval tv;
    static struct timezone tz;
    gettimeofday(&tv,&tz);
    return tv.tv_sec + 1.e-6*tv.tv_usec;
  }
};
#endif
