#ifndef VIENNASHE_UTIL_TIMER_HPP
#define VIENNASHE_UTIL_TIMER_HPP
/* ============================================================================
   Copyright (c) 2011-2022, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
     ViennaSHE - The Vienna Spherical Harmonics Expansion Boltzmann Solver
                            -----------------

                    http://viennashe.sourceforge.net/

   License:         MIT (X11), see file LICENSE in the base directory
=============================================================================== */

/** @file viennashe/util/timer.hpp
    @brief Contains a timer class. Contains one windows implementation and one *NIX implementation
 */


#include <iostream>

#ifdef _WIN32

  #include <windows.h>
  #undef min
  #undef max

#else

  #include <sys/time.h>

#endif

namespace viennashe
{
  namespace util
  {

#ifdef _WIN32 // WINDOWS


    /** @brief A simple timer class */
    class timer
    {
    public:

      timer()  { QueryPerformanceFrequency(&freq_); }

      /** @brief Starts the timer */
      void start()
      {
        QueryPerformanceCounter((LARGE_INTEGER*) &start_time_);
      }

      /** @brief Returns the number of seconds elapsed */
      double get() const
      {
        LARGE_INTEGER  end_time;
        QueryPerformanceCounter((LARGE_INTEGER*) &end_time);
        return (static_cast<double>(end_time.QuadPart) - static_cast<double>(start_time_.QuadPart)) / static_cast<double>(freq_.QuadPart);
      }

    private:
      LARGE_INTEGER freq_;
      LARGE_INTEGER start_time_;
    };

#else // NOT WINDOWS

    /** @brief A simple timer class */
    class timer
    {
    public:

      timer() : ts_(0) { }

      /** @brief Starts the timer */
      void start()
      {
        struct timeval tval;
        gettimeofday(&tval, 0);
        ts_ = tval.tv_sec * 1000000 + tval.tv_usec;
      }

      /** @brief Returns the number of seconds elapsed */
      double get() const
      {
        struct timeval tval;
        gettimeofday(&tval, 0);

        long end_time = tval.tv_sec * 1000000 + tval.tv_usec;

        return static_cast<double>(end_time - ts_) / 1000000.0;
      }

    private:
      long ts_;
    };


#endif

  } // namespace util
} //  namespace viennashe

#endif
