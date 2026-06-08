// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TIMER_H 
#define TIMER_H 

#include <sys/time.h>

/// @brief Keep track of time
class Timer 
{
  public:

    double get_elapsed_time() const
    {
      return get_time() - current_time;
    }

    double get_time() const
    {
      timeval now{};
      gettimeofday(&now, nullptr);
      return static_cast<double>(now.tv_sec) +
             static_cast<double>(now.tv_usec) * 1.0e-6;
    }

    void set_time()
    {
      current_time = get_time();
    }

    double current_time{0.0};
};

#endif
