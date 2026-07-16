// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TIMER_H 
#define TIMER_H 

#include <chrono>

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
      const auto now = std::chrono::steady_clock::now();
      return std::chrono::duration<double>(now.time_since_epoch()).count();
    }

    void set_time()
    {
      current_time = get_time();
    }

    double current_time{0.0};
};

#endif
