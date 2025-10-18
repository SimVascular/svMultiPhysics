// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TIMER_H 
#define TIMER_H 

#include <chrono>
#include <iostream>
#include <string>

/// @brief Keep track of time
class Timer 
{
  public:

    double get_elapsed_time()
    {
      return get_time() - current_time;
    }

    double get_time()
    {
      auto now = std::chrono::system_clock::now();
      auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);

      auto value = now_ms.time_since_epoch();
      auto duration = value.count() / 1000.0;
      return static_cast<double>(duration);
    }

    void set_time()
    {
      current_time = get_time();
    }

    double current_time;
};

#endif

