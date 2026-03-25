// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// This function will return the current time in sec. 

#include <chrono>

namespace fsi_linear_solver {

double fsils_cpu_t()
{
  //auto now = std::chrono::system_clock::now();
  //return static_cast<double>(std::chrono::system_clock::to_time_t(now));

  auto now = std::chrono::system_clock::now();
  auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);

  auto value = now_ms.time_since_epoch();
  auto duration = value.count() / 1000.0;
  return static_cast<double>(duration);
}

};
