// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#ifndef ACTIVE_STRESS_TEST_HELPERS_H
#define ACTIVE_STRESS_TEST_HELPERS_H

#include "active_stress.h"
#include "Vector.h"

#include <cstddef>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// Generic active-stress trajectory runner.
//
// All calls to the node-local virtual interface go through the ActiveStress&
// base reference. This is required because derived-class overrides of
// init_local, advance_time_step_local, and compute_active_tension_local may
// remain protected — the public visibility declared in ActiveStress guarantees
// access only when the static type at the call site is ActiveStress.
//
// The caller must:
//   1. Construct the concrete model and configure its parameters.
//   2. Allocate `state` to `ActiveStress::n_states` elements.
//   3. Call model.init_local(state) before handing off to this function.
//   4. Pass per-step calcium, fiber_stretch, and fiber_stretch_rate via the
//      three callable arguments.
//   5. Provide an `on_step` callback that is invoked *after* each advance
//      with the step index, updated state, and computed active tension.
// ---------------------------------------------------------------------------

inline void run_active_stress_trajectory(
    ActiveStress &model,
    Vector<double> &state,
    int n_steps,
    double dt,
    std::function<double(int)> ca_at,
    std::function<double(int)> lam_at,
    std::function<double(int)> dlam_at,
    std::function<void(int, const Vector<double> &, double)> on_step)
{
  for (int step = 0; step < n_steps; ++step) {
    const double t     = step * dt;
    const double ca    = ca_at(step);
    const double lam   = lam_at(step);
    const double dlam  = dlam_at(step);

    // All calls through ActiveStress& — required for correct access control.
    model.advance_time_step_local(t, dt, ca, lam, dlam, state);

    const double Ta = model.compute_active_tension_local(state, lam);
    on_step(step, state, Ta);
  }
}

// ---------------------------------------------------------------------------
// CSV loader.
//
// Reads a CSV file with an optional leading comment block (lines starting with
// '#'), a mandatory header row, and one data row per checkpoint. Returns the
// rows as a vector of double vectors. The header row is discarded; callers
// index columns by their position in the header (0 = first non-comment column).
//
// Throws std::runtime_error on file-open failure or parse errors.
// ---------------------------------------------------------------------------

inline std::vector<std::vector<double>>
load_csv(const std::string &path)
{
  std::ifstream f(path);
  if (!f.is_open())
    throw std::runtime_error("Cannot open reference CSV: " + path);

  std::vector<std::vector<double>> rows;
  std::string line;
  bool header_seen = false;

  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    if (!header_seen) {
      header_seen = true;   // first non-comment line is the column header
      continue;
    }

    std::vector<double> row;
    std::istringstream ss(line);
    std::string token;
    while (std::getline(ss, token, ',')) {
      size_t pos;
      row.push_back(std::stod(token, &pos));
    }
    rows.push_back(std::move(row));
  }

  return rows;
}

#endif // ACTIVE_STRESS_TEST_HELPERS_H
