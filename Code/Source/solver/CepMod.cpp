// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepMod.h"

#include "ComMod.h"
#include "FE/Common/FEException.h"
#include "Parameters.h"
#include "utils.h"

#include <limits>
#include <math.h>

const std::map<ElectrophysiologyModelType, std::string> cep_model_type_to_name{
    {ElectrophysiologyModelType::NA, "NA"},
    {ElectrophysiologyModelType::AP, "AP"},
    {ElectrophysiologyModelType::BO, "BO"},
    {ElectrophysiologyModelType::FN, "FN"},
    {ElectrophysiologyModelType::TTP, "TTP"}
};

const std::map<std::string,ElectrophysiologyModelType> cep_model_name_to_type
{
  {"aliev-panfilov", ElectrophysiologyModelType::AP},
  {"ap", ElectrophysiologyModelType::AP},
  {"bueno-orovio", ElectrophysiologyModelType::BO},
  {"bo", ElectrophysiologyModelType::BO},
  {"fitzhugh-nagumo", ElectrophysiologyModelType::FN},
  {"fn", ElectrophysiologyModelType::FN},
  {"tentusscher-panfilov", ElectrophysiologyModelType::TTP},
  {"ttp", ElectrophysiologyModelType::TTP}
};

bool stimType::is_active(const double time) const
{
  const double eps = std::numeric_limits<double>::epsilon();

  const double time_relative = std::fmod(time, cycle_length);
  return start_time - eps <= time_relative && time_relative <= start_time + duration + eps;
}

void stimType::SpatialBounds::set_box(const Vector<double>& min, const Vector<double>& max)
{
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      min.size() == 0, "Stimulus box bounds must have at least one coordinate.");
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      min.size() != max.size(), "Stimulus box minimum and maximum must have the same coordinate dimension.");
  box_min = min;
  box_max = max;
  has_box = true;
}

void stimType::SpatialBounds::set_sphere(const Vector<double>& center, const double radius)
{
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      center.size() == 0, "Stimulus sphere center must have at least one coordinate.");
  svmp::throw_if<svmp::FE::InvalidArgumentException>(
      radius < 0.0, "Stimulus sphere radius must be non-negative.");
  sphere_center = center;
  sphere_radius = radius;
  has_sphere = true;
}

bool stimType::SpatialBounds::inside_box(const Vector<double>& x) const
{
  if (x.size() != box_min.size()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        "Point dimension does not match stimulus box dimension.");
  }

  for (int i = 0; i < x.size(); i++) {
    if (x[i] < box_min[i] || x[i] > box_max[i]) {
      return false;
    }
  }

  return true;
}

bool stimType::SpatialBounds::inside_sphere(const Vector<double>& x) const
{
  if (x.size() != sphere_center.size()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        "Point dimension does not match stimulus sphere dimension.");
  }

  double distance_squared = 0.0;

  for (int i = 0; i < x.size(); i++) {
    const double dx = x[i] - sphere_center[i];
    distance_squared += dx * dx;
  }

  return distance_squared <= sphere_radius * sphere_radius;
}

bool stimType::SpatialBounds::contains(const Vector<double>& x) const
{
  if (has_box && !inside_box(x)) return false;
  if (has_sphere && !inside_sphere(x)) return false;
  return true;
}

void stimType::SpatialBounds::distribute(const CmMod& cm_mod, const cmType& cm)
{
  cm.bcast(cm_mod, &has_box);
  cm.bcast(cm_mod, &has_sphere);

  if (has_box) {
    int box_size = box_min.size();
    cm.bcast(cm_mod, &box_size);

    if (cm.slv(cm_mod)) {
      box_min.resize(box_size);
      box_max.resize(box_size);
    }

    cm.bcast(cm_mod, box_min, "SpatialBounds box_min");
    cm.bcast(cm_mod, box_max, "SpatialBounds box_max");
  }

  if (has_sphere) {
    int sphere_center_size = sphere_center.size();
    cm.bcast(cm_mod, &sphere_center_size);

    if (cm.slv(cm_mod)) {
      sphere_center.resize(sphere_center_size);
    }

    cm.bcast(cm_mod, sphere_center, "SpatialBounds sphere_center");
    cm.bcast(cm_mod, &sphere_radius);
  }
}

double stimType::operator()(const double time, const Vector<double>& x) const
{
  if (utils::is_zero(amplitude) ||
      !is_active(time) ||
      !spatial_bounds.contains(x)) {
    return 0.0;
  }

  return amplitude;
}

void stimType::distribute(const CmMod& cm_mod, const cmType& cm)
{
  cm.bcast(cm_mod, &start_time);
  cm.bcast(cm_mod, &duration);
  cm.bcast(cm_mod, &cycle_length);
  cm.bcast(cm_mod, &amplitude);
  spatial_bounds.distribute(cm_mod, cm);
}

void stimType::read_parameters(const StimulusParameters& params, const int nsd, const double default_cycle_length)
{
  amplitude = params.amplitude.value();

  if (!utils::is_zero(amplitude)) {
    start_time = params.start_time.value();
    duration = params.duration.value();
    cycle_length = params.cycle_length.defined() ? params.cycle_length.value() : default_cycle_length;
  }

  const auto& spatial_bounds_params = params.spatial_bounds;

  const auto& box_params = spatial_bounds_params.box;
  const bool box_is_defined = box_params.defined();
  const bool box_min_defined = box_params.minimum.defined();
  const bool box_max_defined = box_params.maximum.defined();

  if (box_is_defined && !(box_min_defined && box_max_defined)) {
    svmp::raise<svmp::ParseException>(
        "Both Minimum and Maximum must be specified for a CEP stimulus box.");
  }

  if (box_is_defined) {
    const auto& box_min_vals = box_params.minimum.value();
    const auto& box_max_vals = box_params.maximum.value();

    if (box_min_vals.size() != box_max_vals.size()) {
      svmp::raise<svmp::ParseException>(
          "Stimulus box Minimum and Maximum must have the same coordinate dimension.");
    }

    if (box_min_vals.size() != static_cast<std::size_t>(nsd)) {
      svmp::raise<svmp::ParseException>(
          "Stimulus box coordinate dimension must match the simulation spatial dimension.");
    }

    Vector<double> box_min(box_min_vals.size());
    Vector<double> box_max(box_max_vals.size());

    for (int i = 0; i < static_cast<int>(box_min_vals.size()); i++) {
      if (box_min_vals[i] > box_max_vals[i]) {
        svmp::raise<svmp::ParseException>(
            "Stimulus box Minimum values must be less than or equal to Maximum values.");
      }

      box_min[i] = box_min_vals[i];
      box_max[i] = box_max_vals[i];
    }

    spatial_bounds.set_box(box_min, box_max);
  }

  const auto& sphere_params = spatial_bounds_params.sphere;
  const bool sphere_is_defined = sphere_params.defined();
  const bool sphere_center_defined = sphere_params.center.defined();
  const bool sphere_radius_defined = sphere_params.radius.defined();

  if (sphere_is_defined && !(sphere_center_defined && sphere_radius_defined)) {
    svmp::raise<svmp::ParseException>(
        "Both Center and Radius must be specified for a CEP stimulus sphere.");
  }

  if (sphere_is_defined) {
    const auto& sphere_center_vals = sphere_params.center.value();
    const double sphere_radius_val = sphere_params.radius.value();

    if (sphere_center_vals.size() != static_cast<std::size_t>(nsd)) {
      svmp::raise<svmp::ParseException>(
          "Stimulus sphere center coordinate dimension must match the simulation spatial dimension.");
    }

    if (sphere_radius_val < 0.0) {
      svmp::raise<svmp::ParseException>(
          "Stimulus sphere Radius must be non-negative.");
    }

    Vector<double> sphere_center(sphere_center_vals.size());

    for (int i = 0; i < static_cast<int>(sphere_center_vals.size()); i++) {
      sphere_center[i] = sphere_center_vals[i];
    }

    spatial_bounds.set_sphere(sphere_center, sphere_radius_val);
  }
}

cepModelType::cepModelType()
{
}

cepModelType::~cepModelType()
{
}

