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

  const int icl = static_cast<int>(fmax(floor(time / CL), 0.0));
  const double Ts_cycle = Ts + static_cast<double>(icl) * CL;
  const double Te_cycle = Ts_cycle + Td;

  return Ts_cycle - eps <= time && time <= Te_cycle + eps;
}

void SpatialBounds::set_box(const Vector<double>& min, const Vector<double>& max)
{
  box_min_ = min;
  box_max_ = max;
  region_type_ = (region_type_ == RegionType::sphere) ? RegionType::both : RegionType::box;
}

void SpatialBounds::set_sphere(const Vector<double>& center, const double radius)
{
  sphere_center_ = center;
  sphere_radius_ = radius;
  region_type_ = (region_type_ == RegionType::box) ? RegionType::both : RegionType::sphere;
}

bool SpatialBounds::inside_box(const Vector<double>& x) const
{
  if (x.size() != box_min_.size()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        SVMP_HERE, "Point dimension does not match stimulus box dimension.");
  }

  for (int i = 0; i < x.size(); i++) {
    if (x(i) < box_min_(i) || x(i) > box_max_(i)) {
      return false;
    }
  }

  return true;
}

bool SpatialBounds::inside_sphere(const Vector<double>& x) const
{
  if (x.size() != sphere_center_.size()) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        SVMP_HERE, "Point dimension does not match stimulus sphere dimension.");
  }

  double distance_squared = 0.0;

  for (int i = 0; i < x.size(); i++) {
    const double dx = x(i) - sphere_center_(i);
    distance_squared += dx * dx;
  }

  return distance_squared <= sphere_radius_ * sphere_radius_;
}

bool SpatialBounds::contains(const Vector<double>& x) const
{
  switch (region_type_) {
    case RegionType::none:   return true;
    case RegionType::box:    return inside_box(x);
    case RegionType::sphere: return inside_sphere(x);
    case RegionType::both:   return inside_box(x) && inside_sphere(x);
  }
  return true;
}

void SpatialBounds::distribute(const CmMod& cm_mod, const cmType& cm)
{
  int region_type_int = static_cast<int>(region_type_);
  cm.bcast(cm_mod, &region_type_int);
  region_type_ = static_cast<RegionType>(region_type_int);

  if (region_type_ == RegionType::box || region_type_ == RegionType::both) {
    int box_size = box_min_.size();
    cm.bcast(cm_mod, &box_size);

    if (box_size <= 0) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          SVMP_HERE, "Stimulus box has invalid coordinate dimension.");
    }

    if (cm.slv(cm_mod)) {
      box_min_.resize(box_size);
      box_max_.resize(box_size);
    }

    cm.bcast(cm_mod, box_min_, "SpatialBounds box_min");
    cm.bcast(cm_mod, box_max_, "SpatialBounds box_max");
  }

  if (region_type_ == RegionType::sphere || region_type_ == RegionType::both) {
    int sphere_center_size = sphere_center_.size();
    cm.bcast(cm_mod, &sphere_center_size);

    if (sphere_center_size <= 0) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          SVMP_HERE, "Stimulus sphere center has invalid coordinate dimension.");
    }

    if (cm.slv(cm_mod)) {
      sphere_center_.resize(sphere_center_size);
    }

    cm.bcast(cm_mod, sphere_center_, "SpatialBounds sphere_center");
    cm.bcast(cm_mod, &sphere_radius_);
  }
}

double stimType::operator()(const double time, const Vector<double>& x) const
{
  if (utils::is_zero(A)) {
    return 0.0;
  }

  if (!is_active(time)) {
    return 0.0;
  }

  if (!spatial_bounds_.contains(x)) {
    return 0.0;
  }

  return A;
}

void stimType::distribute(const CmMod& cm_mod, const cmType& cm)
{
  cm.bcast(cm_mod, &Ts);
  cm.bcast(cm_mod, &Td);
  cm.bcast(cm_mod, &CL);
  cm.bcast(cm_mod, &A);
  spatial_bounds_.distribute(cm_mod, cm);
}

void stimType::read_parameters(const StimulusParameters& params, const int nsd, const double default_cycle_length)
{
  A = params.amplitude.value();

  if (!utils::is_zero(A)) {
    Ts = params.start_time.value();
    Td = params.duration.value();
    CL = params.cycle_length.defined() ? params.cycle_length.value() : default_cycle_length;
  }

  const auto& spatial_bounds_params = params.spatial_bounds;

  const auto& box_params = spatial_bounds_params.box;
  const bool box_is_defined = box_params.defined();
  const bool box_min_defined = box_params.minimum.defined();
  const bool box_max_defined = box_params.maximum.defined();

  if (box_is_defined && !(box_min_defined && box_max_defined)) {
    svmp::raise<svmp::ParseException>(
        SVMP_HERE, "Both Minimum and Maximum must be specified for a CEP stimulus box.");
  }

  if (box_is_defined) {
    const auto& box_min_vals = box_params.minimum.value();
    const auto& box_max_vals = box_params.maximum.value();

    if (box_min_vals.size() != box_max_vals.size()) {
      svmp::raise<svmp::ParseException>(
          SVMP_HERE, "Stimulus box Minimum and Maximum must have the same coordinate dimension.");
    }

    if (box_min_vals.size() < static_cast<std::size_t>(nsd)) {
      svmp::raise<svmp::ParseException>(
          SVMP_HERE, "Stimulus box dimension is smaller than the simulation spatial dimension.");
    }

    Vector<double> box_min(box_min_vals.size());
    Vector<double> box_max(box_max_vals.size());

    for (int i = 0; i < static_cast<int>(box_min_vals.size()); i++) {
      if (box_min_vals[i] > box_max_vals[i]) {
        svmp::raise<svmp::ParseException>(
            SVMP_HERE, "Stimulus box Minimum values must be less than or equal to Maximum values.");
      }

      box_min(i) = box_min_vals[i];
      box_max(i) = box_max_vals[i];
    }

    spatial_bounds_.set_box(box_min, box_max);
  }

  const auto& sphere_params = spatial_bounds_params.sphere;
  const bool sphere_is_defined = sphere_params.defined();
  const bool sphere_center_defined = sphere_params.center.defined();
  const bool sphere_radius_defined = sphere_params.radius.defined();

  if (sphere_is_defined && !(sphere_center_defined && sphere_radius_defined)) {
    svmp::raise<svmp::ParseException>(
        SVMP_HERE, "Both Center and Radius must be specified for a CEP stimulus sphere.");
  }

  if (sphere_is_defined) {
    const auto& sphere_center_vals = sphere_params.center.value();
    const double sphere_radius_val = sphere_params.radius.value();

    if (sphere_center_vals.size() < static_cast<std::size_t>(nsd)) {
      svmp::raise<svmp::ParseException>(
          SVMP_HERE, "Stimulus sphere center dimension is smaller than the simulation spatial dimension.");
    }

    if (sphere_radius_val < 0.0) {
      svmp::raise<svmp::ParseException>(
          SVMP_HERE, "Stimulus sphere Radius must be non-negative.");
    }

    Vector<double> sphere_center(sphere_center_vals.size());

    for (int i = 0; i < static_cast<int>(sphere_center_vals.size()); i++) {
      sphere_center(i) = sphere_center_vals[i];
    }

    spatial_bounds_.set_sphere(sphere_center, sphere_radius_val);
  }
}

cepModelType::cepModelType()
{
}

cepModelType::~cepModelType()
{
}

