// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepMod.h"

#include "ComMod.h"
#include "FE/Common/FEException.h"

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

bool stimType::inside_box(const ComMod& com_mod, const int Ac) const
{
  if (box_min.size() < com_mod.nsd || box_max.size() < com_mod.nsd) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        SVMP_HERE, "Stimulus box dimension is smaller than the simulation spatial dimension.");
  }

  for (int i = 0; i < com_mod.nsd; i++) {
    const double xi = com_mod.x(i, Ac);

    if (xi < box_min(i) || xi > box_max(i)) {
      return false;
    }
  }

  return true;
}

bool stimType::inside_sphere(const ComMod& com_mod, const int Ac) const
{
  if (sphere_center.size() < com_mod.nsd) {
    svmp::raise<svmp::FE::InvalidArgumentException>(
        SVMP_HERE, "Stimulus sphere center dimension is smaller than the simulation spatial dimension.");
  }

  double distance_squared = 0.0;

  for (int i = 0; i < com_mod.nsd; i++) {
    const double dx = com_mod.x(i, Ac) - sphere_center(i);
    distance_squared += dx * dx;
  }

  return distance_squared <= sphere_radius * sphere_radius;
}

bool stimType::contains_node(const ComMod& com_mod, const int Ac) const
{
  if (box_defined && !inside_box(com_mod, Ac)) {
    return false;
  }

  if (sphere_defined && !inside_sphere(com_mod, Ac)) {
    return false;
  }

  return true;
}

double stimType::operator()(const double time, const ComMod& com_mod, const int Ac) const
{
  if (!is_active(time)) {
    return 0.0;
  }

  if (!contains_node(com_mod, Ac)) {
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
  cm.bcast(cm_mod, &box_defined);

  if (box_defined) {
    int box_size = box_min.size();
    cm.bcast(cm_mod, &box_size);

    if (box_size <= 0) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          SVMP_HERE, "Stimulus box has invalid coordinate dimension.");
    }

    if (cm.slv(cm_mod)) {
      box_min.resize(box_size);
      box_max.resize(box_size);
    }

    cm.bcast(cm_mod, box_min, "Stimulus box_min");
    cm.bcast(cm_mod, box_max, "Stimulus box_max");
  }

  cm.bcast(cm_mod, &sphere_defined);

  if (sphere_defined) {
    int sphere_center_size = sphere_center.size();
    cm.bcast(cm_mod, &sphere_center_size);

    if (sphere_center_size <= 0) {
      svmp::raise<svmp::FE::InvalidArgumentException>(
          SVMP_HERE, "Stimulus sphere center has invalid coordinate dimension.");
    }

    if (cm.slv(cm_mod)) {
      sphere_center.resize(sphere_center_size);
    }

    cm.bcast(cm_mod, sphere_center, "Stimulus sphere_center");
    cm.bcast(cm_mod, &sphere_radius);
  }
}

cepModelType::cepModelType()
{
}

cepModelType::~cepModelType()
{
}

