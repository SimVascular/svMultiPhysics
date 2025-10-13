// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepMod.h"
#include <math.h>

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

const std::map<std::string,TimeIntegratioType> cep_time_int_to_type 
{
  {"cn", TimeIntegratioType::CN2},
  {"cn2", TimeIntegratioType::CN2},
  {"implicit", TimeIntegratioType::CN2},

  {"fe", TimeIntegratioType::FE},
  {"euler", TimeIntegratioType::FE},
  {"explicit", TimeIntegratioType::FE},

  {"rk", TimeIntegratioType::RK4},
  {"rk4", TimeIntegratioType::RK4},
  {"runge", TimeIntegratioType::RK4}

};


cepModelType::cepModelType()
{
}

cepModelType::~cepModelType()
{
}

