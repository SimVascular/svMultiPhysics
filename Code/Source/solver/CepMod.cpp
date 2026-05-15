// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "CepMod.h"
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

cepModelType::cepModelType()
{
}

cepModelType::~cepModelType()
{
}

