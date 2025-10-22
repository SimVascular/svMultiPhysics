// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TXT_H 
#define TXT_H 

#include "Simulation.h"
#include "Array.h"
#include "ComMod.h"

#include "consts.h"

namespace txt_ns {

void create_boundary_integral_file(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const std::string& file_name);

void create_volume_integral_file(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const std::string& file_name);

void txt(Simulation* simulation, const bool flag);

void write_boundary_integral_data(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const int m, 
    const std::string file_name, const Array<double>& tmpV, const bool div, const bool pFlag);

void write_volume_integral_data(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const int m, 
    const std::string file_name, const Array<double>& tmpV, const bool div, const bool pFlag);

};

#endif

