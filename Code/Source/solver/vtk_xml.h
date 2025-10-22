// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "Simulation.h"
#include "ComMod.h"
#include "Array.h"

#ifndef VTK_XML_H
#define VTK_XML_H

namespace vtk_xml {

void do_test();

void int_msh_data(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, dataType& d, const int outDof, const int nOute);

void read_vtp(const std::string& file_name, faceType& face);

void read_vtp_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, const int idx, faceType& face);

void read_vtu(const std::string& file_name, mshType& mesh);

void read_precomputed_solution_vtu(const std::string& file_name, const std::string& field_name, mshType& mesh);

void read_vtu_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, const int idx, mshType& mesh);

void read_vtus(Simulation* simulation, Array<double>& lA, Array<double>& lY, Array<double>& lD, const std::string& fName);

void write_vtp(ComMod& com_mod, faceType& lFa, const std::string& fName);

void write_vtu(ComMod& com_mod, mshType& lM, const std::string& fName);

void write_vtu_debug(ComMod& com_mod, mshType& lM, const std::string& fName);

void write_vtus(Simulation* simulation, const Array<double>& lA, const Array<double>& lY, const Array<double>& lD, const bool lAve);

};

#endif


