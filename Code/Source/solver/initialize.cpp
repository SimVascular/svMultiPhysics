/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// The code here replicates the Fortran code in DISTRIBUTE.f.

#include "initialize.h"

#include "distribute.h"

#include "all_fun.h"
#include "baf_ini.h"
#include "cep_ion.h"
#include "consts.h"
#include "fs.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "output.h"
#include "post.h"
#include "set_bc.h"
#include "txt.h"
#include "utils.h"
#include "vtk_xml.h"

#include "commu.h"
#include "fsils.hpp"
#include "lhs.h"

#include "CmMod.h"

//#include "set_equation_dof.h"

#include "mpi.h"

#include <fstream>
#include <iostream>
#include <math.h>

/// @brief This seems to deallocate a bunch of arrys.
///
/// \todo [TODO:DaveP] not implemented.
//
void finalize(Simulation* simulation)
{
}

/// @brief Using the svFSI specific format binary file for initialization
///
/// Reprodices 'SUBROUTINE INITFROMBIN(fName, timeP)' defined in INITIALIZE.f.
//
void init_from_bin(Simulation* simulation, const std::string& fName, std::array<double,3>& timeP)
{
  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  int task_id = cm.idcm();
  #ifdef debug_init_from_bin 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto const recLn = com_mod.recLn;
  auto& cm_mod = simulation->cm_mod;
  auto& cep_mod = simulation->cep_mod;
  auto& cem = cep_mod.cem;

  bool ibFlag = com_mod.ibFlag;
  bool dFlag = com_mod.dFlag;
  bool sstEq = com_mod.sstEq;
  bool pstEq = com_mod.pstEq;
  bool cepEq = cep_mod.cepEq;

  com_mod.timeP[0] = timeP[0];
  com_mod.timeP[1] = timeP[1];
  com_mod.timeP[2] = timeP[2];

  #ifdef debug_init_from_bin 
  dmsg << "dFlag: " << dFlag;
  dmsg << "sstEq: " << sstEq;
  dmsg << "pstEq: " << pstEq;
  dmsg << "cepEq: " << cepEq;
  dmsg << "cm.tF(): " << cm.tF(cm_mod); 
  #endif

  // Open file and position at the location in the file
  // for the current process.
  //
  std::ifstream bin_file(fName, std::ios::binary | std::ios::in);
  int process_id = cm.tF(cm_mod);
  std::streampos write_pos = (process_id  - 1) * recLn;
  bin_file.seekg(write_pos);
  //OPEN(fid, FILE=fName, ACCESS='DIRECT', RECL=recLn)

  std::array<int,7> tStamp;
  auto& cplBC = com_mod.cplBC;
  auto& Ad = com_mod.Ad;
  auto& Ao = com_mod.Ao;
  auto& Do = com_mod.Do;
  auto& pS0 = com_mod.pS0;
  auto& Xion = cep_mod.Xion;
  auto& Yo = com_mod.Yo;

  output::read_restart_header(com_mod, tStamp, timeP[0], bin_file);
  bin_file.read((char*)cplBC.xo.data(), cplBC.xo.msize());
  bin_file.read((char*)Yo.data(), Yo.msize());
  bin_file.read((char*)Ao.data(), Ao.msize());

  if (!ibFlag) {

    // Update mesh and Dn-Do variables
    //
    if (dFlag) {
      bin_file.read((char*)Do.data(), Do.msize());

      if (sstEq) {
        if (pstEq) {
          bin_file.read((char*)pS0.data(), pS0.msize());
          bin_file.read((char*)Ad.data(), Ad.msize());
        } else if (cepEq) {
          bin_file.read((char*)Ad.data(), Ad.msize());
          bin_file.read((char*)Xion.data(), Xion.msize());
          bin_file.read((char*)cem.Ya.data(), cem.Ya.msize());
        } else {
          bin_file.read((char*)Ad.data(), Ad.msize());
        }

      } else {
        if (pstEq) {
          bin_file.read((char*)pS0.data(), pS0.msize());
        } else if (cepEq) {
          bin_file.read((char*)Xion.data(), Xion.msize());
          bin_file.read((char*)cem.Ya.data(), cem.Ya.msize());
        } else {
          //READ(fid,REC=cm.tF()) tStamp, cTS, time, timeP(1), eq.iNorm, cplBC.xo, Yo, Ao, Do
        }
      }

    // Mesh and Dn-Do variables not updated.
    //
    } else {
      if (cepEq) {
        bin_file.read((char*)Xion.data(), Xion.msize());
      } else {
        //READ(fid,REC=cm.tF()) tStamp, cTS, time, timeP(1), eq.iNorm, cplBC.xo, Yo, Ao
      }
    }

  // ibFlag
  //
  // [TODO:DaveP] not implemented.
  //
  } else {
    if (dFlag) {
      if (pstEq) {
        //READ(fid,REC=cm.tF()) tStamp, cTS, time, timeP(1), eq.iNorm, cplBC.xo, Yo, Ao, Do, pS0, ib.Yb, ib.Auo, ib.Ubo
      } else {
         //READ(fid,REC=cm.tF()) tStamp, cTS, time, timeP(1), eq.iNorm, cplBC.xo, Yo, Ao, Do, ib.Yb, ib.Auo, ib.Ubo
      }
    } else {
      //READ(fid,REC=cm.tF()) tStamp, cTS, time, timeP(1), eq.iNorm, cplBC.xo, Yo, Ao, ib.Yb, ib.Auo, ib.Ubo
    }

    /*
    if (com_mod.ib.cpld == ibCpld_I) {
      com_mod.ib.Aun = com_mod.ib.Auo;
      com_mod.ib.Ubn = com_mod.ib.Ubo;
    }
    */
  }

  bin_file.close();

  // First checking all variables on master processor, since on the
  // other processor data will be shifted due to any change on the
  // sizes
  //
  if (cm.mas(cm_mod)) {
    auto stamp = com_mod.stamp;
    std::vector<std::string> error_msgs = {
      "The number of processors ",
      "The number of equations ",
      "The number of meshes ",
      "The number of nodes ",
      "The number of ncplBC.x ",
      "The number of dof ",
      "The dFlag specification "
    };

    std::string error_msg;

    for (int i = 0; i < stamp.size(); i++) {
      if (tStamp[i] != stamp[i]) {
        error_msg = error_msgs[i];
        error_msg += " " + std::to_string(tStamp[i]) + " ";
        break;
      }
    }

    if (error_msg != "") { 
      error_msg += "read from '" + fName + "' don't match.";
      throw std::runtime_error(error_msg);
    }
  }

  int i = 0;
  cm.bcast(cm_mod, &i);

  #ifdef debug_init_from_bin 
  dmsg << "cTS: " << com_mod.cTS; 
  dmsg << "time: " << com_mod.time; 
  dmsg << "timeP(1): " << timeP[0]; 
  #endif

  // if (ANY(tStamp != stamp)) err = "Simulation stamp"// " does not match with "//fName
}

/// @brief Using the saved VTU files for initialization
///
/// Reproduces the Fortran 'INITFROMVTU' subroutine.
//
void init_from_vtu(Simulation* simulation, const std::string& fName, std::array<double,3>& timeP)
{
  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;

  #ifdef debug_init_from_vtu
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  com_mod.cTS = 0;
  com_mod.time = 0.0;
  com_mod.timeP[0] = timeP[0];

  for (auto& eq : com_mod.eq) {
    eq.iNorm = 0.0;
  }

  const int gtnNo = com_mod.gtnNo;
  const int tDof = com_mod.tDof;

  Array<double> tmpA, tmpY, tmpD;

  if (cm.mas(cm_mod)) {
    tmpA.resize(tDof,gtnNo); 
    tmpY.resize(tDof,gtnNo); 
    tmpD.resize(tDof,gtnNo);
    vtk_xml::read_vtus(simulation, tmpA, tmpY, tmpD, fName);
  } else { 
    //ALLOCATE(tmpA(0,0), tmpY(0,0), tmpD(0,0))
  } 

  com_mod.Ao = all_fun::local(com_mod, cm_mod, cm, tmpA);
  com_mod.Yo = all_fun::local(com_mod, cm_mod, cm, tmpY);
  com_mod.Do = all_fun::local(com_mod, cm_mod, cm, tmpD);
}

/// @brief Initialize or finalize svFSI variables/structures.
///
/// \code {.cpp}
/// Modifies the following for each com_mod.eq[].
///   eq.am
///   eq.dof
///   eq.pNorm
///   eq.af
///   eq.beta
//   eq.gam
///   eq.s
///   eq.e
///
///   com_mod.Ao.resize(tDof,tnNo);
///   com_mod.An.resize(tDof,tnNo);
///   com_mod.Yo.resize(tDof,tnNo);
///   com_mod.Yn.resize(tDof,tnNo);
///   com_mod.Do.resize(tDof,tnNo);
///   com_mod.Dn.resize(tDof,tnNo);
///   com_mod.Bf.resize(nsd,tnNo);
///
///   com_mod.An = com_mod.Ao;
///   com_mod.Yn = com_mod.Yo;
///   com_mod.Dn = com_mod.Do;
/// \endcode
///
/// This function replicates the Fortran 'SUBROUTINE INITIALIZE(timeP)' in INITIALIZE.f.
//
void initialize(Simulation* simulation, Vector<double>& timeP)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& cep_mod = simulation->cep_mod; 
  auto& chnl_mod = simulation->chnl_mod;
  int nsd = com_mod.nsd;
  #include "set_equation_dof.h"

  #define n_debug_initialize
  #ifdef debug_initialize
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Setup logging to history file.
  if (!simulation->com_mod.cm.slv(simulation->cm_mod)) {
    std::string hist_file_name;

    if (chnl_mod.appPath != "") { 
      auto mkdir_arg = std::string("mkdir -p ") + chnl_mod.appPath;
      std::system(mkdir_arg.c_str());
      hist_file_name = chnl_mod.appPath + "/" + simulation->history_file_name;
    } else {
      hist_file_name = simulation->history_file_name;
    }

    bool output_to_cout = true;
    simulation->logger.initialize(hist_file_name, output_to_cout);
  }

  #ifdef debug_initialize
  dmsg << "Set time " << std::endl;
  dmsg << "com_mod.timeP[0]: " << com_mod.timeP[0] << std::endl;
  #endif
  com_mod.timeP[0] = timeP[0];
  com_mod.timeP[1] = timeP[1];
  com_mod.timeP[2] = timeP[2];

  // These are global variables that are set later in this function.
  #ifdef debug_initialize
  dmsg << "Get globals " << std::endl;
  #endif
  auto& tDof = com_mod.tDof;
  auto& dFlag = com_mod.dFlag;
  auto& nFacesLS = com_mod.nFacesLS;
  auto& recLn = com_mod.recLn;

  tDof = 0;
  dFlag = false;

  // Set faces for linear solver
  //
  #ifdef debug_initialize
  dmsg << "Set faces for linear solver " << std::endl;
  #endif
  nFacesLS = 0;
  for (auto& eq : com_mod.eq) {
    nFacesLS += eq.nBc;
  }

  // Remove LS pointer for faces with weakly applied Dir. BC
  #ifdef debug_initialize
  dmsg << "Remove LS pointer for faces" << std::endl;
  #endif
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    for (int i = 0; i < eq.nBc; i++) {
      if (eq.bc[i].weakDir) {
        nFacesLS = nFacesLS - 1;
      }
    }
  }

  // For FSI simulations, LS pointer to structure nodes
  if (com_mod.mvMsh) {
    nFacesLS = nFacesLS + 1;
  }

  // For initializing CMM, LS pointer to fixed edge nodes
  if (com_mod.cmmInit) {
    nFacesLS = nFacesLS + 1;
  }

  for (auto& bc : com_mod.eq[0].bc) {
    if (bc.cplBCptr != -1) { 
      com_mod.cplBC.coupled = true;
      break; 
    }
  }

  #ifdef debug_initialize
  dmsg << "nFacesLS: " << nFacesLS;
  dmsg << "com_mod.cplBC.coupled: " << com_mod.cplBC.coupled;
  #endif
  bool flag = false;

  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];

    // Set eq dof, sym and possibly am data memebers.
    //
    // This replaces the Fortran 'SELECT CASE (eq(iEq)%phys)' code.
    //
    std::tie(eq.dof, eq.sym) = equation_dof_map.at(eq.phys);

    if (std::set<EquationType>{Equation_fluid, Equation_heatF, Equation_heatS, Equation_CEP, Equation_stokes}.count(eq.phys) == 0) {
      dFlag = true;
    }

    // For second order eqs. 
    if (std::set<EquationType>{Equation_lElas, Equation_struct, Equation_shell, Equation_mesh}.count(eq.phys) != 0) {
      eq.am = (2.0 - eq.roInf) / (1.0 + eq.roInf);

    // first order equations.
    } else {
      eq.am = 0.5 * (3.0 - eq.roInf) / (1.0 + eq.roInf);
    }

    if ((eq.phys == Equation_CMM) && com_mod.cmmInit) { 
      eq.dof = nsd;
    }

    // This code checks to see if the heatF equation is accompanied
    // by a fluid equation. If not, it adds the appropriate number of
    // degrees of freedom based on the precomputed state-variable (velocity)
    // data.
    if (eq.phys == Equation_heatF) {
        bool fflag = false;
        for (int jEq = 0; jEq < com_mod.nEq; jEq++) {
            if (std::set < EquationType >
                {Equation_fluid, Equation_FSI, Equation_CMM, Equation_stokes}.count(com_mod.eq[jEq].phys)) {
                fflag = true;
            }
        }
        if (com_mod.usePrecomp) {
            if (!fflag) {
                tDof = tDof + nsd;
            }
        } else {
            if (!fflag) {
                throw std::runtime_error(
                        "HeatF equation must be accompanied by a fluid equation or precomputed velocity data.");
            }
        }
    }
    eq.pNorm = std::numeric_limits<double>::max();
    eq.af = 1.0 / (1.0 + eq.roInf);
    eq.beta = 0.25 * pow((1.0 + eq.am - eq.af), 2.0);
    eq.gam = 0.5 + eq.am - eq.af;

    // These are indexes into arrays so need to be zero-based.

    eq.s = tDof;
    eq.e = tDof + eq.dof - 1;
    tDof = eq.e + 1;
    if (eq.useTLS) {
      flag = true;
    }

    #ifdef debug_initialize
    dmsg << "dFlag: " << dFlag;
    dmsg << "eq.dof: " << eq.dof;
    dmsg << "eq.sym: " << eq.sym;
    dmsg << "eq.s: " << eq.s;
    dmsg << "eq.pNorm: " << eq.pNorm;
    dmsg << "eq.af: " << eq.af;
    dmsg << "eq.beta: " << eq.beta;
    dmsg << "eq.gam: " << eq.gam;
    dmsg << "tDof: " << tDof;
    #endif
  }

  int ierr = 0;
  if (dFlag) {
    ierr = 1;
  }

  int i = 0;
  if (com_mod.cplBC.coupled) {
    i = com_mod.cplBC.nX;
  }

  com_mod.stamp = {cm.np(), com_mod.nEq, com_mod.nMsh, com_mod.tnNo, i, tDof, ierr};

  // Calculating the record length
  //

  i = 2*tDof;
  if (dFlag) i = 3*tDof;
  if (com_mod.pstEq) i = i + com_mod.nsymd;
  if (com_mod.sstEq) i = i + nsd;
  if (cep_mod.cepEq) {
    i = i + cep_mod.nXion;
    if (cep_mod.cem.cpld) i = i + 1;
  }

  i = sizeof(int)*(1+com_mod.stamp.size()) + sizeof(double)*(2 + com_mod.nEq + com_mod.cplBC.nX + i*com_mod.tnNo);

  if (com_mod.ibFlag) {
    i = i + sizeof(double)*(3*nsd + 1) * com_mod.ib.tnNo;
  }

  if (cm.seq()) {
    recLn = i;
  } else { 
    MPI_Allreduce(&i, &recLn, 1, cm_mod::mpint, MPI_MAX, cm.com());
  }

  // Initialize shell eIEN data structure. Used later in LHSA.
  //

  // [TODO:DaveP] this is not implemented.
  //
  if (com_mod.shlEq) {
    for (int iM = 0; iM < com_mod.nMsh; iM++) { 
      auto& msh = com_mod.msh[iM];
      if (msh.lShl) {
        if (msh.eType == ElementType::NRB) {
          msh.sbc.resize(msh.eNoN, msh.nEl);
          //ALLOCATE(msh(iM)%eIEN(0,0))
          //ALLOCATE(msh(iM)%sbc(msh(iM)%eNoN,msh(iM)%nEl))
          msh.sbc = 0;
        } else if (msh.eType == ElementType::TRI3) {
          baf_ini_ns::set_shl_xien(simulation, msh);
          //CALL SETSHLXIEN(msh(iM))
        }
      }
    }
  }

  // Initialize tensor operations
  //
  mat_fun::ten_init(nsd);

  // Constructing stiffness matrix.
  //
  int nnz = 0;
  lhsa_ns::lhsa(simulation, nnz);

  int gnnz = nnz;
  MPI_Allreduce(&nnz, &gnnz, 1, cm_mod::mpint, MPI_SUM, cm.com());

  // Initialize FSILS structures
  //
  fsi_linear_solver::FSILS_commuType communicator;

  if (com_mod.resetSim) {
    if (communicator.foC) {
      fsils_commu_free(communicator);
    }

    if (com_mod.lhs.foC) {
      fsils_lhs_free(com_mod.lhs);
    }
  }

  fsi_linear_solver::fsils_commu_create(communicator, cm.com());

  fsi_linear_solver::fsils_lhs_create(com_mod.lhs, communicator, com_mod.gtnNo, com_mod.tnNo, nnz, 
      com_mod.ltg, com_mod.rowPtr, com_mod.colPtr, nFacesLS);

  // Variable allocation and initialization
  int tnNo = com_mod.tnNo; 
  com_mod.Ao.resize(tDof,tnNo); 
  com_mod.An.resize(tDof,tnNo); 
  com_mod.Yo.resize(tDof,tnNo); 
  com_mod.Yn.resize(tDof,tnNo); 
  com_mod.Do.resize(tDof,tnNo); 
  com_mod.Dn.resize(tDof,tnNo); 
  com_mod.Bf.resize(nsd,tnNo);

  // [TODO] DaveP not implemented?
  // IF (ibFlag) CALL IB_MEMALLOC() 

  // Additional physics dependent variables
  // USTRUCT phys
  //
  if (com_mod.sstEq) { 
    com_mod.Ad.resize(nsd,tnNo); 
    com_mod.Rd.resize(nsd,tnNo); 
    com_mod.Kd.resize((nsd+1)*nsd, nnz); 
  }

  // PRESTRESS
  if (com_mod.pstEq) {
    int nsymd = com_mod.nsymd; 
    com_mod.pS0.resize(nsymd,tnNo); 
    com_mod.pSn.resize(nsymd,tnNo); 
    com_mod.pSa.resize(tnNo); 
  } 

  // Electrophysiology
  //
  if (cep_mod.cepEq) {
    cep_mod.Xion.resize(cep_mod.nXion,tnNo);
    cep_ion::cep_init(simulation);

    // Electro-Mechanics
    if (cep_mod.cem.cpld) {
      cep_mod.cem.Ya.resize(tnNo);
    }
  }

  // Setup data for remeshing.
  //
  auto& rmsh = com_mod.rmsh;
  int nMsh = com_mod.nMsh;

  if (!com_mod.resetSim) {
    if (rmsh.flag.size() == 0) {
      rmsh.flag.resize(nMsh);
    }
    std::fill(rmsh.flag.begin(), rmsh.flag.end(), false);
    rmsh.fTS = rmsh.freq;

    if (rmsh.isReqd) {
      rmsh.A0.resize(com_mod.tDof,com_mod.tnNo);
      rmsh.Y0.resize(com_mod.tDof,com_mod.tnNo);
      rmsh.D0.resize(com_mod.tDof,com_mod.tnNo);
      rmsh.iNorm.resize(com_mod.nEq);
    }

    if (FILE *file = fopen(com_mod.iniFilePath.c_str(), "r")) {
      fclose(file);
      flag = true;
    } else {
      flag = false;
    }

    if (flag) { 
      auto& iniFilePath = com_mod.iniFilePath;
      auto& timeP = com_mod.timeP;

      if (iniFilePath.find(".bin") != std::string::npos) { 
        init_from_bin(simulation, iniFilePath, timeP);
      } else {
        init_from_vtu(simulation, iniFilePath, timeP);
     }

    } else {
      if (com_mod.stFileFlag) {
        std::string fTmp = com_mod.stFileName + "_last.bin";

        if (FILE *file = fopen(fTmp.c_str(), "r")) {
          fclose(file);
          auto& timeP = com_mod.timeP;
          init_from_bin(simulation, fTmp, timeP);
        } else {
          if (cm.mas(cm_mod)) {
            std::cout << "WARNING: No '" + fTmp + "' file to restart simulation from; Resetting restart flag to false";
          }
          com_mod.stFileFlag = false;
          zero_init(simulation);
        }

        if (rmsh.isReqd) {
          auto& cTS = com_mod.cTS;
          rmsh.fTS = (cTS/rmsh.fTS + 1)*rmsh.freq;
          rmsh.rTS = cTS;
          rmsh.time = com_mod.time;
          for (int i = 0; i < rmsh.iNorm.size(); i++) {
            rmsh.iNorm(i) = com_mod.eq[i].iNorm;
          }
          rmsh.A0 = com_mod.Ao;
          rmsh.Y0 = com_mod.Yo;
          rmsh.D0 = com_mod.Do;
        }

      } else {
        zero_init(simulation);
      } 
    }

    com_mod.rsTS = com_mod.cTS;

  } else {
    com_mod.cTS  = rmsh.rTS;
    com_mod.time = rmsh.time;
    for (int i = 0; i < com_mod.eq.size(); i++) {
      com_mod.eq[i].iNorm = rmsh.iNorm[i];
    }

    com_mod.Ao = all_fun::local(com_mod, cm_mod, cm, rmsh.A0);
    com_mod.Yo = all_fun::local(com_mod, cm_mod, cm, rmsh.Y0);
    com_mod.Do = all_fun::local(com_mod, cm_mod, cm, rmsh.D0);

    rmsh.A0.resize(tDof,tnNo); 
    rmsh.A0 = com_mod.Ao;

    rmsh.Y0.resize(tDof,tnNo); 
    rmsh.Y0 = com_mod.Yo;

    rmsh.D0.resize(tDof,tnNo); 
    rmsh.D0 = com_mod.Do;
  } // resetSim

  // Initialize new variables
  //
  com_mod.An = com_mod.Ao;
  com_mod.Yn = com_mod.Yo;
  com_mod.Dn = com_mod.Do;

  for (int iM = 0; iM < nMsh; iM++) { 
    if (cm.mas(cm_mod)) {
      std::string fTmp = chnl_mod.appPath + "partitioning_" + com_mod.msh[iM].name + ".bin";
      std::string sTmp = chnl_mod.appPath + "partitioning_" + com_mod.msh[iM].name + "_" + std::to_string(com_mod.cTS) + ".bin";

      if (FILE *file = fopen(fTmp.c_str(), "r")) {
        fclose(file);
        std::ifstream src(fTmp, std::ios::binary);
        std::ofstream dest(sTmp, std::ios::binary);
        dest << src.rdbuf();
      }
    }
  }

  // Initialize function spaces 
  //
  for (int iM = 0; iM < nMsh; iM++) { 
    auto& mesh = com_mod.msh[iM];
    fs::init_fs_msh(com_mod, mesh);
    for (int iFa = 0; iFa < com_mod.msh[iM].nFa; iFa++) { 
      fs::init_fs_face(com_mod, mesh, mesh.fa[iFa]);
    }
  }

  // Initialize Immersed Boundary data structures
  // [TODO:DaveP] not implemented but still need to allocate iblank.
  //
  com_mod.iblank.resize(tnNo);
  /* 
  ALLOCATE(iblank(tnNo))
  iblank = 0
  IF (ibFlag) CALL IB_INIT(Do)
  */

  // Calculating the volume of each domain
  //
  Array<double> s(1,tnNo);
  s = 1.0;
  
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    if (!com_mod.shlEq && !com_mod.cmmInit) {
      //std = " Eq. <"// CLR(eq(iEq)%sym, iEq)//">"
    }

    for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
      int i = eq.dmn[iDmn].Id;
      eq.dmn[iDmn].v = all_fun::integ(com_mod, cm_mod, i, s, 0, 0);
      if (!com_mod.shlEq && !com_mod.cmmInit) {
        //std = "    Volume of domain <"//STR(i)//"> is "// 2            STR(eq(iEq)%dmn(iDmn)%v)
        //IF (ISZERO(eq(iEq)%dmn(iDmn)%v)) wrn = "<< Volume of "// "domain "//iDmn//" of equation "//iEq//" is zero >>"
      }     
    } 
  }

  // Preparing faces and BCs
  //
  baf_ini_ns::baf_ini(simulation);

  // As all the arrays are allocated, call BIN to VTK for conversion
  //
  // This is set using the 'Convert_BIN_to_VTK_format' option.
  //
  if (com_mod.bin2VTK) {
    post::ppbin2vtk(simulation);
  }

  // Making sure the old solution satisfies BCs
  //
  // Modifes
  //  com_mod.Ao - Old time derivative of variables (acceleration)
  //  com_mod.Yo - Old variables (velocity)
  //  com_mod.Do - Old integrated variables (dissplacement)
  //
  set_bc::set_bc_dir(com_mod, com_mod.Ao, com_mod.Yo, com_mod.Do);

  // Preparing TXT files
  txt_ns::txt(simulation, true);

  // Printing the first line and initializing timeP
  int co = 1;
  int iEq = 0;
  output::output_result(simulation, com_mod.timeP, co, iEq);

  std::fill(com_mod.rmsh.flag.begin(), com_mod.rmsh.flag.end(), false);
  com_mod.resetSim = false;
}

//-----------
// zero_init
//-----------
// Initialize state variables Yo, Ao and Do.
//
void zero_init(Simulation* simulation)
{
  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  const int nsd = com_mod.nsd;

  #define n_debug_zero_init
  #ifdef debug_zero_init
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Initialize precomputed state variables
  //

  if (com_mod.usePrecomp) {
    for (int l = 0; l < com_mod.nMsh; l++) {
      auto& msh = com_mod.msh[l];
      for (int a = 0; a < com_mod.tnNo; a++) {
        // In the future this should depend on the equation type.
        for (int i = 0; i < nsd; i++) {
          com_mod.Yo(i,a) = msh.Ys(i,a,0);
        }
      }
    }
  }

  // Load any explicitly provided solution variables
  //
  if (com_mod.Vinit.size() != 0) {
     #ifdef debug_zero_init
     dmsg << "Initialize Yo to provided Vel solution";
     #endif
     for (int a = 0; a < com_mod.tnNo; a++) {
       for (int i = 0; i < nsd; i++) {
         com_mod.Yo(i,a) = com_mod.Vinit(i,a);
       }
     }
  }

  if (com_mod.Pinit.size() != 0) {
     #ifdef debug_zero_init
     dmsg << "Initialize Yo to provided P solution";
     #endif
     for (int a = 0; a < com_mod.tnNo; a++) {
       for (int i = 0; i < nsd; i++) {
         com_mod.Yo(nsd,a) = com_mod.Pinit(a);
       }
     }
  }

  if (com_mod.Dinit.size() != 0) {
     #ifdef debug_zero_init
     dmsg << "Initialize Do to provided D solution";
     #endif
     for (int a = 0; a < com_mod.tnNo; a++) {
       for (int i = 0; i < nsd; i++) {
         com_mod.Do(i,a) = com_mod.Dinit(i,a);
       }
     }
  }
}

