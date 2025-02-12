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

#include "set_bc.h"

#include "all_fun.h"
#include "cmm.h"
#include "consts.h"
#include "eq_assem.h"
#include "fft.h"
#include "fluid.h"
#include "fs.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "ustruct.h"
#include "utils.h"
#include <math.h>
#include "svZeroD_subroutines.h"

namespace set_bc {

/// @brief This function calculates updated cplBC pressures or flowrates from 0D,
/// as well as the resistance matrix M ~ dP/dQ from 0D using finite difference.
/// Updates the pressure or flowrates stored in cplBC.fa[i].y and the resistance
/// matrix M ~ dP/dQ stored in eq.bc[iBc].r.
/// @param com_mod 
/// @param cm_mod 
void calc_der_cpl_bc(ComMod& com_mod, const CmMod& cm_mod)
{
  using namespace consts;

  #define n_debug_calc_der_cpl_bc 
  #ifdef debug_calc_der_cpl_bc 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int iEq = 0;
  // NOTE: For coupling with svZeroDPlus, absTol needs to be > 1e-8 to be compatible with the default convergence tolerance of svZeroDPlus (1e-8)
  // If this is not true, the finite difference computation of bc.r below results in zero because the perturbation is below the svZeroDPlus tolerance 
  const double absTol = 1.0e-7;
  const double relTol = 1.0e-5;

  int nsd = com_mod.nsd;
  auto& eq = com_mod.eq[iEq];
  auto& cplBC = com_mod.cplBC;

  // If coupling is all with Dirichlet faces, no derivative calculation is needed
  // (see Moghadam et al. 2013 Section 2.2.2)
  if (std::count_if(cplBC.fa.begin(),cplBC.fa.end(),[](cplFaceType& fa){return fa.bGrp == CplBCType::cplBC_Dir;}) == cplBC.fa.size()) { 
    #ifdef debug_calc_der_cpl_bc 
    dmsg << "all cplBC_Dir " << std::endl;
    #endif
    return;
  }
  // if (ALL(cplBC.fa.bGrp .EQ. cplBC_Dir)) RETURN

  // Determine current physics
  auto cPhys = eq.phys;

  // Mechanical configuration in which to compute flowrate
  auto cfg_o = MechanicalConfigurationType::reference;
  auto cfg_n = MechanicalConfigurationType::reference;

  bool RCRflag = false;

  // Loop over BCs
  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    #ifdef debug_calc_der_cpl_bc 
    dmsg << "----- iBc " << iBc+1 << " -----" << std::endl;
    #endif
    auto& bc = eq.bc[iBc];
    int iFa = bc.iFa;
    int iM = bc.iM;
    int ptr = bc.cplBCptr;
    #ifdef debug_calc_der_cpl_bc 
    dmsg << "iFa: " << iFa;
    dmsg << "iM: " << iM;
    dmsg << "ptr: " << ptr;
    #endif

    if (utils::btest(bc.bType, iBC_RCR)) {
      if (!RCRflag) {
        RCRflag = true;
      }
    }

    if (ptr != -1) {
      auto& fa = com_mod.msh[iM].fa[iFa];
      // Compute flowrates at 3D Neumann boundaries at timesteps n and n+1
      if (utils::btest(bc.bType, iBC_Neu)) {
        // If struct or ustruct, use old and new configurations to compute flowrate integral
        if ((cPhys == EquationType::phys_struct) || (cPhys == EquationType::phys_ustruct)) {

          // Must use follower pressure load for 0D coupling with struct/ustruct
          if (!bc.flwP) { 
            throw std::runtime_error("[calc_der_cpl_bc]  Follower pressure load must be used for 0D coupling with struct/ustruct");
          }
          cfg_o = MechanicalConfigurationType::old_timestep;
          cfg_n = MechanicalConfigurationType::new_timestep;
        }
        // If fluid, FSI, or CMM, use reference configuration to compute flowrate integral
        // Note that for FSI, mvMsh will modify geometry in gnnb()
        else if ((cPhys == EquationType::phys_fluid) || (cPhys == EquationType::phys_FSI) || (cPhys == EquationType::phys_CMM)) {
          cfg_o = MechanicalConfigurationType::reference;
          cfg_n = MechanicalConfigurationType::reference;
        }
        else {
          throw std::runtime_error("[calc_der_cpl_bc]  Invalid physics type for 0D coupling");
        }
        cplBC.fa[ptr].Qo = all_fun::integ(com_mod, cm_mod, fa, com_mod.Yo, 0, nsd-1, false, cfg_o);
        cplBC.fa[ptr].Qn = all_fun::integ(com_mod, cm_mod, fa, com_mod.Yn, 0, nsd-1, false, cfg_n);
        cplBC.fa[ptr].Po = 0.0;
        cplBC.fa[ptr].Pn = 0.0;
        #ifdef debug_calc_der_cpl_bc 
        dmsg << "iBC_Neu ";
        dmsg << "cplBC.fa[ptr].Qo: " << cplBC.fa[ptr].Qo;
        dmsg << "cplBC.fa[ptr].Qn: " << cplBC.fa[ptr].Qn;
        #endif

      }
      // Compute avg pressures at 3D Dirichlet boundaries at timesteps n and n+1 
      else if (utils::btest(bc.bType, iBC_Dir)) {
        double area = fa.area;
        cplBC.fa[ptr].Po = all_fun::integ(com_mod, cm_mod, fa, com_mod.Yo, nsd) / area;
        cplBC.fa[ptr].Pn = all_fun::integ(com_mod, cm_mod, fa, com_mod.Yn, nsd) / area;
        cplBC.fa[ptr].Qo = 0.0;
        cplBC.fa[ptr].Qn = 0.0;
        #ifdef debug_calc_der_cpl_bc 
        dmsg << "iBC_Dir" << std::endl;
        dmsg << "cplBC.fa[ptr].Po: " << cplBC.fa[ptr].Po;
        dmsg << "cplBC.fa[ptr].Pn: " << cplBC.fa[ptr].Pn;
        #endif
      }
    }
  }

  #ifdef debug_calc_der_cpl_bc 
  dmsg << "RCRflag: " << RCRflag;
  #endif

  // Call genBC or cplBC to get updated pressures or flowrates.
  if (cplBC.useGenBC) {
     set_bc::genBC_Integ_X(com_mod, cm_mod, "D");
   } else if (cplBC.useSvZeroD) {
     svZeroD::calc_svZeroD(com_mod, cm_mod, 'D');
   } else {
     set_bc::cplBC_Integ_X(com_mod, cm_mod, RCRflag);
  }

  // Compute the epsilon parameter (diff) for the finite difference calculation
  // of the resistance matrix M ~ dP/dQ. Slightly different from Eq. 30 in Moghadam et al. 2013
  int j = 0;
  double diff = 0.0;

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int i = bc.cplBCptr;
    if (i != -1 && utils::btest(bc.bType, iBC_Neu)) {
      diff = diff + (cplBC.fa[i].Qn * cplBC.fa[i].Qn);
      j = j + 1;
    }
  }

  diff = sqrt(diff / static_cast<double>(j));
  if (diff*relTol < absTol) {
     diff = absTol;
   } else {
     diff = diff*relTol;
  }

  // Store the original pressures and flowrates
  std::vector<double> orgY(cplBC.fa.size());
  std::vector<double> orgQ(cplBC.fa.size());

  for (size_t i = 0; i < cplBC.fa.size(); i++) {
    orgY[i] = cplBC.fa[i].y;
    orgQ[i] = cplBC.fa[i].Qn;
  }

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int i = bc.cplBCptr;

    if (i != -1 && utils::btest(bc.bType, iBC_Neu)) {

        // Finite difference perturbation in flowrate
        cplBC.fa[i].Qn = orgQ[i] + diff;

        // Call genBC or cplBC again with perturbed flowrate
        if (cplBC.useGenBC) {
          set_bc::genBC_Integ_X(com_mod, cm_mod, "D");
        } else if (cplBC.useSvZeroD) {
          svZeroD::calc_svZeroD(com_mod, cm_mod, 'D');
        } else {
          set_bc::cplBC_Integ_X(com_mod, cm_mod, RCRflag);
        }

        // Finite difference calculation of the resistance dP/dQ
        bc.r = (cplBC.fa[i].y - orgY[i]) / diff;

        // Restore the original pressures and flowrates
        for (size_t j = 0; j < cplBC.fa.size(); j++) {
          cplBC.fa[j].y = orgY[j];
          cplBC.fa[j].Qn = orgQ[j];
}
     }
  }
}

/// @brief Interface to call 0D code (cplBC)
// 
void cplBC_Integ_X(ComMod& com_mod, const CmMod& cm_mod, const bool RCRflag)
{
  using namespace consts;

  int nsd = com_mod.nsd;
  auto& cplBC = com_mod.cplBC;
  auto& cm = com_mod.cm;
  int istat = 0;

  if (cm.mas(cm_mod)) {
    istat = 0;

    if (RCRflag) {
      RCR_Integ_X(com_mod, cm_mod, istat);
    } else {
      throw std::runtime_error("Interface to 0D code is not implemented.");
      //int fid = 1;
      //OPEN(fid, FILE=cplBC.commuName, FORM='UNFORMATTED')
      //WRITE(fid) cplBC.nFa
      //WRITE(fid) cplBC.nX
      //WRITE(fid) cplBC.nXp
      //WRITE(fid) dt
      //WRITE(fid) MAX(time-dt, 0._RKIND)
      //WRITE(fid) cplBC.xo

      for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
        //WRITE(fid) cplBC.fa(iFa).bGrp
        //WRITE(fid) cplBC.fa(iFa).Qo
        //WRITE(fid) cplBC.fa(iFa).Qn
        //WRITE(fid) cplBC.fa(iFa).Po
        //WRITE(fid) cplBC.fa(iFa).Pn
        //WRITE(fid) cplBC.fa(iFa).name
      }
      //CLOSE(fid)

      //CALL SYSTEM(TRIM(cplBC.binPath)//" "//TRIM(cplBC.commuName))

      //OPEN(fid,FILE=TRIM(cplBC.commuName),STATUS='OLD', FORM='UNFORMATTED')
      //READ(fid) istat
      //READ(fid) cplBC.xn
      //READ(fid) cplBC.xp

      for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
        //READ(fid) cplBC.fa(iFa).y
      }
      //CLOSE(fid)
    }
  }

  cm.bcast(cm_mod, &istat);

  if (istat != 0) {
    if (RCRflag) {
      throw std::runtime_error("RCR integration error detected, Aborting!");
    } else {
      throw std::runtime_error("CPLBC Error detected, Aborting!");
    }
  }

  if (!cm.seq()) {
    Vector<double> y(cplBC.nFa);

    if (cm.mas(cm_mod)) {
      for (int i = 0; i < cplBC.nFa; i++) {
        y(i) = cplBC.fa[i].y;
      }
    }

    cm.bcast(cm_mod, cplBC.xn);
    cm.bcast(cm_mod, y);

    if (cm.slv(cm_mod)) { 
      for (int i = 0; i < cplBC.nFa; i++) {
        cplBC.fa[i].y = y(i);
      }
    }
  }
}


//---------------
// genBC_Integ_X
//---------------
// Interface to call 0D code (genBC/gcode)
//
//
void genBC_Integ_X(ComMod& com_mod, const CmMod& cm_mod, const std::string& genFlag)
{
  using namespace consts;

  int nDir = 0;
  int nNeu = 0;
  double dt = com_mod.dt;
  auto& cplBC = com_mod.cplBC;
  auto& cm = com_mod.cm;

  int len = cplBC.binPath.size() + cplBC.commuName.size() + 1;
  char command[len];
  strcpy(command, cplBC.binPath.c_str());
  strcat(command, " ");
  strcat(command, cplBC.commuName.c_str());

  // If this process is the master process on the communicator
  if (cm.mas(cm_mod)) {
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      auto& fa = cplBC.fa[iFa];

      if (fa.bGrp == CplBCType::cplBC_Dir) {
        nDir = nDir + 1;
      } else if (fa.bGrp == CplBCType::cplBC_Neu) {
        nNeu = nNeu + 1;
      }
    }

    // Write coupling info (number of Dirichlet and Neumann surfaces, pressure, 
    // flow rate) from 3D to cplBC communication file (for GenBC, usually 
    // called GenBC.int)
    int int_size = sizeof(int);
    int double_size = sizeof(double);
    int flag_size = genFlag.length();

    std::ofstream genBC_writer;
    genBC_writer.open(cplBC.commuName, std::ios::out|std::ios::binary);
    if (!genBC_writer.is_open()) {
      throw std::runtime_error("Failed to open the genBC initialization file '" + cplBC.commuName + "' to write.");
    }
    // Flag for how genBC behaves (I: Initializing, T: Iteration loop, L: Last iteration, D: Derivative)
    genBC_writer.write( (char*)&flag_size, int_size);
    genBC_writer.write(genFlag.c_str(), flag_size);
    genBC_writer.write( (char*)&flag_size, int_size);
    
    genBC_writer.write( (char*)&double_size, int_size);
    genBC_writer.write( (char*)&dt, double_size);
    genBC_writer.write( (char*)&double_size, int_size);
    
    genBC_writer.write( (char*)&int_size, int_size);
    genBC_writer.write( (char*)&nDir, int_size);
    genBC_writer.write( (char*)&int_size, int_size);
    
    genBC_writer.write( (char*)&int_size, int_size);
    genBC_writer.write( (char*)&nNeu, int_size);
    genBC_writer.write( (char*)&int_size, int_size);

    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      if (cplBC.fa[iFa].bGrp == CplBCType::cplBC_Dir) {
        genBC_writer.write( (char*)&double_size, int_size);
        genBC_writer.write( (char*)&cplBC.fa[iFa].Po, double_size);
        genBC_writer.write( (char*)&double_size, int_size);
        genBC_writer.write( (char*)&double_size, int_size);
        genBC_writer.write( (char*)&cplBC.fa[iFa].Pn, double_size);
        genBC_writer.write( (char*)&double_size, int_size);
      }
    }
    
    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      if (cplBC.fa[iFa].bGrp == CplBCType::cplBC_Neu) {
        genBC_writer.write( (char*)&double_size, int_size);
        genBC_writer.write( (char*)&cplBC.fa[iFa].Qo, double_size);
        genBC_writer.write( (char*)&double_size, int_size);
        genBC_writer.write( (char*)&double_size, int_size);
        genBC_writer.write( (char*)&cplBC.fa[iFa].Qn, double_size);
        genBC_writer.write( (char*)&double_size, int_size);
      }
    }
    genBC_writer.close();
    
    // Call genBC executable which reads the communication file GenBC.int
    system(command);

    // Read outputs from genBC, which are in the same GenBC.int
    std::ifstream genBC_reader(cplBC.commuName, std::ios::out | std::ios::binary);
    if (!genBC_reader.is_open()) {
      throw std::runtime_error("Failed to open the genBC interface file '" + cplBC.commuName + "' to read.");
    }

    int size_buffer;

    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      if (cplBC.fa[iFa].bGrp == CplBCType::cplBC_Dir) {
        genBC_reader.read( (char*)&size_buffer, int_size );
        genBC_reader.read( (char*)&cplBC.fa[iFa].y, size_buffer );
        genBC_reader.read( (char*)&size_buffer, int_size );
      }
    }

    for (int iFa = 0; iFa < cplBC.nFa; iFa++) {
      if (cplBC.fa[iFa].bGrp == CplBCType::cplBC_Neu) {
        genBC_reader.read( (char*)&size_buffer, int_size );
        genBC_reader.read( (char*)&cplBC.fa[iFa].y, size_buffer );
        genBC_reader.read( (char*)&size_buffer, int_size );
      }
    }

    genBC_reader.close();
  }

  // If there are multiple procs (not sequential), broadcast genBC outputs to
  // follower procs
  if (!cm.seq()) {
    Vector<double> y(cplBC.nFa);

    if (cm.mas(cm_mod)) {
      for (int i = 0; i < cplBC.nFa; i++) {
        y(i) = cplBC.fa[i].y;
      }
    }

    cm.bcast(cm_mod, y);

    if (cm.slv(cm_mod)) {
      for (int i = 0; i < cplBC.nFa; i++) {
        cplBC.fa[i].y = y(i);
      }
    }
  }

}

void RCR_Integ_X(ComMod& com_mod, const CmMod& cm_mod, int istat)
{
  using namespace consts;
  const int nTS = 100;

  int nsd = com_mod.nsd;
  auto& cplBC = com_mod.cplBC;
  double time = com_mod.time;
  double dt = com_mod.dt;

  double tt = fmax(time - dt, 0.0);
  double dtt = dt / static_cast<double>(nTS);
  int nX = cplBC.nFa;

  Vector<double> Rp(nX), C(nX), Rd(nX), Pd(nX); 
  Vector<double> X(nX), Xrk(nX); 
  Array<double> frk(nX,4), Qrk(nX,4);

  for (int i = 0; i < nX; i++) {
    Rp(i) = cplBC.fa[i].RCR.Rp;
    C(i) = cplBC.fa[i].RCR.C;
    Rd(i) = cplBC.fa[i].RCR.Rd;
    Pd(i) = cplBC.fa[i].RCR.Pd;
  }
  X = cplBC.xo;

  for (int n = 0; n < nTS; n++) {
    for (int i = 0; i < 4; i++) {
      double r = static_cast<double>(i) / 3.0;
      r = (static_cast<double>(n) + r) / static_cast<double>(nTS);
      for (int j = 0; j < Qrk.nrows(); j++) {
        Qrk(j,i) = cplBC.fa[j].Qo + (cplBC.fa[j].Qn - cplBC.fa[j].Qo) * r;
      }   
    }

    // RK-4 1st pass
    //
    double trk = tt;
    Xrk = X;

    for (int j = 0; j < nX; j++) {
      frk(j,0) = (Qrk(j,0) - (Xrk(j) - Pd(j)) / Rd(j)) / C(j);
    }

    // RK-4 2nd pass
    trk = tt + dtt / 3.0;
    Xrk = X  + dtt * frk.col(0) / 3.0;

    for (int j = 0; j < Qrk.nrows(); j++) {
      frk(j,1) = (Qrk(j,1) - (Xrk(j)-Pd(j)) / Rd(j)) / C(j);
    }

    // RK-4 3rd pass
    trk = tt + 2.0 * dtt / 3.0;
    Xrk = X - dtt * frk.col(0) / 3.0  +  dtt * frk.col(1);
    for (int j = 0; j < Qrk.nrows(); j++) {
      frk(j,2) = (Qrk(j,2) - (Xrk(j) - Pd(j)) / Rd(j)) / C(j);
    }

    // RK-4 4th pass
    trk = tt + dtt;
    Xrk = X  + dtt * frk.col(0)  -  dtt * frk.col(1)  +  dtt * frk.col(2);
    for (int j = 0; j < Qrk.nrows(); j++) {
      frk(j,3) = (Qrk(j,3) - (Xrk(j) - Pd(j)) / Rd(j)) / C(j);
    }

    double r = dtt / 8.0;
    X  = X + r*(frk.col(0) + 3.0*(frk.col(1) + frk.col(2)) + frk.col(3));
    tt = tt + dtt;

    for (int i = 0; i < nX; i++) {
      if (isnan(X(i))) {
        throw std::runtime_error("ERROR: NaN detected in RCR integration");
        istat = -1;
        return;
      }
    }
  }

  cplBC.xn = X;
  cplBC.xp(0) = tt;

  for (int i = 0; i < nX; i++) {
    cplBC.xp(i+1) = Qrk(i,3); //cplBC.fa(i).Qn
    cplBC.fa[i].y = X(i) + (cplBC.fa[i].Qn * Rp(i));
  }

}

/// @brief Initialize RCR variables (Xo) from flow field or using user-
/// provided input. This subroutine is called only when the simulation
/// is not restarted.
///
/// Replaces 'SUBROUTINE RCRINIT()'
//
void rcr_init(ComMod& com_mod, const CmMod& cm_mod)
{
  using namespace consts;

  const int iEq = 0;
  int nsd = com_mod.nsd;
  auto& eq = com_mod.eq[iEq];
  auto& cplBC = com_mod.cplBC;

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int iFa = bc.iFa;
    int iM  = bc.iM;
    int ptr = bc.cplBCptr;

    if (!utils::btest(bc.bType, iBC_RCR)) {
      continue; 
    }

    if (ptr != -1) {
      if (cplBC.initRCR) {
        auto& fa = com_mod.msh[iM].fa[iFa];
        double area = fa.area;
        double Qo = all_fun::integ(com_mod, cm_mod, fa, com_mod.Yo, 0, nsd-1);
        double Po = all_fun::integ(com_mod, cm_mod, fa, com_mod.Yo, nsd)  / area;
        cplBC.xo[ptr] = Po - (Qo * cplBC.fa[ptr].RCR.Rp);
      } else { 
        cplBC.xo[ptr] = cplBC.fa[ptr].RCR.Xo;
      }
    }
  }
}

/// @brief Below defines the SET_BC methods for the Coupled Momentum Method (CMM)
//
void set_bc_cmm(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Ag, const Array<double>& Dg ) 
{
  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];

    if (!utils::btest(bc.bType,iBC_CMM)) {
      continue; 
    }

    int iFa = bc.iFa;
    int iM = bc.iM;

    if (com_mod.msh[iM].eType != ElementType::TET4 && com_mod.msh[iM].fa[iFa].eType != ElementType::TRI3) {
      throw std::runtime_error("[set_bc_cmm] CMM equation is formulated for tetrahedral elements (volume) and triangular (surface) elements");
    }

    set_bc_cmm_l(com_mod, cm_mod, com_mod.msh[iM].fa[iFa], Ag, Dg);
  }
}

void set_bc_cmm_l(ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Array<double>& Ag, const Array<double>& Dg ) 
{
  using namespace consts;

  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  const auto& pS0 = com_mod.pS0;
  const auto& Bf = com_mod.Bf;

  int iM = lFa.iM;
  Array<double> al(tDof,3), dl(tDof,3), xl(3,3), bfl(3,3); 
  Vector<int> ptr(3);

  // Constructing the CMM contributions to the LHS/RHS and
  // assembling them
  //
  for (int e = 0; e < lFa.nEl; e++) {
    cDmn = all_fun::domain(com_mod, com_mod.msh[iM], cEq, lFa.gE(e));
    if (eq.dmn[cDmn].phys != EquationType::phys_CMM) {
      continue;
    } 

    Vector<double> pSl(6), vwp(2); 

    for (int a = 0; a < 3; a++) {
      int Ac = lFa.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) { 
        xl(i,a) = com_mod.x(i,Ac);
        al(i,a) = Ag(i,Ac);
        dl(i,a) = Dg(i,Ac);
        bfl(i,a) = Bf(i,Ac);
      }

      if (com_mod.pS0.size() != 0) {
        pSl = pSl + pS0.col(Ac);
      }

      if (com_mod.cmmVarWall) {
        vwp = vwp + com_mod.varWallProps.col(Ac);
      }
    }

    pSl = pSl / 3.0;
    vwp = vwp / 3.0;

    // Add CMM BCs contributions to the LHS/RHS
    cmm::cmm_b(com_mod, lFa, e, al, dl, xl, bfl, pSl, vwp, ptr);
  }

}

/// @brief Coupled BC quantities are computed here.
/// Reproduces the Fortran 'SETBCCPL()' subrotutine.
//
void set_bc_cpl(ComMod& com_mod, CmMod& cm_mod)
{
  static double absTol = 1.E-8, relTol = 1.E-5;

  using namespace consts;

  const int nsd = com_mod.nsd;
  auto& cplBC = com_mod.cplBC;
  auto& Yo = com_mod.Yo;
  auto& Yn = com_mod.Yn;
  const int iEq = 0;
  auto& eq = com_mod.eq[iEq];

  // Determine current physics
  auto cPhys = eq.phys;

  // Configuration in which to compute flowrate
  auto cfg_o = MechanicalConfigurationType::reference;
  auto cfg_n = MechanicalConfigurationType::reference;

  // If coupling scheme is implicit, calculate updated pressure and flowrate 
  // from 0D, as well as resistance from 0D using finite difference.
  if (cplBC.schm == CplBCType::cplBC_I) { 
    calc_der_cpl_bc(com_mod, cm_mod);

  // If coupling scheme is semi-implicit or explicit, only calculated updated
  // pressure and flowrate from 0D
  } else {
    bool RCRflag = false; 

    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      int iFa = bc.iFa;
      int iM  = bc.iM;
      int ptr = bc.cplBCptr;

      if (utils::btest(bc.bType,iBC_RCR)) {
        if (!RCRflag) {
          RCRflag = true;
        }
      }


      if (ptr != -1) {
        // Compute flowrates at 3D Neumann boundaries at timesteps n and n+1
        if (utils::btest(bc.bType,iBC_Neu)) {
          // If struct or ustruct, use old and new configurations to compute flowrate integral
          if ((cPhys == EquationType::phys_struct) || (cPhys == EquationType::phys_ustruct)) {

            // Must use follower pressure load for 0D coupling with struct/ustruct
            if (!bc.flwP) { 
              throw std::runtime_error("[set_bc_cpl]  Follower pressure load must be used for 0D coupling with struct/ustruct");
            }
            cfg_o = MechanicalConfigurationType::old_timestep;
            cfg_n = MechanicalConfigurationType::new_timestep;
          }
          // If fluid, FSI, or CMM, use reference configuration to compute flowrate integral
          // Note that for FSI, mvMsh will modify geometry in gnnb()
          else if ((cPhys == EquationType::phys_fluid) || (cPhys == EquationType::phys_FSI) || (cPhys == EquationType::phys_CMM)) {
            cfg_o = MechanicalConfigurationType::reference;
            cfg_n = MechanicalConfigurationType::reference;
          }
          else {
            throw std::runtime_error("[set_bc_cpl]  Invalid physics type for 0D coupling");
          }
        
          cplBC.fa[ptr].Qo = all_fun::integ(com_mod, cm_mod, com_mod.msh[iM].fa[iFa], Yo, 0, nsd-1, false, cfg_o);
          cplBC.fa[ptr].Qn = all_fun::integ(com_mod, cm_mod, com_mod.msh[iM].fa[iFa], Yn, 0, nsd-1, false, cfg_n);
          cplBC.fa[ptr].Po = 0.0;
          cplBC.fa[ptr].Pn = 0.0;
        } 
        // Compute avg pressures at 3D Dirichlet boundaries at timesteps n and n+1
        else if (utils::btest(bc.bType,iBC_Dir)) {
          double area = com_mod.msh[iM].fa[iFa].area;
          cplBC.fa[ptr].Po = all_fun::integ(com_mod, cm_mod, com_mod.msh[iM].fa[iFa], Yo, nsd) / area;
          cplBC.fa[ptr].Pn = all_fun::integ(com_mod, cm_mod, com_mod.msh[iM].fa[iFa], Yn, nsd) / area;
          cplBC.fa[ptr].Qo = 0.0;
          cplBC.fa[ptr].Qn = 0.0;
        }
      }
    }

    // Call genBC or cplBC to get updated pressures or flowrates.
    // Updates pressure or flowrates stored in cplBC.fa[i].y
    if (cplBC.useGenBC) {
       set_bc::genBC_Integ_X(com_mod, cm_mod, "D");
    } else if (cplBC.useSvZeroD){
      svZeroD::calc_svZeroD(com_mod, cm_mod, 'D');
    } else {
       set_bc::cplBC_Integ_X(com_mod, cm_mod, RCRflag);
    }
  }

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int iFa = bc.iFa;
    int ptr = bc.cplBCptr;
    if (ptr != -1) {
      bc.g = cplBC.fa[ptr].y;
    }
  }
}

/// @brief Apply Dirichlet BCs strongly.
///
/// Parameters
///   lA - New time derivative of variables (An)
///   lY - New variables (Yn)
///   lD - New integrated variables (Dn)
///
/// Modfies:
///   lA(tDof, tnNo)
///   lY(tDof, tnNo)
///   lD(tDof, tnNo)
///   com_mod.Ad - Time derivative of displacement
///
/// Reproduces 'SUBROUTINE SETBCDIR(lA, lY, lD)'
//
void set_bc_dir(ComMod& com_mod, Array<double>& lA, Array<double>& lY, Array<double>& lD)
{
  using namespace consts;

  #define n_set_bc_dir
  #ifdef set_bc_dir
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  int nEq = com_mod.nEq;
  std::vector<bool> eDir(maxNSD);

  for (int iEq = 0; iEq < nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    #ifdef set_bc_dir
    dmsg << ">>>> iEq: " << iEq;
    dmsg << "eq.nBc: " << eq.nBc;
    #endif

    for (int iBc = 0; iBc < eq.nBc; iBc++) {
      auto& bc = eq.bc[iBc];
      #ifdef set_bc_dir
      dmsg << ">> iBc " << iBc+1;
      #endif

      if (utils::btest(bc.bType, iBC_CMM)) {
        int s = eq.s;
        int e = eq.e;
        if (eq.dof == nsd+1) {
          e = e - 1;
        }
        int iFa = bc.iFa;
        int iM = bc.iM;

        for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
          if (utils::is_zero(bc.gx(a))) {
            int Ac = com_mod.msh[iM].fa[iFa].gN(a);
            for (int i = s; i <= e; i++) { 
              lA(i,Ac) = 0.0;
              lY(i,Ac) = 0.0;
            }
          }
        }
      } // END bType_CMM

      if (!utils::btest(bc.bType, iBC_Dir)) {
        continue;
      }

      if (bc.weakDir) {
        continue; 
      }
      int s = eq.s;
      int e = eq.e;
      if (eq.dof == nsd+1) {
        e = e - 1;
      }
      #ifdef set_bc_dir
      dmsg << ">> s: " << s;
      dmsg << ">> e: " << e;
      #endif
      std::fill(eDir.begin(), eDir.end(), false);
      int lDof = 0;

      for (int i = 0; i < nsd; i++) {
        if (bc.eDrn(i) != 0) {
          eDir[i] = true; 
          lDof = lDof + 1;
        }
      }

      if (lDof == 0) {
        lDof = e - s + 1;
      }
      int iFa = bc.iFa;
      int iM = bc.iM;
      int nNo = com_mod.msh[iM].fa[iFa].nNo;
      #ifdef set_bc_dir
      dmsg << ">> lDof: " << lDof;
      dmsg << ">> iM: " << iM;
      dmsg << ">> iFa: " << iFa;
      dmsg << ">> name: " << com_mod.msh[iM].fa[iFa].name ;
      dmsg << ">> nNo: " << nNo;
      #endif

      Array<double> tmpA(lDof,nNo); 
      Array<double> tmpY(lDof,nNo);

      // Modifies: tmpA, tmpY
      set_bc::set_bc_dir_l(com_mod, bc, com_mod.msh[iM].fa[iFa], tmpA, tmpY, lDof);

      if (std::find(eDir.begin(), eDir.end(), true) != eDir.end()) {
        if (utils::btest(bc.bType, enum_int(BoundaryConditionType::bType_impD))) {

          for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
            int Ac = com_mod.msh[iM].fa[iFa].gN(a);
            lDof = 0;

            for (int i = 0; i < nsd; i++) {
              if (eDir[i]) {
                lY(s+i,Ac) = tmpA(lDof,a);
                lD(s+i,Ac) = tmpY(lDof,a);
                lDof = lDof + 1;
              }
            }
          }

        } else {
          for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
            int Ac = com_mod.msh[iM].fa[iFa].gN(a);
            lDof = 0;
            for (int i = 0; i < nsd; i++) {
              if (eDir[i]) {
                lA(s+i,Ac) = tmpA(lDof,a);
                lY(s+i,Ac) = tmpY(lDof,a);
                lDof = lDof + 1;
              }
            }
          }
        }

      // No eDir[] is true. 
      //
      } else {
        if (utils::btest(bc.bType, enum_int(BoundaryConditionType::bType_impD))) {
          for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
            int Ac = com_mod.msh[iM].fa[iFa].gN(a);
            for (int i = 0; i < tmpA.nrows(); i++) {
              lY(i+s,Ac) = tmpA(i,a);
              lD(i+s,Ac) = tmpY(i,a);
            }
          }
        } else {
          for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
            int Ac = com_mod.msh[iM].fa[iFa].gN(a);
            for (int i = 0; i < lDof; i++) {
              lA(i+s,Ac) = tmpA(i,a);
              lY(i+s,Ac) = tmpY(i,a);
            }
          }
        }
      } // if (std::find(eDir.begin(), eDir.end(), true) != eDir.end())

      // if FSI and velocity-pressure based structural dynamics solver is used 
      // of nonlinear structure (v-p).
      //
      if ((eq.phys == EquationType::phys_FSI && com_mod.sstEq) || (eq.phys == EquationType::phys_ustruct)) {
        double c1  = eq.gam * com_mod.dt;
        double c1i = 1.0 / c1;
        double c2  = (eq.gam - 1.0)*com_mod.dt;

        if (std::find(eDir.begin(), eDir.end(), true) != eDir.end()) {
          if (utils::btest(bc.bType, enum_int(BoundaryConditionType::bType_impD))) {

            for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
              int Ac = com_mod.msh[iM].fa[iFa].gN(a);
              for (int i = 0; i < nsd; i++) {
                if (eDir[i]) {
                  int j = s + i;
                  lA(j,Ac) = c1i*(lY(j,Ac) - com_mod.Yo(j,Ac) + c2*com_mod.Ao(j,Ac));
                  com_mod.Ad(i,Ac) = c1i*(lD(j,Ac) - com_mod.Do(j,Ac) + c2*com_mod.Ad(i,Ac));
                }
              }
            }
          } else {
            for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
              int Ac = com_mod.msh[iM].fa[iFa].gN(a);
              for (int i = 0; i < nsd; i++) {
                if (eDir[i]) {
                  int j = s + i;
                  lD(j,Ac) = c1*lY(j,Ac) - c2*com_mod.Ad(i,Ac) + com_mod.Do(j,Ac);
                  com_mod.Ad(i,Ac) = lY(j,Ac);
                }
              }
            }
          }

        } else {
          if (utils::btest(bc.bType, enum_int(BoundaryConditionType::bType_impD))) {
            for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
              int Ac = com_mod.msh[iM].fa[iFa].gN(a);
              for (int i = 0; i < com_mod.Ad.nrows(); i++) {
                lA(i+s,Ac) = c1i*(lY(i+s,Ac) - com_mod.Yo(i+s,Ac) + c2*com_mod.Ao(i+s,Ac));
                com_mod.Ad(i,Ac) = c1i*(lD(i+s,Ac) - com_mod.Do(i+s,Ac) + c2*com_mod.Ad(i,Ac));
              }
            }

          } else {
            for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
              int Ac = com_mod.msh[iM].fa[iFa].gN(a);
              for (int i = 0; i < com_mod.Ad.nrows(); i++) {
                lD(i+s,Ac) = c1*lY(i+s,Ac) - c2*com_mod.Ad(i,Ac) + com_mod.Do(i+s,Ac);
                com_mod.Ad(i,Ac) = lY(i+s,Ac);
              }
            }
          }
        }
      }
    } // iBc
  } // iEq

}

/// Modifies:
///   lA(lDof,lFa.nNo)
///   lY(lDof,lFa.nNo)
///
/// Reproduces 'SUBROUTINE SETBCDIRL(lBc, lFa, lA, lY, lDof)'
//
void set_bc_dir_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa, Array<double>& lA, Array<double>& lY, int lDof)
{
  using namespace consts;

  #define n_debug_set_bc_dir_l
  #ifdef debug_set_bc_dir_l
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  int nEq = com_mod.nEq;

  double dirY, dirA;

  if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_gen))) {
    #ifdef debug_set_bc_dir_l
    dmsg << "bType_gen";
    #endif
    if (lDof != lBc.gm.dof) {
      throw std::runtime_error("[set_bc_dirl] Inconsistent DOF");
    }

    igbc(com_mod, lBc.gm, lY, lA);
    return;
  
  } else if (utils::btest(lBc.bType, enum_int(BoundaryConditionType::bType_ustd))) {
    #ifdef debug_set_bc_dir_l
    dmsg << "bType_ustd";
    #endif
    Vector<double> dirY_v(1), dirA_v(1);
    ifft(com_mod, lBc.gt, dirY_v, dirA_v);
    dirY = dirY_v(0);
    dirA = dirA_v(0);

  } else { 
    dirA = 0.0;
    dirY = lBc.g;
  }

  if (lDof == nsd) {
    for (int a = 0; a < lFa.nNo; a++) {

      for (int i = 0; i < lA.nrows(); i++) {
        double nV = lFa.nV(i,a);
        lA(i,a) = dirA * lBc.gx(a) * nV;
        lY(i,a) = dirY * lBc.gx(a) * nV;
      }
    }

  } else {
    for (int a = 0; a < lFa.nNo; a++) {
      for (int i = 0; i < lDof; i++) {
        lA(i,a) = dirA*lBc.gx(a);
        lY(i,a) = dirY*lBc.gx(a);
      }
    }
  }
}

/// @brief Weak treatment of Dirichlet boundary conditions
//
void set_bc_dir_w(ComMod& com_mod, const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int iM = bc.iM;
    int iFa = bc.iFa;
    if (!bc.weakDir) {
      continue;
    }
    set_bc_dir_wl(com_mod, bc, com_mod.msh[iM], com_mod.msh[iM].fa[iFa], Yg, Dg);
  }
}

/// @brief Reproduces Fortran 'SETBCDIRWL'.
//
void set_bc_dir_wl(ComMod& com_mod, const bcType& lBc, const mshType& lM, const faceType& lFa, const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  #define n_debug_set_bc_dir_wl
  #ifdef debug_set_bc_dir_wl
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;
  const int tnNo = com_mod.tnNo;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  const int nNo   = lFa.nNo;
  const int nEl   = lFa.nEl;
  const int eNoNb = lFa.eNoN;
  const int eNoN  = lM.eNoN;
  const auto& tauB  = lBc.tauB;

  bool flag = false;
  if (lM.nFs == 1) {
    flag = true; 
  } else {
    flag = false;
  }

  // Compute the Dirichlet value to be applied weakly on the face
  //
  int ss = eq.s;
  int ee = eq.e;

  if (eq.dof == nsd+1) {
    ee = ee - 1;
  }

  std::vector<bool> eDir(maxNSD);
  std::fill(eDir.begin(), eDir.end(), false);
  int lDof = 0;

  for (int i = 0; i < nsd; i++) {
    if (lBc.eDrn(i) != 0) {
      eDir[i] = true;
      lDof = lDof + 1;
     }
  }

  if (lDof == 0) {
    lDof = ee - ss + 1;
  }

  #ifdef debug_set_bc_dir_wl
  dmsg << "flag: " << flag;
  dmsg << "ss: " << ss;
  dmsg << "ee: " << ee;
  dmsg << "nsd: " << nsd;
  dmsg << "lDof: " << lDof;
  dmsg << "nNo: " << nNo;
  dmsg << "eNoN: " << eNoN;
  dmsg << "eNoNb: " << eNoNb;
  #endif

  Array<double> tmpA(lDof,nNo), tmpY(lDof,nNo);

  set_bc::set_bc_dir_l(com_mod, lBc, lFa, tmpA, tmpY, lDof);

  if (utils::btest(lBc.bType, iBC_impD)) {
    tmpY = tmpA;
  }

  // Transfer Dirichlet value to global numbering (lFa.nNo -> tnNo).
  // Take effective direction into account if set.
  //
  Array<double> ubg(nsd,tnNo);

  if (std::count(eDir.begin(), eDir.end(), true) != 0) {
    for (int a = 0; a < nNo; a++) {
      int Ac = lFa.gN(a);

      for (int i = 0; i < nsd; i++) {
        int lDof = 0;
        if (eDir[i]) {
          ubg(i,Ac) = tmpY(lDof,a);
          lDof = lDof + 1;
        }
      }
    }

  } else {
    for (int a = 0; a < nNo; a++) {
      int Ac = lFa.gN(a);
      ubg.set_col(Ac, tmpY.col(a));
    }
  }

  Vector<int> ptr(eNoN); 
  Array<double> xl(nsd,eNoN), yl(tDof,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);
  Array<double> xbl(nsd,eNoNb), ubl(nsd,eNoNb);

  // Loop over all the elements of the face and construct residual and
  // tangent matrices
  //

  for (int e = 0; e < nEl; e++) {
    int Ec = lFa.gE(e);
    cDmn = all_fun::domain(com_mod, lM, cEq, Ec);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_fluid) {
      throw std::runtime_error("[set_bc_dir_wl] Weakly applied Dirichlet BC is allowed for fluid phys only");
    }

    // Initialize local residual and stiffness
    lR = 0.0;
    lK = 0.0;

    // Create local copies of fluid velocity and position vector
    //
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,Ec);
      ptr(a) = Ac;

      for (int i = 0; i < tDof; i++) {
        yl(i,a) = Yg(i,Ac);
      }

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        if (com_mod.mvMsh) {
          xl(i,a) = xl(i,a) + Dg(i+nsd+1,Ac);
        }
      }
    }

    // Set function spaces for velocity and pressure on mesh
    std::array<fsType,2> fs;
    fs::get_thood_fs(com_mod, fs, lM, flag, 1);

    Array<double> xwl(nsd,fs[0].eNoN); 
    Vector<double> Nw(fs[0].eNoN); 
    Array<double> Nwx(nsd,fs[0].eNoN), Nwxi(nsd,fs[0].eNoN);

    Array<double> xql(nsd,fs[1].eNoN), Nqx(nsd,fs[1].eNoN);
    Vector<double> Nq(fs[1].eNoN); 

    xwl = xl;

    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < fs[1].eNoN; j++) {
        xql(i,j) = xl(i,j);
      }
    }

    // Create local copies of the wall/solid/interface quantites
    //
    for (int a = 0; a < eNoNb; a++) {
      int Ac = lFa.IEN(a,e);

      for (int i = 0; i < nsd; i++) {
        xbl(i,a) = com_mod.x(i,Ac);
        ubl(i,a) = ubg(i,Ac);

        if (com_mod.mvMsh) {
          // xbl(i,a) = xbl(i,a) + Dg(i+nsd,Ac);
          xbl(i,a) = xbl(i,a) + Dg(i+nsd+1,Ac);
        }
      }
    }

    // Initialize parameteric coordinate for Newton's iterations.
    // Newton method is used to compute derivatives on the face using
    // mesh-based shape functions as an inverse problem.
    //
    Vector<double> xi0(nsd);
    for (int g = 0; g < fs[0].nG; g++) {
      xi0 = xi0 + fs[0].xi.col(g);
    }
    xi0 = xi0 / static_cast<double>(fs[0].nG);

    // Gauss integration 1
    //
    for (int g = 0; g < lFa.nG; g++) {
      Vector<double> nV(nsd);
      auto Nx = lFa.Nx.slice(g);
      nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, eNoNb, Nx, nV);
      double Jac = sqrt(utils::norm(nV));
      nV = nV / Jac;
      double w = lFa.w(g) * Jac;

      Vector<double> ub(nsd), xp(nsd);
      for (int a = 0; a < eNoNb; a++) {
        xp = xp + xbl.col(a)*lFa.N(a,g);
        ub = ub + ubl.col(a)*lFa.N(a,g);
      }

      // Compute Nw and Nwxi of the mesh at the face integration
      // point using Newton method. Then calculate Nwx.
      //
      auto xi = xi0;
      nn::get_nnx(nsd, fs[0].eType, fs[0].eNoN, xwl, fs[0].xib, fs[0].Nb, xp, xi, Nw, Nwxi);

      if (g==0 || !fs[0].lShpF) {
        Array<double> Ks(nsd,nsd);
        nn::gnn(fs[0].eNoN, nsd, nsd, Nwxi, xwl, Nwx, Jac, Ks);
      }

      // Compute Nq of the mesh at the face integration point using
      // Newton method.
      //
      xi = xi0;
      nn::get_nnx(nsd, fs[1].eType, fs[1].eNoN, xql, fs[1].xib, fs[1].Nb, xp, xi, Nq, Nqx);

      if (nsd == 3) {
        fluid::bw_fluid_3d(com_mod, fs[0].eNoN, fs[1].eNoN, w, Nw, Nq, Nwx, yl, ub, nV, tauB, lR, lK);
      } else {
        fluid::bw_fluid_2d(com_mod, fs[0].eNoN, fs[1].eNoN, w, Nw, Nq, Nwx, yl, ub, nV, tauB, lR, lK);
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }
}

/// @brief Set outlet BCs.
//
void set_bc_neu(ComMod& com_mod, const CmMod& cm_mod, const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  #define n_debug_set_bc_neu
  #ifdef debug_set_bc_neu
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  #ifdef debug_set_bc_neu
  dmsg << "cEq: " << cEq;
  #endif

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int iFa = bc.iFa;
    int iM = bc.iM;
    #ifdef debug_set_bc_neu
    dmsg << "----- iBc " << iBc+1;
    #endif

    if (utils::btest(bc.bType, iBC_Neu)) {
      #ifdef debug_set_bc_neu
      dmsg << "iM: " << iM+1;
      dmsg << "iFa: " << iFa+1;
      dmsg << "name: " << com_mod.msh[iM].fa[iFa].name;
      #endif
      set_bc_neu_l(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa], Yg, Dg);

    } else if (utils::btest(bc.bType,iBC_trac)) { 
      set_bc_trac_l(com_mod, cm_mod, bc, com_mod.msh[iM].fa[iFa]); 
    } 
  }
}

/// @brief Set Neumann BC
//
void set_bc_neu_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa, const Array<double>& Yg, const Array<double>& Dg) 
{
  using namespace consts;

  #define n_debug_set_bc_neu_l
  #ifdef debug_set_bc_neu_l
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  int tnNo = com_mod.tnNo;
  int nsd = com_mod.nsd;
  auto& Yn = com_mod.Yn;

  int nNo = lFa.nNo;
  Vector<double> h(1), rtmp(1);
  Vector<double> tmpA(nNo); 

  // Geting the contribution of Neu BC
  //
  if (utils::btest(lBc.bType,iBC_cpl) || utils::btest(lBc.bType,iBC_RCR)) {
    h(0) = lBc.g;

  } else {
    
    if (utils::btest(lBc.bType,iBC_gen)) {
       // [NOTE] The Fortran code was passing a vector to 'igbc' 
       // which treats is as an array dY(gm%dof,SIZE(gm%d,2).
       //
       // So here we create and pass an array and later
       // set the values of the 'tmpA'.
       //
       int nrows = lBc.gm.dof;
       int ncols = lBc.gm.d.ncols();
       Array<double> hg(nrows,ncols);
       Array<double> tmpA_a(nrows,ncols);

       igbc(com_mod, lBc.gm, tmpA_a, hg);

       int n = 0;
       for (int j = 0; j < ncols; j++) {
         for (int i = 0; i < nrows; i++) {
           tmpA(n) = tmpA_a(i,j);
           n += 1;
         }
       }

     } else if (utils::btest(lBc.bType,iBC_res)) {
       h(0) = lBc.r * all_fun::integ(com_mod, cm_mod, lFa, Yn, eq.s, eq.s+nsd-1);

     } else if (utils::btest(lBc.bType,iBC_std)) {
       h(0) = lBc.g;

     } else if (utils::btest(lBc.bType,iBC_ustd)) {
       ifft(com_mod, lBc.gt, h, rtmp);

     } else {
       throw std::runtime_error("[set_bc_neu_l] Correction in SETBCNEU is needed");
     }
  }

  #ifdef debug_set_bc_neu_l
  dmsg << "h(1): " << h(0);
  dmsg << "tnNo: " << tnNo;
  dmsg << "lBc.flwP: " << lBc.flwP;
  #endif

  Vector<double> hg(tnNo);

  // Transforming it to a unified format
  //
  if (utils::btest(lBc.bType,iBC_gen)) {
    for (int a = 0; a < nNo; a++) {
      int Ac = lFa.gN(a);
      hg(Ac) = tmpA(a);
    }
  } else {
    for (int a = 0; a < nNo; a++) {
      int Ac = lFa.gN(a);
      hg(Ac) = -h(0)*lBc.gx(a);
      #ifdef debug_set_bc_neu_l
      dmsg << "hg(Ac): " << hg(Ac);
      #endif
    }
  }

  // Add Neumann BCs contribution to the residual (and tangent if flwP)
  //
  if (lBc.flwP) {
    eq_assem::b_neu_folw_p(com_mod, lBc, lFa, hg, Dg);

  } else {
    eq_assem::b_assem_neu_bc(com_mod, lFa, hg, Yg);
  }


  // Now update surface integrals involved in coupled/resistance BC
  // contribution to stiffness matrix to reflect deformed geometry.
  // We must do this if we have a resistance boundary condition and
  // a follower pressure load (struct/ustruct) or a moving mesh (FSI)
  if (utils::btest(lBc.bType, iBC_res)) {
    if (lBc.flwP || com_mod.mvMsh) {
      eq_assem::fsi_ls_upd(com_mod, lBc, lFa);
    }
  }
  // Now treat Robin BC (stiffness and damping) here
  //
  if (utils::btest(lBc.bType,iBC_Robin)) {
    set_bc_rbnl(com_mod, lFa, lBc.k, lBc.c, lBc.rbnN, Yg, Dg);
  }
}

/// @brief Set Robin BC contribution to residual and tangent
//
void set_bc_rbnl(ComMod& com_mod, const faceType& lFa, const double ks, const double cs, const bool isN, 
  const Array<double>& Yg, const Array<double>& Dg)
{
  using namespace consts;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const double dt = com_mod.dt;
  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  auto& cDmn = com_mod.cDmn;

  int s = eq.s;
  double afv = eq.af * eq.gam * dt;
  double afu = eq.af * eq.beta * dt * dt;
  double afm = afv / eq.am;

  int iM = lFa.iM;
  int eNoN = lFa.eNoN;

  Vector<double> N(eNoN);
  Vector<int> ptr(eNoN);
  Array<double> xl(nsd,eNoN), yl(nsd,eNoN), dl(nsd,eNoN), lR(dof,eNoN); 
  Array3<double> lK(dof*dof,eNoN,eNoN), lKd(nsd*dof,eNoN,eNoN);

  for (int e = 0; e < lFa.nEl; e++) {
    cDmn = all_fun::domain(com_mod, com_mod.msh[iM], cEq, lFa.gE(e));
    auto cPhys = eq.dmn[cDmn].phys;

    for (int a = 0; a < eNoN; a++) {
      int Ac = lFa.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        yl(i,a) = Yg(i+s,Ac);
        dl(i,a) = Dg(i+s,Ac);
      }
    }

    if (lFa.eType == ElementType::NRB) {
    }

    lK  = 0.0;
    lR  = 0.0;
    lKd = 0.0;

    for (int g = 0; g < lFa.nG; g++) {
      Vector<double> nV(nsd);
      auto Nx = lFa.Nx.slice(g);
      nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, eNoN, Nx, nV);
      double Jac = sqrt(utils::norm(nV));
      nV  = nV / Jac;
      double w = lFa.w(g) * Jac; 
      N = lFa.N.col(g);
      Vector<double> u(nsd), ud(nsd);

      for (int a = 0; a < eNoN; a++) {
        for (int i = 0; i < nsd; i++) {
          u(i)  = u(i)  + N(a)*dl(i,a);
          ud(i) = ud(i) + N(a)*yl(i,a);
        }
      }

      auto nDn = mat_fun::mat_id(nsd);
      Vector<double> h;
      h = ks*u + cs*ud;

      if (isN) {
        h = utils::norm(h, nV) * nV;
        for (int a = 0; a < nsd; a++) {
          for (int b = 0; b < nsd; b++) {
            nDn(a,b) = nV(a)*nV(b);
          }
        }
      }

      if (nsd == 3) {
        for (int a = 0; a < eNoN; a++) {
          lR(0,a) = lR(0,a) + w*N(a)*h(0);
          lR(1,a) = lR(1,a) + w*N(a)*h(1);
          lR(2,a) = lR(2,a) + w*N(a)*h(2);
        }

        if (cPhys == EquationType::phys_ustruct) {
          double wl = w * afv;

          for (int a = 0; a < eNoN; a++) {
            for (int b = 0; b < eNoN; b++) {
              double T1 = wl*N(a)*N(b);
              double T2 = (afm*ks + cs)*T1;
              T1 = T1*ks;

              // dM_1/dV_1 + af/am*dM_1/dU_1
              lKd(0,a,b) = lKd(0,a,b) + T1*nDn(0,0);
              lK(0,a,b)  = lK(0,a,b)  + T2*nDn(0,0);

              // dM_1/dV_2 + af/am*dM_1/dU_2
              lKd(1,a,b) = lKd(1,a,b) + T1*nDn(0,1);
              lK(1,a,b)  = lK(1,a,b)  + T2*nDn(0,1);

              // dM_1/dV_3 + af/am*dM_1/dU_3
              lKd(2,a,b) = lKd(2,a,b) + T1*nDn(0,2);
              lK(2,a,b)  = lK(2,a,b)  + T2*nDn(0,2);

              // dM_2/dV_1 + af/am*dM_2/dU_1
              lKd(3,a,b) = lKd(3,a,b) + T1*nDn(1,0);
              lK(4,a,b)  = lK(4,a,b)  + T2*nDn(1,0);

              // dM_2/dV_2 + af/am*dM_2/dU_2
              lKd(4,a,b) = lKd(4,a,b) + T1*nDn(1,1);
              lK(5,a,b)  = lK(5,a,b)  + T2*nDn(1,1);

              // dM_2/dV_3 + af/am*dM_2/dU_3
              lKd(5,a,b) = lKd(5,a,b) + T1*nDn(1,2);
              lK(6,a,b)  = lK(6,a,b)  + T2*nDn(1,2);

              // dM_3/dV_1 + af/am*dM_3/dU_1
              lKd(6,a,b) = lKd(6,a,b) + T1*nDn(2,0);
              lK(8,a,b)  = lK(8,a,b)  + T2*nDn(2,0);

              // dM_3/dV_2 + af/am*dM_3/dU_2
              lKd(7,a,b) = lKd(7,a,b) + T1*nDn(2,1);
              lK(9,a,b) = lK(9,a,b) + T2*nDn(2,1);

              // dM_3/dV_3 + af/am*dM_3/dU_3
              lKd(8,a,b) = lKd(8,a,b) + T1*nDn(2,2);
              lK(10,a,b) = lK(10,a,b) + T2*nDn(2,2);
            }
          }

        // cPhys != EquationType::phys_ustruct
        //
        } else {
          double wl = w * (ks*afu + cs*afv);

          for (int a = 0; a < eNoN; a++) {
            for (int b = 0; b < eNoN; b++) {
              double T1 = N(a)*N(b);
              lK(0,a,b) = lK(0,a,b) + wl*T1*nDn(0,0);
              lK(1,a,b) = lK(1,a,b) + wl*T1*nDn(0,1);
              lK(2,a,b) = lK(2,a,b) + wl*T1*nDn(0,2);

              lK(dof+0,a,b) = lK(dof+0,a,b) + wl*T1*nDn(1,0);
              lK(dof+1,a,b) = lK(dof+1,a,b) + wl*T1*nDn(1,1);
              lK(dof+2,a,b) = lK(dof+2,a,b) + wl*T1*nDn(1,2);

              lK(2*dof+0,a,b) = lK(2*dof+0,a,b) + wl*T1*nDn(2,0);
              lK(2*dof+1,a,b) = lK(2*dof+1,a,b) + wl*T1*nDn(2,1);
              lK(2*dof+2,a,b) = lK(2*dof+2,a,b) + wl*T1*nDn(2,2);
            }
          }
        }

      } else if (nsd == 2) {
        for (int a = 0; a < eNoN; a++) {
          lR(0,a) = lR(0,a) + w*N(a)*h(0);
          lR(1,a) = lR(1,a) + w*N(a)*h(1);
         }

        if (cPhys == EquationType::phys_ustruct) {
          double wl = w*afv;

          for (int a = 0; a < eNoN; a++) {
            for (int b = 0; b < eNoN; b++) {
              double T1 = wl*N(a)*N(b);
              double T2 = (afm*ks + cs)*T1;
              T1 = T1*ks;

              // dM_1/dV_1 + af/am*dM_1/dU_1
              lKd(0,a,b) = lKd(0,a,b) + T1*nDn(0,0);
              lK(0,a,b)  = lK(0,a,b)  + T2*nDn(0,0);

              // dM_1/dV_2 + af/am*dM_1/dU_2
              lKd(1,a,b) = lKd(1,a,b) + T1*nDn(0,1);
              lK(1,a,b)  = lK(1,a,b)  + T2*nDn(0,1);

              // dM_2/dV_1 + af/am*dM_2/dU_1
              lKd(2,a,b) = lKd(2,a,b) + T1*nDn(1,0);
              lK(3,a,b)  = lK(3,a,b)  + T2*nDn(1,0);

              // dM_2/dV_2 + af/am*dM_2/dU_2
              lKd(3,a,b) = lKd(3,a,b) + T1*nDn(1,1);
              lK(4,a,b)  = lK(4,a,b)  + T2*nDn(1,1);
            }
          }

        } else {
          double wl = w * (ks*afu + cs*afv);

          for (int a = 0; a < eNoN; a++) {
            for (int b = 0; b < eNoN; b++) {
              double T1 = N(a)*N(b);
              lK(0,a,b) = lK(0,a,b) + wl*T1*nDn(0,0);
              lK(1,a,b) = lK(1,a,b) + wl*T1*nDn(0,1);
              lK(dof+0,a,b) = lK(dof+0,a,b) + wl*T1*nDn(1,0);
              lK(dof+1,a,b) = lK(dof+1,a,b) + wl*T1*nDn(1,1);
            }
          }
        }
      }
    }

    if (cPhys == EquationType::phys_ustruct) {
      ustruct::ustruct_do_assem(com_mod, eNoN, ptr, lKd, lK, lR);
    } else {
      eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
    }

  } // for int e = 0; e < lFa.nEl; e++

}

/// @brief Set Traction BC
//
void set_bc_trac_l(ComMod& com_mod, const CmMod& cm_mod, const bcType& lBc, const faceType& lFa) 
{
  using namespace consts;

  #define n_debug_set_bc_trac_l 
  #ifdef debug_set_bc_trac_l 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  const int tnNo = com_mod.tnNo;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  int iM = lFa.iM;
  int nNo = lFa.nNo;
  int eNoN = lFa.eNoN;
  #ifdef debug_set_bc_trac_l 
  dmsg << "eNoN: " << eNoN;
  dmsg << "nNo: " << nNo;
  #endif

  // Geting the contribution of traction BC
  //
  Array<double> hg(nsd,tnNo);

  if (utils::btest(lBc.bType,iBC_gen)) {
    Array<double> tmpA(nsd,nNo), hl(nsd,nNo);
    igbc(com_mod, lBc.gm, tmpA, hl);
    for (int a = 0; a < nNo; a++) {
      int Ac = lFa.gN(a);
      hg.set_col(Ac, tmpA.col(a));
    }

  } else if (utils::btest(lBc.bType,iBC_std)) {
    for (int a = 0; a < nNo; a++) {
      int Ac = lFa.gN(a);
      hg.set_col(Ac, lBc.h);
    }

  } else if (utils::btest(lBc.bType,iBC_ustd)) {
     if (lBc.gt.d != nsd) {
       throw std::runtime_error("[set_bc_trac_l]  Traction dof not initialized properly");
     }
     Vector<double> h(nsd);
     Vector<double> tmp(nsd);
     ifft(com_mod, lBc.gt, h, tmp);
     for (int a = 0; a < nNo; a++) {
       int Ac = lFa.gN(a);
       hg.set_col(Ac, h);
      }

  } else {
    throw std::runtime_error("[set_bc_trac_l] Undefined time dependence for traction BC on face '" +  lFa.name + "'.");
  }

  Vector<double> N(eNoN); 
  Vector<int> ptr(eNoN); 
  Array<double> hl(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);

  // Constructing LHS/RHS contribution and assembiling them
  //
  for (int e = 0; e < lFa.nEl; e++) {
    cDmn = all_fun::domain(com_mod, com_mod.msh[iM], cEq, lFa.gE(e));
    auto cPhys = eq.dmn[cDmn].phys;

    if (lFa.eType == ElementType::NRB) {
      //CALL NRBNNXB(msh(iM), lFa, e)
    }

    for (int a = 0; a < eNoN; a++) {
      int Ac = lFa.IEN(a,e);
      ptr(a) = Ac;
      hl.set_col(a, hg.col(Ac));
    }

    lK = 0.0;
    lR = 0.0;

    for (int g = 0; g < lFa.nG; g++) {
      Vector<double> nV(nsd);
      auto Nx = lFa.Nx.slice(g);
      nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, eNoN, Nx, nV);
      double Jac = sqrt(utils::norm(nV));
      double w = lFa.w(g)*Jac;
      N = lFa.N.col(g);

      Vector<double> h(nsd);
      for (int a = 0; a < eNoN; a++) {
        h = h + N(a)*hl.col(a);
      }

      for (int a = 0; a < eNoN; a++) {
        for (int i = 0; i < nsd; i++) {
          lR(i,a) = lR(i,a) - w*N(a)*h(i);
        }
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }
}

/// @brief Treat Neumann boundaries that are not deforming.
///
/// Leave the row corresponding to the master node of the owner
/// process in the LHS matrix and the residual vector untouched. For
/// all the other nodes of the face, set the residual to be 0 for
/// velocity dofs. Zero out all the elements of corresponding rows of
/// the LHS matrix. Make the diagonal elements of the LHS matrix equal
/// to 1 and the column entry corresponding to the master node, -1
//
void set_bc_undef_neu(ComMod& com_mod)
{
  using namespace consts;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];

  for (int iBc = 0; iBc < eq.nBc; iBc++) {
    auto& bc = eq.bc[iBc];
    int iFa = bc.iFa;
    int iM  = bc.iM;

    if (utils::btest(bc.bType,iBC_undefNeu)) {
      set_bc_undef_neu_l(com_mod, bc, com_mod.msh[iM].fa[iFa]);
    }
  }
}

/// Modifies: com_mod.R, com_mod.Val
//
void set_bc_undef_neu_l(ComMod& com_mod, const bcType& lBc, const faceType& lFa)
{
  using namespace consts;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const int nsd = com_mod.nsd;

  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;
  auto& R = com_mod.R;
  auto& Val = com_mod.Val;

  int masN = lBc.masN;

  if (lFa.nNo == 0 || masN == 0) {
    return;
  }

  if (nsd == 2) {
    for (int a = 0; a < lFa.nNo; a++) {
      int rowN = lFa.gN(a);
      if (rowN == masN) {
        continue;
      }

      R(0,rowN) = 0.0;
      R(1,rowN) = 0.0;

      // Diagonalize the stiffness matrix (A)
      //
      for (int i = rowPtr(rowN); i <= rowPtr(rowN+1)-1; i++) {
        int colN = colPtr(i);

        if (colN == rowN) {
          Val(0,i) = 1.0;
          Val(4,i) = 1.0;
        } else if (colN == masN) {
          Val(0,i) = -1.0;
          Val(4,i) = -1.0;
        }
      }
    }

  } else if (nsd == 3) {
    for (int a = 0; a < lFa.nNo; a++) {
      int rowN = lFa.gN(a);
      if (rowN == masN) {
        continue;
      }
      R(0,rowN) = 0.0;
      R(1,rowN) = 0.0;
      R(2,rowN) = 0.0;

      // Diagonalize the stiffness matrix (A)
      //
      for (int i = rowPtr(rowN); i <= rowPtr(rowN+1)-1; i++) {
        int colN = colPtr(i);
        if (colN == rowN) {
          Val(0, i) = 1.0;
          Val(5, i) = 1.0;
          Val(10,i) = 1.0;
        } else if (colN == masN) {
          Val(0, i) = -1.0;
          Val(5, i) = -1.0;
          Val(10,i) = -1.0;
        }
      }
    }
  }
}

};


