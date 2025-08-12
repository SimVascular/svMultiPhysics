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

// This is to calculate average/flux of variables at the boundaries
// and print them into the "B_" files.

#include "txt.h"

#include "all_fun.h"
#include "consts.h"
#include "post.h"
#include "set_bc.h"
#include "utils.h"
#include <math.h>
#include "svZeroD_subroutines.h"

namespace txt_ns {

/// @brief Create a text file for writing boundary integral quantities.
//
void create_boundary_integral_file(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const std::string& file_name)
{
  int fid = 1;
  if (com_mod.cm.slv(cm_mod)) {
    return;
  }

  bool file_exists = false;
  std::string geometry;

  if (FILE *fp = fopen(file_name.c_str(), "r")) {
    file_exists = true;
    fclose(fp);
  }

  // [NOTE] What is this all about?
  if (com_mod.cTS != 0 && file_exists) {
    //CALL TRIMFILE(cTS+3,fName(i))
    return;
  }

  auto fp = fopen(file_name.c_str(), "w");
  fprintf(fp, "# svMultiPhysics boundary integral results file. \n");
  fprintf(fp, "# Quantities represent averaged scalar or flux values over each mesh face.\n");
  fprintf(fp, "#        \n");
  fprintf(fp, "# Format \n");
  fprintf(fp, "# ------ \n");
  fprintf(fp, "# face areas: [list of mesh face areas]  \n");
  fprintf(fp, "# step time [list of mesh face names] \n");
  fprintf(fp, "# [time step] [time] [list of computed values for each mesh face]\n");

  // Write face areas.
  //
  fprintf(fp, "face areas: ");
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    for (int iFa = 0; iFa < msh.nFa; iFa++) {
      auto& fa = msh.fa[iFa];
      fprintf(fp, " %4.3e ", fa.area);
    }
  }
  fprintf(fp, "\n");

  // Write face names.
  //
  fprintf(fp, "step   time  ");

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    for (int iFa = 0; iFa < msh.nFa; iFa++) {
      auto& fa = msh.fa[iFa];
      auto stmp = fa.name;
      fprintf(fp, " %s ", stmp.c_str());
    }
  }

  fprintf(fp, "\n");
  fclose(fp);
}

/// @brief Create a text file for writing volume integral quantities.
//
void create_volume_integral_file(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const std::string& file_name)
{
  int fid = 1;
  if (com_mod.cm.slv(cm_mod)) {
    return;
  }

  bool file_exists = false;
  std::string geometry;

  if (FILE *fp = fopen(file_name.c_str(), "r")) {
    file_exists = true;
    fclose(fp);
  }

  // [NOTE] What is this all about?
  if (com_mod.cTS != 0 && file_exists) {
    //CALL TRIMFILE(cTS+3,fName(i))
    return;
  }

  auto fp = fopen(file_name.c_str(), "w");
  fprintf(fp, "# svMultiPhysics volume integral results file. \n");
  fprintf(fp, "# Quantities represent averaged scalar or flux values over each volume mesh domain.\n");
  fprintf(fp, "#        \n");
  fprintf(fp, "# Format \n");
  fprintf(fp, "# ------ \n");
  fprintf(fp, "# domain volumes: [list of mesh domain volumes]  \n");
  fprintf(fp, "# step time [list of mesh domain names] \n");
  fprintf(fp, "# [time step] [time] [list of computed values for each mesh domain]\n");

  // Write domain volumes.
  //
  fprintf(fp, "domain volumes: ");
  for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
    fprintf(fp, " %4.3e ", lEq.dmn[iDmn].v);
  }
  fprintf(fp, "\n");

  // Write domain names.
  //
  fprintf(fp, "step   time  ");

  for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
    std::string stmp = "domain-" + std::to_string(lEq.dmn[iDmn].Id);
    if (lEq.dmn[iDmn].Id == -1) {
      stmp = "domain";
    }
    fprintf(fp, "%s ", stmp.c_str());
  }

  fprintf(fp, "\n");
  fclose(fp);
}

/// \todo [NOTE] This is not fully implemented.
//
// init_write - if true then this is the start of the simulation and 
//   so create a new file to initialize output.
//
void txt(Simulation* simulation, const bool init_write) 
{
  using namespace consts;
  using namespace utils;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& chnl_mod = simulation->chnl_mod;
  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;
  auto& cplBC = com_mod.cplBC;

  #define n_debug_txt
  #ifdef debug_txt
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "maxNSD: " << maxNSD;
  dmsg << "cplBC.coupled: " << cplBC.coupled;
  #endif

  int fid = 1;

  // Writing data related to cplBC
  //
  if (!com_mod.resetSim) {
    if (cplBC.coupled) {
      bool ltmp = false;

      if (!init_write) {
        if (cplBC.useGenBC) {
          set_bc::genBC_Integ_X(com_mod, cm_mod, "L");

        } else if (cplBC.useSvZeroD){
          svZeroD::calc_svZeroD(com_mod, cm_mod, 'L');
          
        } else {

          for (auto& bc : com_mod.eq[0].bc) {
            if (utils::btest(bc.bType, iBC_RCR)) {
              ltmp = true;
              break;
            }
          }

          set_bc::cplBC_Integ_X(com_mod, cm_mod, ltmp);
        }
      }

      if (com_mod.cm.mas(cm_mod) && !cplBC.useGenBC) {
        if (init_write) {
          if (com_mod.cTS == 0 || !ltmp) {
            //OPEN(fid, FILE=cplBC.saveName)
            //CLOSE(fid, STATUS='DELETE')
          } else {
            //CALL TRIMFILE(cTS,cplBC.saveName)
          }
        } else {
          //OPEN(fid, FILE=cplBC.saveName, POSITION='APPEND')
          //WRITE(fid,'(ES14.6E2)',ADVANCE='NO') cplBC.xp(1)

          for (int i = 0; i < cplBC.nX; i++) {
            // WRITE(fid,'(2(1X,ES14.6E2))',ADVANCE='NO') cplBC.xn(i), cplBC.fa(i).y
          }

          for (int i = 1; i < cplBC.nXp; i++) {
            //WRITE(fid,'(ES14.6E2)',ADVANCE='NO') cplBC.xp(i)
          }

          //WRITE(fid,*)
          //CLOSE(fid)
        }
      }
    }
  } // resetSim

  // Process eq outputs
  //
  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    #ifdef debug_txt
    dmsg << "---------- iEq " << iEq+1;
    #endif
    Array<double> tmpV(maxNSD,tnNo);
    auto& eq = com_mod.eq[iEq];

    for (int iOut = 0; iOut < eq.nOutput; iOut++) {
      auto& output = eq.output[iOut];
      // Don't write it when it doesn't suppose to be written
      #ifdef debug_txt
      dmsg;
      dmsg << "----- iOut " << iOut+1;
      #endif

      // Get options for computing boundary and volume integral quantities.
      // 
      auto output_options = eq.output[iOut].options;
      if (!output_options.boundary_integral && !output_options.volume_integral) { 
        continue;
      } 

      // eq.s - Pointer to start of unknown Yo(:,s:e)
      // eq.output[iOut].l - Length of the outputed variable
      // eq.output[iOut].o - Offset from the first index
      //
      int l = eq.output[iOut].l;
      int s = eq.s + eq.output[iOut].o;  
      int e = s + l - 1;
      #ifdef debug_txt
      dmsg << "l: " << l;
      dmsg << "s: " << s;
      dmsg << "e: " << e;
      #endif

      auto oGrp = eq.output[iOut].grp;
      bool div = true;
      bool pflag = false;
      #ifdef debug_txt
      dmsg << oGrp;
      dmsg << "oGrp: " << oGrp;
      dmsg << "name: " << eq.output[iOut].name;
      #endif

      switch (oGrp) {
        case OutputNameType::outGrp_NA:
          throw std::runtime_error("NA outGrp in TXT");
        break;

        case OutputNameType::outGrp_A:
          for (int j = 0; j < tnNo; j++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,j) = com_mod.An(i+s,j);
            }
          }
        break;

        case OutputNameType::outGrp_Y:
          for (int j = 0; j < tnNo; j++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,j) = com_mod.Yn(i+s,j);
            }
          }

          if (l == 1) {
            pflag = true;
          }
        break;

        case OutputNameType::outGrp_D:
          for (int j = 0; j < tnNo; j++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,j) = com_mod.Dn(i+s,j);
            }
          }
        break;

        case OutputNameType::outGrp_WSS: 
        case OutputNameType::outGrp_vort: 
        case OutputNameType::outGrp_trac:
          post::all_post(simulation, tmpV, com_mod.Yn, com_mod.Dn, oGrp, iEq);
          for (int a = 0; a < tnNo; a++) {
            auto vec = tmpV.col(a, {0,nsd-1});
            tmpV(0,a) = sqrt(norm(vec));
          }
          l = 1;
        break;

        case OutputNameType::outGrp_eFlx: 
        case OutputNameType::outGrp_hFlx: 
        case OutputNameType::outGrp_divV: 
        case OutputNameType::outGrp_J: 
        case OutputNameType::outGrp_mises:
          post::all_post(simulation, tmpV, com_mod.Yn, com_mod.Dn, oGrp, iEq);
        break;

        case OutputNameType::outGrp_absV:
          for (int a = 0; a < tnNo; a++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,a) = com_mod.Yn(i,a) - com_mod.Yn(i+nsd+1,a);
            }
          }
        break;

        case OutputNameType::outGrp_I:
          div = false;
          for (int a = 0; a < tnNo; a++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,a) = 1.0;
            }
          }
        break;

        default:
          throw std::runtime_error("Undefined output '" + std::to_string(static_cast<int>(oGrp)) + "'");
      } 

      // Write boundary and volume integral quantities.
      //
      #ifdef debug_txt
      dmsg << "Write average/flux of variables ...";
      #endif
      const auto& appPath = chnl_mod.appPath;
      std::string file_name = eq.sym + "_" + output.name;
      std::string boundary_file_name;

      if (l == nsd) {
        boundary_file_name = appPath + "B_" + file_name + "_flux.txt";
      } else {
        boundary_file_name = appPath + "B_" + file_name + "_average.txt";
      }

      std::string volume_file_name = appPath + "V_" + file_name + "_flux.txt";

      #ifdef debug_txt
      dmsg << "init_write: " << init_write;
      dmsg << "boundary_file_name: " << boundary_file_name;
      dmsg << "volume_file_name: " << volume_file_name;
      #endif

      if (init_write) {
        if (output_options.boundary_integral) {
          create_boundary_integral_file(com_mod, cm_mod, eq, boundary_file_name);
        }
        if (output_options.volume_integral) {
          create_volume_integral_file(com_mod, cm_mod, eq, volume_file_name);
        }

      } else {
        if (output_options.boundary_integral) {
          write_boundary_integral_data(com_mod, cm_mod, eq, l, boundary_file_name, tmpV, div, pflag);
        }
        if (output_options.volume_integral) {
          write_volume_integral_data(com_mod, cm_mod, eq, l, volume_file_name, tmpV, div, pflag);
        }
      }
    }

    // ECG leads output
    if (com_mod.cm.idcm() == cm_mod.master) {
      auto& cep_mod = simulation->get_cep_mod();

      double time = com_mod.time;
      for (int index = 0; index < cep_mod.ecgleads.num_leads; index++) {
        FILE *fp = fopen(cep_mod.ecgleads.out_files[index].c_str(), "a+");

        /*if (fp == nullptr)
        {
          fp = fopen(cep_mod.ecgleads.out_files[index].c_str(), "wb");
        }*/

        fprintf(fp, "%g,", time);
        fprintf(fp, "%g", cep_mod.ecgleads.pseudo_ECG[index]);
        fprintf(fp, "\n");

        fclose(fp);
      }
    }
  }

  #ifdef debug_txt
  dmsg << "Done txt";
  #endif
}

/// @brief Compute boundary integral for quantities at the current time step 
/// and write them to a text file.
///
/// NOTE: Be carefu of a potential indexing problem here because 'm' is a length and not an index.
//
void write_boundary_integral_data(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const int m, 
    const std::string file_name, const Array<double>& tmpV, const bool div, const bool pFlag)
{
  #define n_debug_write_boundary_integral_data
  #ifdef debug_write_boundary_integral_data
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "file_name: " << file_name;
  dmsg << "m: " << m;
  dmsg << "div: " << div;
  dmsg << "pFlag: " << pFlag;
  #endif

  double time = com_mod.time;
  int time_step = com_mod.cTS;
  #ifdef debug_wtxt
  dmsg << "time_step: " << time_step;
  dmsg << "time: " << time;
  #endif

  int fid = 1;
  int nsd = com_mod.nsd;
  FILE* fp = nullptr;

  if (com_mod.cm.mas(cm_mod)) {
    fp = fopen(file_name.c_str(), "a");
    fprintf(fp, " %d   %.10e ", time_step, time);
  }

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    bool lTH = false; 

    if (msh.nFs == 2) {
      lTH = true;
    }

    for (int iFa = 0; iFa < msh.nFa; iFa++) {
      auto& fa = msh.fa[iFa];
      double tmp = 0.0;

      if (m == 1) {
        if (div) {
          tmp = fa.area;
          tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 0) / tmp;
        } else {
          if (pFlag && lTH) {
            tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 0, true);
          } else {
            tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 0);
          }
        }

      } else if (m == nsd) {
        tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 0, m-1);
      } else {
        throw std::runtime_error("WTXT only accepts 1 and nsd");
      }

      if (com_mod.cm.mas(cm_mod)) {
        fprintf(fp, " %.10e ", tmp);
      }
    }
  }
  if (com_mod.cm.mas(cm_mod)) {
    fprintf(fp, "\n");
    fclose(fp);
  }
}

/// @brief Compute volume integral for quantities at the current time step 
/// and write them to a text file.
///
/// NOTE: Be carefu of a potential indexing problem here because 'm' is a length and not an index.
//
void write_volume_integral_data(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const int m, 
    const std::string file_name, const Array<double>& tmpV, const bool div, const bool pFlag)
{
  #define n_debug_write_volume_integral_data
  #ifdef debug_write_volume_integral_data
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "file_name: " << file_name;
  dmsg << "m: " << m;
  dmsg << "div: " << div;
  dmsg << "pFlag: " << pFlag;
  #endif

  double time = com_mod.time;
  int time_step = com_mod.cTS;
  #ifdef debug_wtxt
  dmsg << "time_step: " << time_step;
  dmsg << "time: " << time;
  #endif

  int fid = 1;
  int nsd = com_mod.nsd;
  FILE* fp = nullptr;

  if (com_mod.cm.mas(cm_mod)) {
    fp = fopen(file_name.c_str(), "a");
    fprintf(fp, " %d   %.10e ", time_step, time);
  }

  for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
    auto& dmn = lEq.dmn[iDmn];
    double tmp = 0.0;

    if (div) {
      tmp = dmn.v;
      tmp = all_fun::integ(com_mod, cm_mod, dmn.Id, tmpV, 0, m-1) / tmp;
    } else {
      tmp = all_fun::integ(com_mod, cm_mod, dmn.Id, tmpV, 0, m-1, pFlag);
    }

    if (com_mod.cm.mas(cm_mod)) {
      fprintf(fp, " %.10e ", tmp);
    }
  }

  if (com_mod.cm.mas(cm_mod)) {
    fprintf(fp, "\n");
    fclose(fp);
  }
}

};
