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

// This routine contains multiple functions that are generally
// desined to interface with user.

#include "output.h"
#include "utils.h"

#include <math.h>

namespace output {

/// @brief Prepares the output of svFSI to the standard output.
///
/// Modifies: timeP
///
/// \todo NOTE: This is not fully implemented.
//
void output_result(Simulation* simulation,  std::array<double,3>& timeP, const int co, const int iEq)
{
  #ifdef debug_output_result
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& eq = com_mod.eq[iEq];
  auto cTS = com_mod.cTS;

  // Writes to history file and optionally to cout.
  auto& logger = simulation->logger;

  if (com_mod.cm.slv(cm_mod)) {
    return;
  }

  int fid = 1;
  double tmp = utils::cput();
  std::string sepLine(69,'-');

  if (co == 1) {
     timeP[0] = tmp - timeP[0];
     timeP[1] = 0.0;
     logger << sepLine << std::endl;
     logger << " Eq     N-i     T       dB  Ri/R1   Ri/R0    R/Ri     lsIt   dB  %t" << std::endl;
     //std::cout << sepLine << std::endl;
     //std::cout << " Eq     N-i     T       dB  Ri/R1   Ri/R0    R/Ri     lsIt   dB  %t" << std::endl;
     if (com_mod.nEq == 1) {
       logger << sepLine << std::endl;
       //std::cout << sepLine << std::endl;
     }
     return;
  }

  if ((com_mod.nEq > 1) && (iEq == 0) && (eq.itr == 1)) {
    logger << sepLine << std::endl;
  }

  std::string c1 = " ";
  std::string c2 = " ";

  if (co == 3) {
    c1 = "s";
  }

  // NS     1-2  3.82E1  [ -62 7.92E-4 7.92E-4 3.60E-4]  [   5  -15  23]
  //-------------------
  //
  timeP[2] = tmp - timeP[0];
  char time_str[20];
  sprintf(time_str, "%4.3e", timeP[2]);
  std::string sOut = " " + eq.sym + " " + std::to_string(cTS) + "-" + std::to_string(eq.itr) + c1 + " " + time_str;

  int i;
  double tmp1 = 1.0;
  double tmp2 = 1.0;

  if (utils::is_zero(eq.iNorm)) {
    tmp  = 1.0;
    tmp1 = 1.0;
    tmp2 = 1.0;
    i = 0;
  } else {
    tmp = eq.FSILS.RI.iNorm / eq.iNorm;
    tmp1 = tmp / eq.pNorm;
    tmp2 = eq.FSILS.RI.fNorm / eq.FSILS.RI.iNorm;
    i = static_cast<int>(20.0*log10(tmp1));
  }

  if (i > 20) {
    c1 = "!"; 
    c2 = "!";
  } else {
    c1 = "["; 
    c2 = "]";
  }

  char norm1_str[20], norm2_str[20], norm3_str[20];
  sprintf(norm1_str, "%4.3e", tmp);
  sprintf(norm2_str, "%4.3e", tmp1);
  sprintf(norm3_str, "%4.3e", tmp2);

  // NS     1-2  3.82E1  [ -62 7.92E-4 7.92E-4 3.60E-4]  [   5  -15  23] 
  //                       ----------------------------
  sOut += "  " + c1 + std::to_string(i) + " " + norm2_str  + " " + norm1_str + " " + norm3_str + c2;
  double eps = std::numeric_limits<double>::epsilon();

  if (utils::is_zero(timeP[2],timeP[1])) {
    timeP[2] = (1.0 + eps) * timeP[1] + eps;
  }

  // Percent of time in solver?
  //
  tmp = 100.0 * eq.FSILS.RI.callD / (timeP[2] - timeP[1]);
  timeP[1] = timeP[2];
  if (fabs(tmp) > 100.0) {
    tmp = 100.0;
  }

  // Add a warning if the solution to the linear system did not converge.
  std::string convergence_msg;
  if (eq.FSILS.RI.suc) {
     c1 = "["; 
     c2 = "]";
  } else {
     c1 = "!";  
     c2 = "!";
     convergence_msg = "  WARNING: The linear system solution has not converged";
  }

  // NS     1-2  3.82E1  [ -62 7.92E-4 7.92E-4 3.60E-4]  [   5  -15  23] 
  //                                                      -------------
  auto db_str = std::to_string(static_cast<int>(round(eq.FSILS.RI.dB)));
  auto calld_str = std::to_string(static_cast<int>(round(tmp)));
  sOut += "  " + c1 + std::to_string(eq.FSILS.RI.itr) + " " + db_str + " " + calld_str + c2;
  sOut += convergence_msg; 

  if (com_mod.nEq > 1) {
    logger << sOut << std::endl;
  } else {
    logger << sOut << std::endl;
  }

  // Print a warning message if the maximum number of nonlinear iterations has been exceeded.
  if (eq.itr > eq.maxItr) {
    auto msg = "[svMultiPhysics] WARNING: The number of nonlinear iterations (" + std::to_string(eq.itr) + 
        ") has exceeded the maximum number set by the value of the Add_equation/Max_iterations parameter in the svMultiPhysics solver input file.";
    if (!eq.FSILS.RI.suc) {
      msg += " This may be due to the failure of the linear system solution to converge.";
    }
    msg += "\n";
    logger << msg << std::endl; 
  }
}

void read_restart_header(ComMod& com_mod, std::array<int,7>& tStamp, double& timeP, std::ifstream& restart_file)
{
  auto& cTS = com_mod.cTS;
  auto& time = com_mod.time;

  restart_file.read((char*)tStamp.data(), sizeof(tStamp));
  restart_file.read((char*)&cTS, sizeof(cTS));
  restart_file.read((char*)&time, sizeof(time));
  restart_file.read((char*)&timeP, sizeof(timeP));

  for (auto& eq : com_mod.eq) {
    restart_file.read((char*)&eq.iNorm, sizeof(eq.iNorm));
  }
}

/// @brief Reproduces the Fortran 'WRITERESTART' subroutine.
//
void write_restart(Simulation* simulation, std::array<double,3>& timeP)
{
  auto& com_mod = simulation->com_mod;
  #define n_debug_write_restart
  #ifdef debug_write_restart
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "timeP: " << timeP[0] << " " << timeP[1] << " " << timeP[2];
  #endif

  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;
  auto& cep_mod = simulation->cep_mod;
  auto const cTS = com_mod.cTS;
  auto const time = com_mod.time;
  auto const stFileRepl = com_mod.stFileRepl;
  auto const recLn = com_mod.recLn;
  auto& stamp = com_mod.stamp;

  const bool ibFlag = com_mod.ibFlag;
  const bool dFlag = com_mod.dFlag;
  const bool sstEq = com_mod.sstEq; 
  const bool pstEq = com_mod.pstEq;
  const bool cepEq = cep_mod.cepEq;
  const auto& stFileName = com_mod.stFileName;

  auto& cplBC = com_mod.cplBC;
  auto& Ad = com_mod.Ad;
  auto& An = com_mod.An;
  auto& Dn = com_mod.Dn;
  auto& pS0 = com_mod.pS0;
  auto& Yn = com_mod.Yn;
  auto& Xion = cep_mod.Xion;
  auto& cem = cep_mod.cem;

  #ifdef debug_write_restart
  dmsg << "stFileName: " << stFileName;
  dmsg << "dFlag: " << dFlag;
  dmsg << "sstEq: " << sstEq;
  dmsg << "pstEq: " << pstEq;
  dmsg << "cepEq: " << cepEq;
  dmsg << "stFileName: " << stFileName;
  dmsg << "stFileRepl: " << stFileRepl;
  #endif 

  int fid = 27;
  int myID = cm.tF(cm_mod);

  auto fName = stFileName + "_last.bin";
  auto tmpS = fName;
  #ifdef debug_write_restart
  dmsg;
  dmsg << "cTS: " << cTS;
  dmsg << "time: " << time;
  dmsg << "recLn: " << recLn;
  dmsg << "myID: " << myID;
  dmsg << "fName: " << fName;
  dmsg << "stamp: ";
  for (auto value : com_mod.stamp) {
    std::cout << value << " ";
  }
  std::cout;
  #endif 

  if (!com_mod.stFileRepl) {
    char fName_num[100];
    if (cTS >= 1000) {
      sprintf(fName_num, "%d", cTS);
    } else {
      sprintf(fName_num, "%03d", cTS);
    }
    fName = stFileName + "_" + fName_num + ".bin";
  }

  // Create the file.
  //
  if (cm.mas(cm_mod)) {
    int np = cm.np();
    std::ofstream restart_file(fName, std::ios::out | std::ios::binary);
    char data{0};
    for (int i = 0; i < np * recLn; i++) {
      //restart_file.write((char*)&data, sizeof(char));
    }
    restart_file.close();
  }

  // This call is to block all processors
  cm.bcast(cm_mod, &fid);

  std::ofstream restart_file(fName, std::ios::out | std::ios::binary | std::ios::in);
  std::streampos write_pos = (myID - 1) * recLn;
  restart_file.seekp(write_pos);

  write_restart_header(com_mod, timeP, restart_file);
  restart_file.write((char*)cplBC.xn.data(), cplBC.xn.msize());
  restart_file.write((char*)Yn.data(), Yn.msize());
  restart_file.write((char*)An.data(), An.msize());

  if (!ibFlag) {
    if (dFlag) {
      restart_file.write((char*)Dn.data(), Dn.msize());

      if (sstEq) {
        if (pstEq) {
          restart_file.write((char*)pS0.data(), pS0.msize());
          restart_file.write((char*)Ad.data(), Ad.msize());
        } else if (cepEq) {
          restart_file.write((char*)Ad.data(), Ad.msize());
          restart_file.write((char*)Xion.data(), Xion.msize());
          restart_file.write((char*)cem.Ya.data(), cem.Ya.msize());
        } else {
          restart_file.write((char*)Ad.data(), Ad.msize());
        }

      } else {
        if (pstEq) {
          restart_file.write((char*)pS0.data(), pS0.msize());
        } else if (cepEq) {
          restart_file.write((char*)Xion.data(), Xion.msize());
          restart_file.write((char*)cem.Ya.data(), cem.Ya.msize());
        }
      }

    } else {
      if (cepEq) {
        restart_file.write((char*)Xion.data(), Xion.msize());
      } else {
        //WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1), eq%iNorm, cplBC%xn, Yn, An
      }
    }

  // ibFlag = true
  //
  // [NOTE] not implemented.
  //
  } else {
    if (dFlag) {
      if (pstEq) {
        //WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1), eq.iNorm, cplBC.xn, Yn, An, Dn, pS0, ib.Yb, ib.Auo, ib.Ubo
      } else {
        //WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1), eq.iNorm, cplBC.xn, Yn, An, Dn, ib.Yb, ib.Auo, ib.Ubo
      }
    } else {
      //WRITE(fid, REC=myID) stamp, cTS, time, CPUT()-timeP(1), eq.iNorm, cplBC.xn, Yn, An, ib.Yb, ib.Auo, ib.Ubo
    }
  }

  restart_file.close();

  // Create a soft link to the bin file for the last time step.
  //
  if (!com_mod.stFileRepl && cm.mas(cm_mod)) {
    std::string cmd = "ln -f " + fName + " " + tmpS;
    std::system(cmd.c_str());
  }
}

void write_restart_header(ComMod& com_mod, std::array<double,3>& timeP, std::ofstream& restart_file)
{
  auto const cTS = com_mod.cTS;
  auto const time = com_mod.time;
  auto& stamp = com_mod.stamp;
  double cpu_time = utils::cput() - timeP[0];

  restart_file.write((char*)stamp.data(), sizeof(stamp));
  restart_file.write((char*)&cTS, sizeof(cTS));
  restart_file.write((char*)&time, sizeof(time));
  restart_file.write((char*)&cpu_time, sizeof(cpu_time));

  for (auto& eq : com_mod.eq) {
    restart_file.write((char*)&eq.iNorm, sizeof(eq.iNorm));
  }
}

/// \todo [NOTE] not fully implemented.
///
/// Reproduces: WRITE(fid, REC=myID) stamp, cTS, time,CPUT()-timeP(1), eq.iNorm, cplBC.xn, Yn, An, Dn
//
void write_results(ComMod& com_mod, const std::array<double,3>& timeP, const std::string& fName, const bool sstEq)
{
  int cTS = com_mod.cTS;

  auto& An = com_mod.An;
  auto& Dn = com_mod.Dn;
  auto& Yn = com_mod.Yn;

  auto& stamp = com_mod.stamp;

  FILE* fp = fopen(fName.c_str(), "w");
  for (auto value : stamp) {
    fprintf(fp, " %d ", value);
  }
  fprintf(fp, "\n");

  fclose(fp);
}

};

