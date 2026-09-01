// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// This routine contains multiple functions that are generally
// desined to interface with user.

#include "output.h"
#include "utils.h"

#include <cstdio>
#include <iomanip>
#include <limits>
#include <math.h>
#include <sstream>

namespace output {

namespace {
// Field widths of the history table. The rows and the header use the same
// widths, so that they are vertically aligned. A value wider than its field
// shifts the rest of the row rather than being truncated.
constexpr int eq_width        = 2;  // equation symbol
constexpr int time_step_width = 6;  // time step
constexpr int number_width    = 10; // numbers in scientific notation
constexpr int db_width        = 4;  // dB columns, three digits and their sign
constexpr int ls_iter_width   = 5;  // linear solver iterations
constexpr int pct_width       = 4;  // percentage of the time spent in the linear solver

// Width of the time step and nonlinear iteration column: the time step, the
// dash separating it from the two-digit nonlinear iteration, and the 's'
// flagging a time step whose results are written to a file.
constexpr int iter_width = time_step_width + 4;

// Number of digits printed after the decimal point in scientific notation.
constexpr int number_precision = 3;

// Number of characters in a row of the history table: the fields above,
// plus the ten spaces and the four brackets separating them.
constexpr int table_width = eq_width + iter_width + 4 * number_width +
                            2 * db_width + ls_iter_width + pct_width + 14;

// Separator line of the history table.
std::string separator_line() { return std::string(table_width, '-'); }
} // namespace

void output_header(const Simulation *simulation, std::array<double, 3> &timeP) {
  auto& com_mod = simulation->com_mod;
  auto &cm_mod = simulation->cm_mod;

  // Writes to history file and optionally to cout.
  auto& logger = simulation->logger;

  if (com_mod.cm.slv(cm_mod)) {
    return;
  }

  timeP[0] = utils::cput() - timeP[0];
  timeP[1] = 0.0;

  // The dash of the "N-i" title sits above the dash separating the time step
  // from the nonlinear iteration.
  std::ostringstream header;
  header << " " << std::left << std::setw(eq_width) << "Eq"
         << " " << std::setw(iter_width) << "     N-i"
         << " " << std::right << std::setw(number_width) << "T"
         << "  " << std::setw(db_width) << "dB"
         << " " << std::setw(number_width) << "Ri/R1"
         << " " << std::setw(number_width) << "Ri/R0"
         << " " << std::setw(number_width) << "R/Ri"
         << "   " << std::setw(ls_iter_width) << "lsIt"
         << " " << std::setw(db_width) << "dB"
         << " " << std::setw(pct_width) << "%t";

  logger << separator_line() << std::endl;
  logger << header.str() << std::endl;
  if (com_mod.nEq == 1) {
    logger << separator_line() << std::endl;
  }
}

void output_result(const Simulation *simulation, std::array<double, 3> &timeP,
                   const bool save_results, const int iEq) {
  auto &com_mod = simulation->com_mod;
  auto &cm_mod = simulation->cm_mod;
  auto &eq = com_mod.eq[iEq];
  auto cTS = com_mod.cTS;

#ifdef debug_output_result
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
#endif

  // Writes to history file and optionally to cout.
  auto &logger = simulation->logger;

  if (com_mod.cm.slv(cm_mod)) {
    return;
  }

  double tmp = utils::cput();

  if ((com_mod.nEq > 1) && (iEq == 0) && (eq.itr == 1)) {
    logger << separator_line() << std::endl;
  }

  // The time step and the nonlinear iteration, flagged with an 's' when the
  // results of this time step are written to a file.
  const char *save_flag = save_results ? "s" : "";
  std::ostringstream iter;
  iter << std::setw(time_step_width) << cTS << "-" << eq.itr << save_flag;

  timeP[2] = tmp - timeP[0];

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

  // The residuals are bracketed by '!' when the nonlinear residual has grown by
  // more than 20 dB.
  const char *nl_open = (i > 20) ? "!" : "[";
  const char *nl_close = (i > 20) ? "!" : "]";

  double eps = std::numeric_limits<double>::epsilon();

  if (utils::is_zero(timeP[2],timeP[1])) {
    timeP[2] = (1.0 + eps) * timeP[1] + eps;
  }

  // Percent of time in solver?
  //
  double solver_pct = 100.0 * eq.FSILS.RI.callD / (timeP[2] - timeP[1]);
  timeP[1] = timeP[2];
  if (fabs(solver_pct) > 100.0) {
    solver_pct = 100.0;
  }

  // Add a warning if the solution to the linear system did not converge.
  const char *ls_open = eq.FSILS.RI.success ? "[" : "!";
  const char *ls_close = eq.FSILS.RI.success ? "]" : "!";
  std::string convergence_msg;
  if (!eq.FSILS.RI.success) {
    convergence_msg = "  WARNING: The linear system solution has not converged";
  }

  std::ostringstream row;
  row << " " << std::left << std::setw(eq_width) << eq.sym << " "
      << std::setw(iter_width) << iter.str() << " " << std::right
      << std::scientific << std::setprecision(number_precision)
      << std::setw(number_width) << timeP[2] << " " << nl_open
      << std::setw(db_width) << i << " " << std::setw(number_width) << tmp1
      << " " << std::setw(number_width) << tmp << " " << std::setw(number_width)
      << tmp2 << nl_close << " " << ls_open << std::setw(ls_iter_width)
      << eq.FSILS.RI.itr << " " << std::setw(db_width)
      << static_cast<int>(round(eq.FSILS.RI.dB)) << " " << std::setw(pct_width)
      << static_cast<int>(round(solver_pct)) << ls_close;

  logger << row.str() << convergence_msg << std::endl;

  // Print a warning message if the maximum number of nonlinear iterations has been exceeded.
  if (eq.itr > eq.maxItr) {
    auto msg = "[svMultiPhysics] WARNING: The number of nonlinear iterations (" + std::to_string(eq.itr) + 
        ") has exceeded the maximum number set by the value of the Add_equation/Max_iterations parameter in the svMultiPhysics solver input file.";
    if (!eq.FSILS.RI.success) {
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
void write_restart(Simulation* simulation, std::array<double,3>& timeP, const SolutionStates& solutions)
{
  const auto& An = solutions.current.get_acceleration();
  const auto& Yn = solutions.current.get_velocity();
  const auto& Dn = solutions.current.get_displacement();

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
  const bool risFlag = com_mod.risFlag;
  const bool urisFlag = com_mod.urisFlag;
  const auto& stFileName = com_mod.stFileName;

  auto& cplBC = com_mod.cplBC;
  auto& Ad = com_mod.Ad;
  auto& pS0 = com_mod.pS0;
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
          restart_file.write((char*)cem.Ya_f.data(), cem.Ya_f.msize());
          restart_file.write((char*)cem.Ya_s.data(), cem.Ya_s.msize());
          restart_file.write((char*)cem.Ya_n.data(), cem.Ya_n.msize());

        } else if (risFlag) {
          restart_file.write((char*)Ad.data(), Ad.msize());
          write_ris_data(com_mod, restart_file);

        } else if (urisFlag) {
          restart_file.write((char*)Ad.data(), Ad.msize());
          write_uris_data(com_mod, restart_file);

        } else {
          restart_file.write((char*)Ad.data(), Ad.msize());
        }

      // If not sstEq.

      } else {

        if (pstEq) {
          restart_file.write((char*)pS0.data(), pS0.msize());

        } else if (cepEq) {
          restart_file.write((char*)Xion.data(), Xion.msize());
          restart_file.write((char*)cem.Ya_f.data(), cem.Ya_f.msize());
          restart_file.write((char*)cem.Ya_s.data(), cem.Ya_s.msize());
          restart_file.write((char*)cem.Ya_n.data(), cem.Ya_n.msize());

        } else if (risFlag) {
          write_ris_data(com_mod, restart_file);

        } else if (urisFlag) {
          write_uris_data(com_mod, restart_file);
        }
      }

    // If not dFlag.

    } else {
      if (cepEq) {
        restart_file.write((char*)Xion.data(), Xion.msize());

      } else if (risFlag) {
        write_ris_data(com_mod, restart_file);

      } else if (urisFlag) {
        write_uris_data(com_mod, restart_file);
      }
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

void write_ris_data(ComMod& com_mod, std::ofstream& restart_file)
{
  std::vector<char> clsFlagChar(com_mod.ris.clsFlg.size());

  for (int i = 0; i < com_mod.ris.clsFlg.size(); i++) {
    clsFlagChar[i] = com_mod.ris.clsFlg[i];
  }

  restart_file.write(clsFlagChar.data(), clsFlagChar.size()*sizeof(char));
}

void write_uris_data(ComMod& com_mod, std::ofstream& restart_file)
{
  Vector<int> urisCnt(com_mod.nUris);
  std::vector<char> urisClsFlagChar(com_mod.nUris);

  for (int i = 0; i < com_mod.nUris; i++) {
    urisCnt(i) = com_mod.uris[i].cnt;
    urisClsFlagChar[i] = com_mod.uris[i].clsFlg;
  }

  restart_file.write((char*)urisCnt.data(), urisCnt.msize());
  restart_file.write(urisClsFlagChar.data(), urisClsFlagChar.size()*sizeof(char));
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
void write_results(ComMod& com_mod, const std::array<double,3>& timeP, const std::string& fName, const bool sstEq, const SolutionStates& solutions)
{
  const auto& An = solutions.current.get_acceleration();
  const auto& Yn = solutions.current.get_velocity();
  const auto& Dn = solutions.current.get_displacement();

  int cTS = com_mod.cTS;

  auto& stamp = com_mod.stamp;

  FILE* fp = fopen(fName.c_str(), "w");
  for (auto value : stamp) {
    fprintf(fp, " %d ", value);
  }
  fprintf(fp, "\n");

  fclose(fp);
}

};

