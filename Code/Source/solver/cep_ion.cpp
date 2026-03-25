// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "cep_ion.h"

#include "all_fun.h"
#include "post.h"
#include "utils.h"
#include <math.h>

namespace cep_ion {

/// @brief Modifies:
/// \code {.cpp}
///   cep_mod.Xion
/// \endcode
//
void cep_init(Simulation* simulation)
{
  using namespace consts;
  auto& com_mod = simulation->com_mod;

  #define n_debug_cep_init 
  #ifdef debug_cep_init 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cm = com_mod.cm;
  auto& cep_mod = simulation->cep_mod;
  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const int nXion = cep_mod.nXion;
  #ifdef debug_cep_init 
  dmsg << "tnNo: " << tnNo;
  dmsg << "nXion: " << nXion;
  #endif

  for (auto& eq : com_mod.eq) {
    if (eq.phys != EquationType::phys_CEP) {
      continue;
    }

    if (com_mod.dmnId.size() != 0) {
      Vector<double> sA(tnNo); 
      Array<double> sF(nXion,tnNo);

      for (int a = 0; a < tnNo; a++) {
        if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_CEP)) {
          continue;
        }
        for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
          auto cPhys = eq.dmn[iDmn].phys;
          int dID = eq.dmn[iDmn].Id;
          if ((cPhys != EquationType::phys_CEP) || (dID >= 0 && !utils::btest(com_mod.dmnId(a),dID))) {
            continue;
          }
          int nX = eq.dmn[iDmn].cep.nX;
          int nG = eq.dmn[iDmn].cep.nG;
          int imyo = eq.dmn[iDmn].cep.imyo;

          Vector<double> Xl(nX); 
          Vector<double> Xgl(nG);

          cep_init_l(eq.dmn[iDmn].cep, nX, nG, Xl, Xgl);

          sA(a) = sA(a) + 1.0;

          for (int i = 0; i < nX; i++) {
            sF(i,a)  = sF(i,a) + Xl(i);
          }

          for (int i = 0; i < nG; i++) {
            sF(i+nX,a) = sF(i+nX,a) + Xgl(i);
          }
        }
      }

      all_fun::commu(com_mod, sA);
      all_fun::commu(com_mod, sF);

      for (int a = 0; a < tnNo; a++) {
        if (!utils::is_zero(sA(a))) {
          for (int i = 0; i < cep_mod.Xion.nrows(); i++) {
            cep_mod.Xion(i,a) = sF(i,a) / sA(a);
          }
        }
      }

    } else {
      for (int a = 0; a < tnNo; a++) { 
        if (!all_fun::is_domain(com_mod, eq, a, EquationType::phys_CEP)) {
          continue;
        }
        int nX = eq.dmn[0].cep.nX;
        int nG = eq.dmn[0].cep.nG;
        Vector<double> Xl(nX); 
        Vector<double> Xgl(nG);

        cep_init_l(eq.dmn[0].cep, nX, nG, Xl, Xgl);

        for (int i = 0; i < nX; i++) {
          cep_mod.Xion(i,a) = Xl(i);
        }
        for (int i = 0; i < nG; i++) {
          cep_mod.Xion(i+nX,a) = Xgl(i);
        }
      }
    }
  }
}

//------------
// cep_init_l
//------------
//
void cep_init_l(cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg)
{
  switch (cep.cepType) {

    case ElectrophysiologyModelType::AP:
      cep.ap.init(nX, X);
    break;

    case ElectrophysiologyModelType::BO:
      cep.bo.init(nX, X);
    break;

    case ElectrophysiologyModelType::FN:
      cep.fn.init(nX, X);
    break;

    case ElectrophysiologyModelType::TTP:
      cep.ttp.init(cep.imyo, nX, nG, X, Xg);
    break;
  }
}

//-----------
// cep_integ
//-----------
// State variable integration.
//
void cep_integ(Simulation* simulation, const int iEq, const int iDof, const Array<double>& Dg)
{
  static bool IPASS = true;

  using namespace consts;

  auto& com_mod = simulation->com_mod;

  #define n_debug_cep_integ 
  #ifdef debug_cep_integ
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cm = com_mod.cm;
  int tnNo = com_mod.tnNo;
  double dt = com_mod.dt;
  double time = com_mod.time;

  auto& cep_mod = simulation->cep_mod;
  auto& cem = cep_mod.cem;
  auto& eq = com_mod.eq[iEq];

  auto& Yo = com_mod.Yo;
  auto& Xion = cep_mod.Xion;
  int nXion = cep_mod.nXion;

  Vector<double> I4f(tnNo);

  #ifdef debug_cep_integ
  dmsg << "cem.cpld: " << cem.cpld;
  dmsg << "time: " << time;
  #endif

  // Electromechanics: get fiber stretch for stretch activated currents
  //
  if (cem.cpld) {
    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      auto& msh = com_mod.msh[iM];

      if (msh.nFn != 0) {
        Vector<double> sA(msh.nNo);
        post::fib_strech(simulation, iEq, msh, Dg, sA);
        for (int a = 0; a < msh.nNo; a++) {
          int Ac = msh.gN(a);
          I4f(Ac) = sA(a);
        }
      }
    }
  }

  //  Ignore first pass as Xion is already initialized
  if (IPASS) {
    IPASS = false;

  // Copy action potential after diffusion as first state variable
  } else {
    for (int Ac = 0; Ac < tnNo; Ac++) {
      Xion(0,Ac) = Yo(iDof,Ac);
    }
  }

  // Integrate electric potential based on cellular activation model
  //
  if (com_mod.dmnId.size() != 0) {
    Vector<double> sA(tnNo); 
    Array<double> sF(nXion,tnNo); 
    Vector<double> sY(tnNo);

    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!all_fun::is_domain(com_mod, eq, Ac, Equation_CEP)) {
        continue;
      }

      for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
        auto& dmn = eq.dmn[iDmn];
        auto cPhys = dmn.phys;
        int dID = dmn.Id;

        if (cPhys != Equation_CEP || (dID >= 0 && !utils::btest(com_mod.dmnId(Ac),dID))) {
          continue;
	      }

        int nX = dmn.cep.nX;
        int nG = dmn.cep.nG;
        #ifdef debug_cep_integ
        dmsg << "nX: " << nX ;
        dmsg << "nG: " << nG ;
        #endif

        auto Xl = Xion.rows(0,nX-1,Ac);

        // [NOTE] nG can be 0.
        Vector<double> Xgl;
        if (nG != 0) {
          Xgl.resize(nG);
          for (int i = 0; i < nG; i++) {
            Xgl(i) = Xion(i+nX,Ac);
          }
        }

        double yl = 0.0;
        if (cem.cpld) {
          yl = cem.Ya(Ac);
        }

        cep_integ_l(cep_mod, dmn.cep, nX, nG, Xl, Xgl, time-dt, yl, I4f(Ac), dt);

        sA(Ac) = sA(Ac) + 1.0;
        for (int i = 0; i < nX; i++) {
          sF(i,Ac) += Xl(i);
        }

        for (int i = 0; i < nG; i++) {
          sF(nX+i,Ac) += Xgl(i);
        }

        if (cem.cpld) {
          sY(Ac) = sY(Ac) + yl;
        }
      }
    }

    all_fun::commu(com_mod, sA);
    all_fun::commu(com_mod, sF);

    if (cem.cpld) {
      all_fun::commu(com_mod, sY);
    }

    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!utils::is_zero(sA(Ac))) {
        Xion.set_col(Ac, sF.col(Ac) / sA(Ac));
        if (cem.cpld) {
          cem.Ya(Ac) = sY(Ac) / sA(Ac);
        }
      }
    }

  } else {
    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (!all_fun::is_domain(com_mod, eq, Ac, Equation_CEP)) {
        continue;
      }

      int nX = eq.dmn[0].cep.nX;
      int nG = eq.dmn[0].cep.nG;
      auto Xl = Xion.rows(0,nX-1,Ac);
      auto Xgl = Xion.rows(nX,nX+nG-1,Ac);

      double yl = 0.0;
      if (cem.cpld) {
        yl = cem.Ya(Ac);
      }

      cep_integ_l(cep_mod, eq.dmn[0].cep, nX, nG, Xl, Xgl, time-dt, yl, I4f(Ac), dt);

      for (int i = 0; i < nX; i++) {
        Xion(i,Ac) = Xl(i);
      }

      for (int i = 0; i < nG; i++) {
        Xion(nX+i,Ac) = Xgl(i);
      }

      if (cem.cpld) {
        cem.Ya(Ac) = yl;
      }
    }
  }

  for (int Ac = 0; Ac < tnNo; Ac++) {
    Yo(iDof,Ac) = Xion(0,Ac);
  }
}

//-------------
// cep_integ_l
//-------------
// Integrate local electrophysiology variables from t1 to t1+dt. Also
// integrate excitation-activation variables form coupled electro-
// mechanics. The equations are integrated at domain nodes.
//
void cep_integ_l(CepMod& cep_mod, cepModelType& cep, int nX, int nG, Vector<double>& X, Vector<double>& Xg, 
    const double t1, double& yl, const double I4f, const double dt)
{
  using namespace consts;

  #define n_debug_cep_integ_l
  #ifdef debug_cep_integ_l
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cem = cep_mod.cem;

  // Feedback coefficient for stretch-activated-currents
  double Ksac = 0.0;
  if (I4f > 1.0) {
     Ksac = cep.Ksac * (sqrt(I4f) - 1.0);
  } else {
     Ksac = 0.0;
  }

  // Total time steps
  int nt = static_cast<int>(dt/cep.dt);

  // External stimulus duration
  int icl = static_cast<int>(fmax(floor(t1/cep.Istim.CL),0.0));
  double Ts = cep.Istim.Ts + static_cast<double>(icl)*cep.Istim.CL;
  double Te = Ts + cep.Istim.Td;
  double eps = std::numeric_limits<double>::epsilon();

  #ifdef debug_cep_integ_l
  dmsg << "nt: " << nt;
  dmsg << "Ksac: " << Ksac;
  dmsg << "icl: " << icl;
  dmsg << "Ts: " << Ts;
  dmsg << "Te: " << Te;
  dmsg << "cep.cepType: " << cep.cepType;
  dmsg << "cep.odes.tIntTyp: " << cep.odes.tIntType;
  #endif

  switch (cep.cepType) {
    case ElectrophysiologyModelType::AP: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(2);
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep.ap.integ_fe(nX, X, t, cep.dt, Istim, Ksac);
  
            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.ap.actv_strs(X(0), cep.dt, yl, epsX);
            }
          }
        } break; 

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep.ap.integ_rk(nX, X, t, cep.dt, Istim, Ksac);

            //  Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.ap.actv_strs(X(0), cep.dt, yl, epsX);
            }
          }
        } break; 

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep.ap.integ_cn2(nX, X, t, cep.dt, Istim, Ksac, IPAR, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.ap.actv_strs(X(0), cep.dt, yl, epsX);
            }
          }
        } break; 
      } 
    } break; 

    case ElectrophysiologyModelType::BO: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(5);
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0.0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep.bo.integ_fe(cep.imyo, nX, X, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.bo.actv_strs(X(0), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep.bo.actv_strn(X(3), I4f, cep.dt, yl);
            }
          }
        } break;

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep.bo.integ_rk(cep.imyo, nX, X, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.bo.actv_strs(X(0), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep.bo.actv_strn(X(3), I4f, cep.dt, yl);
            }
          }
        } break;

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep.bo.integ_cn2(cep.imyo, nX, X, t, cep.dt, Istim, Ksac, IPAR, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.bo.actv_strs(X(0), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep.bo.actv_strn(X(3), I4f, cep.dt, yl);
            }
          }
        } break;
      } 
    } break; 

    case ElectrophysiologyModelType::FN: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(2);
      IPAR(0) = cep.odes.maxItr;
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep.fn.integ_fe(nX, X, t, cep.dt, Istim);
          }
        } break;

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep.fn.integ_rk(nX, X, t, cep.dt, Istim);
          }
        } break;

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep.fn.integ_cn2(nX, X, t, cep.dt, Istim, IPAR, RPAR);
           }
        } break;
      }
    } break; 

    case ElectrophysiologyModelType::TTP: {
      Vector<int> IPAR(2); 
      Vector<double> RPAR(18);
      IPAR(0) = cep.odes.maxItr;
      IPAR(1) = 0;
      RPAR(0) = cep.odes.absTol;
      RPAR(1) = cep.odes.relTol;

      switch (cep.odes.tIntType) {
        case TimeIntegratioType::FE: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }
            cep.ttp.integ_fe(cep.imyo, nX, nG, X, Xg, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.ttp.actv_strs(X(3), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep.ttp.actv_strn(X(3), I4f, cep.dt, yl);
            }
          }
        } break;

        case TimeIntegratioType::RK4: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep.ttp.integ_rk(cep.imyo, nX, nG, X, Xg, t, cep.dt, Istim, Ksac, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.ttp.actv_strs(X(3), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep.ttp.actv_strn(X(3), I4f, cep.dt, yl);
            }
          }
        } break;

        case TimeIntegratioType::CN2: {
          for (int i = 0; i < nt; i++) {
            double t = t1 + static_cast<double>(i) * cep.dt;
            double Istim;
            if (t >= Ts-eps &&  t <= Te+eps) {
              Istim = cep.Istim.A;
            } else {
              Istim = 0.0;
            }

            cep.ttp.integ_cn2(cep.imyo, nX, nG, X, Xg, t, cep.dt, Istim, Ksac, IPAR, RPAR);

            // Electromechanics excitation-activation
            if (cem.aStress) {
              double epsX;
              cep.ttp.actv_strs(X(3), cep.dt, yl, epsX);
            } else if (cem.aStrain) {
              cep.ttp.actv_strn(X(3), I4f, cep.dt, yl);
            }
          }
        } break;
      }
    } break; 
  } 

  if (isnan(X(0)) ||  isnan(yl)) {
    throw std::runtime_error("[cep_integ_l] A NaN has been computed during time integration of electrophysiology variables.");
  }
}


};
