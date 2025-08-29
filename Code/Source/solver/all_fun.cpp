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

#include "all_fun.h"

#include "fsils_api.hpp"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"
#include "consts.h"

#include <bitset>
#include <math.h>

namespace all_fun {

  using namespace consts;

//--------------
// aspect_ratio
//--------------
//
double aspect_ratio(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x)
{
  Array<int> rowM(eNoN,eNoN-1); 
  Array<int> colM(nDim,nDim-1);

  Vector<double> s(eNoN);
  Array<double> Dsub(nDim,nDim); 
  Vector<double> detD(nDim);

  if (nDim == 2) {
    for (int a = 0; a < eNoN; a++) {
      int ap = a + 1;
      if (a == eNoN-1) {
        ap = 0;
      }
      auto x_diff = x.col(a) - x.col(ap);
      s(a) = sqrt( x_diff * x_diff );
    }

  // This is only for tri and tets so for nDim=3 eNoN=4.
  //
  } else if (nDim == 3) {
    Array<int> rowM{ {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3} }; 
    Array<int> colM{ {0, 1}, {1, 2}, {2, 0} };

    for (int a = 0; a < eNoN; a++) {
      for (int b = 0; b < nDim; b++) {
        Dsub = 1.0;

        for (int i = 0; i < eNoN-1; i++) {
          int irow = rowM(a,i);
          for (int j = 0; j < nDim-1; j++) {
            int icol = colM(b,j);
            Dsub(i,j) = x(icol,irow);
          } 
        } 
        detD(b) = mat_fun::mat_det(Dsub,nDim);
      } 
      s(a) = 0.5 * sqrt(detD * detD);
    } 
  }

  return s.max() / s.min();
}

//-------
// commu
//-------
//
void commu(const ComMod& com_mod, Vector<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  if (U.size() != com_mod.lhs.nNo) {
    throw std::runtime_error("COMMU is only specified for vector with size nNo");
  }

  U = mkc(com_mod, U);

  fsi_linear_solver::fsils_commus(com_mod.lhs, U);
  //CALL FSILS_COMMUS(lhs, U)

  mkci(com_mod, U);
  //CALL MKCI(U)
}

//-------
// commu
//-------
//
void commu(const ComMod& com_mod, Array<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  int m = U.nrows();
  if (U.ncols() != com_mod.lhs.nNo) {
    throw std::runtime_error("COMMU is only specified for vector with size nNo");
  }

  U = mkc(com_mod, U);

  fsi_linear_solver::fsils_commuv(com_mod.lhs, m, U);

  mkci(com_mod, U);

  //U = MKC(U)
  //CALL FSILS_COMMUV(lhs, m, U)
  //CALL MKCI(U)

}

/// @brief This function returns the domain that an element of a mesh belongs to
///
/// Reproduces 'FUNCTION DOMAIN(lM, iEq, e)' defined in ALLFUN.f.
//
int domain(const ComMod& com_mod, const mshType& lM, const int iEq, const int e)
{
  int domain_id = -1;
  auto& eq = com_mod.eq[iEq];

  // Domain Id of -1 counts for the entire domain
  //
  for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
    domain_id = iDmn;
    if (eq.dmn[iDmn].Id == -1) {
      return domain_id;
    }
  }

  if (lM.eId.size() == 0) { 
    throw std::runtime_error("eId is not allocated");
  }

  for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
    domain_id = iDmn;
    if (utils::btest(lM.eId[e], eq.dmn[iDmn].Id)) {
      return domain_id;
    }
  }

  return domain_id;
}

/// @brief Find the face ID and mesh ID based on the face name.
//
void find_face(const std::vector<mshType>& mesh_list, const std::string& faceName, int& iM, int& iFa)
{
  iFa = -1;
  iM = -1;

  for (int mi = 0; mi < mesh_list.size(); mi++) { 
    auto& mesh = mesh_list[mi];
    for (int fi = 0; fi < mesh.nFa; fi++) {
      auto& face = mesh.fa[fi];
      if (face.name == faceName) {
        iM = mi;
        iFa = fi;
        return;
      }
    }
  }

  if (iM == -1) {
    throw std::runtime_error("Can't find face named '" + faceName + "' from defined mesh names."); 
  }
}

/// @brief Find the mesh ID based on the mesh name.
//
void find_msh(const std::vector<mshType>& mesh_list, const std::string& mesh_name, int& iM)
{
  iM = -1;
  
  for (int i = 0; i < mesh_list.size(); i++) {
    if (mesh_list[i].name == mesh_name) { 
      iM = i;
      break;
    }
  }
}

/// @brief Reproduces 'FUNCTION GLOBALRV(lM, U)' defined in ALLFUN.f.
//
Array<double> 
global(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, const Array<double>& U)
{
  auto& cm = com_mod.cm;

  #define n_debug_global_rv
  #ifdef debug_global_rv 
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  int m = U.nrows();
  #ifdef debug_global_rv 
  dmsg << "m: " << m;
  dmsg << "U.ncols(): " << U.ncols();
  dmsg << "lM.nNo: " << lM.nNo;
  dmsg << "lM.gnNo: " << lM.gnNo;
  dmsg << "lM.nEl: " << lM.nEl;
  #endif 

  if (U.ncols() != lM.nNo) {
    throw std::runtime_error("GLOBAL is only specified for array with columns size nNo");
  }

  if (cm.seq()) {
    return U;
  }

  Array<double> result;
  Vector<int> sCount(cm.np()); 
  Vector<int> disp(cm.np());
  Array<double> ienU(m*lM.eNoN, lM.nEl);
  Array<double> gienU;

  if (cm.mas(cm_mod)) {
    gienU.resize(m*lM.eNoN, lM.gnEl); 
    result.resize(m,lM.gnNo);
   } else {
  }

  for (int e = 0; e < lM.nEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      Ac = lM.lN(Ac);
      for (int i = 0; i < m; i++) {
        ienU(m*a+i,e) = U(i,Ac);
      }
    }
  }

  int a = lM.eNoN*m;

  for (int i = 0; i < cm.np(); i++) {
    disp(i) = lM.eDist(i)*a;
    sCount(i) = lM.eDist(i+1)*a - disp(i);
    #ifdef debug_global_rv 
    dmsg << ">>> i: " << i;
    dmsg << "  disp(i): " << disp(i);
    dmsg << "  sCount(i): " << sCount(i);
    #endif
  }

  MPI_Gatherv(ienU.data(), lM.nEl*a, cm_mod::mpreal, gienU.data(), sCount.data(), disp.data(), cm_mod::mpreal, cm_mod.master, cm.com());

  // If a slave process return an empty result.
  if (cm.slv(cm_mod)) {
    return result;
  }

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      for (int i = 0; i < m; i++) {
        result(i,Ac) = gienU(m*a+i,e);
      }
    }
  }

  return result;
}

///////////////////////////////////////////////////////////////////////////
/// Convert the vIntegM function from svFSI to c++.
/// @brief This routine integrated a scalar field over a particular domain.
///
/// @param iM mesh index
/// @param s an array containing a scalar value for each node in the mesh
/// Replicates 'FUNCTION vIntegM(dId, s, l, u, pFlag)' defined in ALLFUN.f.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, int iM, const Array<double>& s)
{
  using namespace consts;

  #define n_debug_integ_v
  #ifdef debug_integ_v
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vIntegM " << " ";
  dmsg << "iM: " << iM;
  #endif

  int nNo = s.ncols();
  int tnNo = com_mod.tnNo;
  bool ibFlag = com_mod.ibFlag;

  if (nNo != tnNo) {
    std::string msg = "Incompatible vector size in vIntegM in domain: ";
    msg += std::to_string(iM);
    throw std::runtime_error(msg);
  }

  double result = 0.0; 
  int nsd = com_mod.nsd;
  fsType fs;

  auto& msh = com_mod.msh[iM];
  int insd = nsd;
  if (msh.lShl) insd = nsd-1;
  if (msh.lFib) insd = 1;

  // Update pressure function space for Taylor-Hood type element
  if (com_mod.msh[iM].nFs == 2) {
    fs.nG    = msh.fs[1].nG;
    fs.eType = msh.fs[1].eType;
    fs.lShpF = msh.fs[1].lShpF;
    fs.eNoN  = msh.fs[1].eNoN;

    fs.w.resize(fs.nG); 
    fs.N.resize(fs.eNoN,fs.nG); 
    fs.Nx.resize(nsd,fs.eNoN,fs.nG);

    if (fs.eType != ElementType::NRB) {
      fs.w  = msh.fs[1].w;
      fs.N  = msh.fs[1].N;
      fs.Nx = msh.fs[1].Nx;
    }

  } else { 
    fs.nG    = msh.fs[0].nG;
    fs.eType = msh.fs[0].eType;
    fs.lShpF = msh.fs[0].lShpF;
    fs.eNoN  = msh.fs[0].eNoN;

    fs.w.resize(fs.nG); 
    fs.N.resize(fs.eNoN,fs.nG); 
    fs.Nx.resize(nsd,fs.eNoN,fs.nG);

    if (fs.eType != ElementType::NRB) {
      fs.w  = msh.fs[0].w;
      fs.N  = msh.fs[0].N;
      fs.Nx = msh.fs[0].Nx;
    }
  }
  int eNoN = fs.eNoN;

  Array<double> xl(nsd,eNoN); 
  Array<double> Nxi(insd,eNoN); 
  Array<double> Nx(insd,eNoN); 
  Vector<double> sl(eNoN); 
  Array<double> tmps(nsd,insd);
  Array<double> tmp(nsd,nsd);

  for (int e = 0; e < msh.nEl; e++) {
    // Updating the shape functions, if this is a NURB
    //
    // [TODO:DaveP] not implemented. 
    //
    if (msh.eType == ElementType::NRB) {
      //CALL NRBNNX(msh(iM), e)
      fs.w  = msh.w;
      fs.N  = msh.N;
      fs.Nx = msh.Nx;
    }

    int ibl = 0;
    for (int a = 0; a < eNoN; a++) { 
      int Ac = msh.IEN(a,e);
      xl.set_col(a, com_mod.x.col(Ac));

      if (com_mod.mvMsh) {
        for (int i = 0; i < nsd; i++) { 
          xl(i,a) += com_mod.Do(i+nsd+1,Ac);
        }
      }
      sl(a) = s(0,Ac);
      ibl = ibl + com_mod.iblank(Ac);
    }

    if (ibl == eNoN) {
      continue;
    }

    double Jac = 0.0;

    for (int g = 0; g < fs.nG; g++) {
      Nxi = fs.Nx.slice(g);

      if (g == 0 || !fs.lShpF) {
        if (msh.lShl) {
          Vector<double> nV(nsd);
          nn::gnns(nsd, eNoN, Nxi, xl, nV, tmps, tmps);
          Jac  = sqrt(utils::norm(nV));
        } else { 
          nn::gnn(eNoN, nsd, insd, Nxi, xl, Nx, Jac, tmp);
        }
      }

      if (utils::is_zero(Jac)) {
        throw std::runtime_error("Jac < 0 for element: " + std::to_string(e) + ")");
      }
      double sHat = 0.0;

      for (int a = 0; a < eNoN; a++) {
        int Ac = msh.IEN(a,e);
        sHat = sHat + sl(a)*fs.N(a,g);
      }
      result += fs.w(g) * Jac * sHat;
    }
  }

  if (com_mod.cm.seq()) {
    return result;
  }
  result = com_mod.cm.reduce(cm_mod, result);
  return result;
}
///////////////////////////////////////////////////////////////////////////

/// @brief This routine integrated a scalar field over a particular domain.
///
/// Note that 'l' and 'u' should be 0-based and are used to index into 's'.
/// @param dId domain id
/// @param s an array containing a scalar value for each node in the mesh
/// @param l lower index of s
/// @param u upper index of s (must be equal to l)
/// @param pFlag flag for using Taylor-Hood function space for pressure
/// Replicates 'FUNCTION vInteg(dId, s, l, u, pFlag)' defined in ALLFUN.f.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, int dId, const Array<double>& s, int l, int u, bool pFlag)
{
  using namespace consts;

  #define n_debug_integ_v
  #ifdef debug_integ_v
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "vInteg " << " ";
  dmsg << "dId: " << dId;
  dmsg << "l: " << l;
  dmsg << "pFlag: " << pFlag;
  #endif

  int nNo = s.ncols();
  int tnNo = com_mod.tnNo;
  bool ibFlag = com_mod.ibFlag;

  if (nNo != tnNo) {
    if (ibFlag) {
      if (nNo != com_mod.ib.tnNo) {
          std::string msg = "Incompatible vector size in vInteg in domain: ";
          msg += std::to_string(dId);
          msg += "\nNumber of nodes in s must be equal to total number of nodes (immersed boundary).\n";
          throw std::runtime_error(msg);
      }
    } else { 
      std::string msg = "Incompatible vector size in vInteg in domain: ";
      msg += std::to_string(dId);
      msg += "\nNumber of nodes in s must be equal to total number of nodes.\n";
      throw std::runtime_error(msg);
    } 
  }

  #ifdef debug_integ_v
  dmsg << "tnNo: " << tnNo;
  dmsg << "nNo: " << nNo;
  #endif

  bool flag = pFlag; 
  if (l != u) {
    std::string msg = " l=" + std::to_string(l) + " u=" + std::to_string(u);
    throw std::runtime_error("Incompatible spatial output setting and element type." + msg);
  }

  bool isIB = false;
  if (ibFlag) {
    if (nNo == com_mod.ib.tnNo) {
      isIB = true;
    }
  } 

  double result = 0.0; 
  int nsd = com_mod.nsd;
  auto& ib = com_mod.ib;
  fsType fs;

  if (!isIB) {
    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      auto& msh = com_mod.msh[iM];
      int insd = nsd;
      if (msh.lShl) insd = nsd-1;
      if (msh.lFib) insd = 1;

      // Update pressure function space for Taylor-Hood type element

      if (flag && (com_mod.msh[iM].nFs == 2)) {
        fs.nG    = msh.fs[1].nG;
        fs.eType = msh.fs[1].eType;
        fs.lShpF = msh.fs[1].lShpF;
        fs.eNoN  = msh.fs[1].eNoN;

        fs.w.resize(fs.nG); 
        fs.N.resize(fs.eNoN,fs.nG); 
        fs.Nx.resize(nsd,fs.eNoN,fs.nG);

        if (fs.eType != ElementType::NRB) {
          fs.w  = msh.fs[1].w;
          fs.N  = msh.fs[1].N;
          fs.Nx = msh.fs[1].Nx;
        }

      } else { 
        fs.nG    = msh.fs[0].nG;
        fs.eType = msh.fs[0].eType;
        fs.lShpF = msh.fs[0].lShpF;
        fs.eNoN  = msh.fs[0].eNoN;

        fs.w.resize(fs.nG); 
        fs.N.resize(fs.eNoN,fs.nG); 
        fs.Nx.resize(nsd,fs.eNoN,fs.nG);

        if (fs.eType != ElementType::NRB) {
          fs.w  = msh.fs[0].w;
          fs.N  = msh.fs[0].N;
          fs.Nx = msh.fs[0].Nx;
        }
      }
      int eNoN = fs.eNoN;

      Array<double> xl(nsd,eNoN); 
      Array<double> Nxi(insd,eNoN); 
      Array<double> Nx(insd,eNoN); 
      Vector<double> sl(eNoN); 
      Array<double> tmps(nsd,insd);
      Array<double> tmp(nsd,nsd);

      for (int e = 0; e < msh.nEl; e++) {
        if (dId > 0 &&  msh.eId.size() != 0) {
          if (!utils::btest(msh.eId(e),dId)) {
            continue;
          }
        }

        // Updating the shape functions, if this is a NURB
        //
        // [TODO:DaveP] not implemented. 
        //
        if (msh.eType == ElementType::NRB) {
          //CALL NRBNNX(msh(iM), e)
          fs.w  = msh.w;
          fs.N  = msh.N;
          fs.Nx = msh.Nx;
        }

        int ibl = 0;
        for (int a = 0; a < eNoN; a++) { 
          int Ac = msh.IEN(a,e);
          xl.set_col(a, com_mod.x.col(Ac));

          if (com_mod.mvMsh) {
            for (int i = 0; i < nsd; i++) { 
              xl(i,a) += com_mod.Do(i+nsd+1,Ac);
            }
          }

          if (l == u) {
            sl(a) = s(l,Ac);
          } else { 
            auto rows = s.col(Ac, {l,u});
            sl(a) = sqrt(utils::norm(rows));
          }
          ibl = ibl + com_mod.iblank(Ac);
        }

        if (ibl == eNoN) {
          continue;
        }

        double Jac = 0.0;

        for (int g = 0; g < fs.nG; g++) {
          Nxi = fs.Nx.slice(g);

          if (g == 0 || !fs.lShpF) {
            if (msh.lShl) {
              Vector<double> nV(nsd);
              nn::gnns(nsd, eNoN, Nxi, xl, nV, tmps, tmps);
              Jac  = sqrt(utils::norm(nV));
            } else { 
              nn::gnn(eNoN, nsd, insd, Nxi, xl, Nx, Jac, tmp);
            }
          }

          if (utils::is_zero(Jac)) {
            throw std::runtime_error("Jac < 0 for element: " + std::to_string(e) + ")");
          }
          double sHat = 0.0;

          for (int a = 0; a < eNoN; a++) {
            int Ac = msh.IEN(a,e);
            sHat = sHat + sl(a)*fs.N(a,g);
          }
          result += fs.w(g) * Jac * sHat;
        }
      }
    }

  } else { 

    for (int iM = 0; iM < ib.nMsh; iM++) {
      int eNoN = ib.msh[iM].eNoN;
      int insd = nsd;

      Array<double> xl(nsd,eNoN); 
      Array<double> Nxi(insd,eNoN); 
      Array<double> Nx(insd,eNoN); 
      Vector<double> sl(eNoN); 
      Array<double> tmps(nsd,insd);
      Array<double> tmp(nsd,nsd);

      for (int e = 0; e < ib.msh[iM].nEl; e++) {
        if (dId > 0 && ib.msh[iM].eId.size() != 0) {
          if (!utils::btest(ib.msh[iM].eId(e),dId)) {
            continue;
          }
        }

        // Updating the shape functions, if this is a NURB
        // [TODO:DaveP] not implemented. 
        if (ib.msh[iM].eType == ElementType::NRB) {
          //CALL NRBNNX(ib.msh(iM), e)
        }

        for (int a = 0; a < eNoN; a++) {
          int Ac = ib.msh[iM].IEN(a,e);
          for (int i = 0; i < nsd; i++) { 
            xl(i,a) = ib.x(i,Ac) + ib.Ubo(i,Ac);
          }

          if (l == u) {
            sl(a) = s(l,Ac);
          } else { 
            auto rows = s.col(Ac, {l,u});
            sl(a) = sqrt(utils::norm(rows));
          }
        }

        for (int g = 0; g < ib.msh[iM].nG; g++) {
          double Jac = 0.0;
          Nxi = ib.msh[iM].Nx.slice(g);
          if (g == 0 ||  !ib.msh[iM].lShpF) {
            nn::gnn(eNoN, nsd, insd, Nxi, xl, Nx, Jac, tmp);
          }
          double sHat = 0.0;

          for (int a = 0; a < eNoN; a++) {
            int Ac = ib.msh[iM].IEN(a,e);
            sHat = sHat + sl(a)*ib.msh[iM].N(a,g);
          }
          result += ib.msh[iM].w(g)*Jac*sHat;
        }
      }
    }
  }

  if (com_mod.cm.seq() || isIB) {
    return result;
  }

  result = com_mod.cm.reduce(cm_mod, result);
  return result;
}

/// @brief This routine integrate a scalar field s over the face lFa.
///
/// Reproduces 'FUNCTION IntegS(lFa, s, pflag)'.
///
/// @param lFa face type, representing a face on the computational mesh
/// @param s an array containing a scalar value for each node in the mesh
/// @param pFlag flag for using Taylor-Hood function space for pressure
/// @param cfg denotes which mechanical configuration (reference/timestep 0, old/timestep n, or new/timestep n+1). Default reference.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Vector<double>& s, bool pFlag, MechanicalConfigurationType cfg)
{
  using namespace consts;
  #define n_debug_integ_s
  #ifdef debug_integ_s
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "IntegS " << " ";
  dmsg << "lFa.iM: " << lFa.iM+1;
  dmsg << "lFa.name: " << lFa.name;
  dmsg << "lFa.eType: " << lFa.eType;
  #endif 

  bool flag = pFlag; 
  int nsd = com_mod.nsd;
  int insd = nsd - 1;

  if (com_mod.msh[lFa.iM].lShl) {
    insd = insd - 1;
  }

  if (com_mod.msh[lFa.iM].lFib) {
    insd = 0;
  }

  int nNo = s.size(); // Total number of nodes on a processor
  #ifdef debug_integ_s
  dmsg << "nNo: " << nNo;
  dmsg << "insd: " << insd;
  dmsg << "flag: " << flag;
  #endif

  if (nNo != com_mod.tnNo) {
    if (com_mod.ibFlag) {
      if (nNo != com_mod.ib.tnNo) {
        std::string msg = "Incompatible vector size in integS on face: ";
        msg += lFa.name;
        msg +=  "\nNumber of nodes in s must be equal to total number of nodes (immersed boundary).\n";
        throw std::runtime_error(msg);
      }
    } else {
      std::string msg = "Incompatible vector size in integS on face: ";
      msg += lFa.name;
      msg +=  "\nNumber of nodes in s must be equal to total number of nodes.\n";
      throw std::runtime_error(msg);
    }
  }

  bool isIB = false;

  if (com_mod.ibFlag) {
    if (nNo == com_mod.ib.tnNo) {
      isIB = true;
    }
  }

  // Update pressure function space for Taylor-Hood element
  //
  fsType fs;

  if (flag) {
    if (lFa.nFs != 2) {
      throw std::runtime_error("Incompatible boundary integral function call and face element type");
    }

    fs.nG    = lFa.fs[1].nG;
    fs.eType = lFa.fs[1].eType;
    fs.lShpF = lFa.fs[1].lShpF;
    fs.eNoN  = lFa.fs[1].eNoN;

    fs.w.resize(fs.nG); 
    fs.N.resize(fs.eNoN,fs.nG); 
    fs.Nx.resize(insd,fs.eNoN,fs.nG);

    if (fs.eType != ElementType::NRB) {
      fs.w  = lFa.fs[1].w;
      fs.N  = lFa.fs[1].N;
      fs.Nx = lFa.fs[1].Nx;
    }
  } else {
    fs.nG    = lFa.fs[0].nG;
    fs.eType = lFa.fs[0].eType;
    fs.lShpF = lFa.fs[0].lShpF;
    fs.eNoN  = lFa.fs[0].eNoN;

    fs.w.resize(fs.nG); 
    fs.N.resize(fs.eNoN,fs.nG); 
    fs.Nx.resize(insd,fs.eNoN,fs.nG);

    if (fs.eType != ElementType::NRB) {
      fs.w  = lFa.fs[0].w;
      fs.N  = lFa.fs[0].N;
      fs.Nx = lFa.fs[0].Nx;
    }
  }

  #ifdef debug_integ_s
  dmsg << "fs.nG: " << fs.nG;
  dmsg << "fs.eType: " << fs.eType;
  dmsg << "fs.w: " << fs.w;
  #endif

  // Initialize integral to 0
  double result = 0.0;

  // Loop over elements on face
  for (int e = 0; e < lFa.nEl; e++) {
    // [TODO:DaveP] not implemented.
    if (lFa.eType == ElementType::NRB) {
      if (!isIB) {
        //CALL NRBNNXB(msh(lFa.iM), lFa, e)
      } else {
        //CALL NRBNNXB(ib.msh(lFa.iM), lFa, e)
      }
      fs.w  = lFa.w;
      fs.N  = lFa.N;
      fs.Nx = lFa.Nx;
    }

    // Loop over the Gauss points
    for (int g = 0; g < fs.nG; g++) {
      Vector<double> n(nsd);
      if (!isIB) {
        // Get normal vector in cfg configuration
        auto Nx = fs.Nx.slice(g);
        nn::gnnb(com_mod, lFa, e, g, nsd, insd, fs.eNoN, Nx, n, cfg);
      }

      // Calculating the Jacobian (encodes area of face element)
      double Jac = sqrt(utils::norm(n));

      // Calculating the function value at Gauss point
      double sHat = 0.0;
      for (int a = 0; a < fs.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        sHat = sHat + s(Ac)*fs.N(a,g);
      }

      // Now integrating
      result = result + Jac*fs.w(g)*sHat;
     }
  }

  // If using multiple processors, add result from all processors
  if (com_mod.cm.seq() || isIB) {
    return result; 
  }

  result = com_mod.cm.reduce(cm_mod, result);
  return result; 
}

/// @brief This routine integrates vector field s dotted with the face normal n
/// over the face lFa. For example, if s contains the velocity at each node on 
/// the face, this function computed the velocity flux through the face.
///
/// Reproduces 'FUNCTION IntegV(lFa, s)'
///
/// @param lFa face type, representing a face on the computational mesh
/// @param s an array containing a vector value for each node in the mesh
/// @param pFlag flag for using Taylor-Hood function space for pressure
/// @param cfg denotes which configuration (reference/timestep 0, old/timestep n, or new/timestep n+1). Default reference.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, 
            const Array<double>& s, MechanicalConfigurationType cfg)
{
  using namespace consts;

  #define n_debug_integ_V
  #ifdef debug_integ_V
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "IntegV " << " ";
  dmsg << "lFa.name: " << lFa.name;
  #endif

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int insd = nsd - 1;
  int tnNo = com_mod.tnNo;

  #ifdef debug_integ_V
  dmsg << "s.nrows(): " << s.nrows();
  dmsg << "nsd: " << nsd;
  #endif

  if (s.nrows() != nsd) {
    std::string msg = "Incompatible vector size in integV on face: "; 
    msg += lFa.name;
    msg += "\nNumber of rows in s must be equal to number of spatial dimensions.\n";
    throw std::runtime_error(msg);
  }

  int nNo = s.ncols();
  #ifdef debug_integ_V
  dmsg << "nNo: " << nNo;
  dmsg << "nsd: " << nsd;
  #endif

  if (nNo != tnNo) {
    if (com_mod.ibFlag) {
      if (nNo != com_mod.ib.tnNo) {
        std::string msg = "Incompatible vector size in integV on face: ";
        msg += lFa.name;
        msg += "\nNumber of nodes in s must be equal to total number of nodes (immersed boundary).\n";
        throw std::runtime_error(msg);
      }
    } else {
      std::string msg = "Incompatible vector size in integV on face: ";
      msg += lFa.name; 
      msg += "\nNumber of nodes in s must be equal to total number of nodes.\n";
      throw std::runtime_error(msg);
    }
  }

  // If using Immersed Boundary Method
  bool isIB = false;
  if (com_mod.ibFlag) {
    if (nNo ==  com_mod.ib.tnNo) {
      isIB = true;
    }
  }

  // Initialize integral to 0
  double result =  0.0;

  // Loop over elements on face
  for (int e = 0; e < lFa.nEl; e++) {
    //dmsg << "----- e " << e+1 << " -----";
    //  Updating the shape functions, if this is a NURB
    if (lFa.eType == ElementType::NRB) {
      if (!isIB) {
         //CALL NRBNNXB(msh(lFa.iM), lFa, e)
      } else {
         //CALL NRBNNXB(ib.msh(lFa.iM), lFa, e)
      }
    }

    // Loop over the Gauss points
    for (int g = 0; g < lFa.nG; g++) {
      //dmsg << ">>> g: " << g+1;
      Vector<double> n(nsd);
      if (!isIB) {
        // Get normal vector in cfg configuration
        auto Nx = lFa.Nx.slice(g);
        nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, lFa.eNoN, Nx, n, cfg);
        //CALL GNNB(lFa, e, g, nsd-1, lFa.eNoN, lFa.Nx(:,:,g), n)
      } else {
        //CALL GNNIB(lFa, e, g, n)
      }

      //  Calculating the function value (s dot n)dA at this Gauss point
      double sHat = 0.0;

      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        // Compute s dot n
        for (int i = 0; i < nsd; i++) {
          sHat = sHat + lFa.N(a,g) * s(i,Ac) * n(i);
          //dmsg << "s(i,Ac): " << s(i,Ac);
        }
      }

      //  Now integrating
      //dmsg << "sHat: " << sHat;
      //dmsg << "lFa.w(g): " << lFa.w(g);
      result = result + lFa.w(g) * sHat;
    }
  }

  // If using multiple processors, add result from all processors
  if (cm.seq() || isIB) {
    return result; 
  }

  result = cm.reduce(cm_mod, result);

  return result; 
}

/// @brief This routine integrate s(l:u,:) over the surface faId, where s is an
/// array of scalars or an array of nsd-vectors. This routine calls overloaded
/// functions to integrate scalars, if s is scalar (i.e. l=u), or vectors if s 
/// is vector (i.e. l<u).
///
/// Note that 'l' and 'u' should be 0-based and are used to index into 's'.
///
/// Reproduces 'FUNCTION IntegG(lFa, s, l, uo, THflag)'.
///
/// @param lFa face type, representing a face on the computational mesh.
/// @param s an array containing a vector value for each node in the mesh.
/// @param l lower index of s
/// @param uo optional: upper index of s. Default u = l.
/// @param THlag flag for using Taylor-Hood function space for pressure.
/// @param cfg denotes which configuration (reference/timestep 0, old/timestep n, or new/timestep n+1). Default reference.
///
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, 
    const Array<double>& s, const int l, std::optional<int> uo, bool THflag, MechanicalConfigurationType cfg)
{
  using namespace consts;

  #define n_debug_integ_g
  #ifdef debug_integ_g
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "IntegG " << " ";
  dmsg << "l: " << l;
  dmsg << "uo: " << uo;
  #endif

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int insd = nsd - 1;
  int tnNo = com_mod.tnNo;

  // Set u if uo is given. Else, set u = l.
  int u{0};
  if (uo.has_value()) { 
    u = uo.value();
  } else {
    u = l;
  }

  bool flag = THflag; 
  int nNo = s.ncols();

  if (nNo != tnNo) {
    if (com_mod.ibFlag) {
      if (nNo != com_mod.tnNo) { 
        std::string msg = "Incompatible vector size in integG on face: ";
        msg += lFa.name; 
        msg += "\nNumber of nodes in s must be equal to total number of nodes (immersed boundary).\n";
        throw std::runtime_error(msg);
      }
    } else {
      std::string msg = "Incompatible vector size in integG on face: ";
      msg += lFa.name;
      msg += "\nNumber of nodes in s must be equal to total number of nodes.\n";
      throw std::runtime_error(msg);
    }
  }

  // Initialize integral to 0 (not strictly necessary to initialize to 0)
  double result = 0.0; 
  
  // If s vector, integrate as vector (dot with surface normal)
  if (u-l+1 == nsd) { 
     Array<double> vec(nsd,nNo);
     for (int a = 0; a < nNo; a++) {
       for (int i = l, n = 0; i <= u; i++, n++) {
         vec(n,a) = s(i,a);                 
       }
     }
     result = integ(com_mod, cm_mod, lFa, vec, cfg);
  // If s scalar, integrate as scalar
  } else if (l == u) {
     Vector<double> sclr(nNo);
     for (int a = 0; a < nNo; a++) {
        sclr(a) = s(l,a);
     }
     result = integ(com_mod, cm_mod, lFa, sclr, flag, cfg);
  } else {
    throw std::runtime_error("Unexpected dof in integ");
  }

  return result; 
}


bool is_domain(const ComMod& com_mod, const eqType& eq, const int node, const consts::EquationType phys)
{
  bool result = false;

  // Single domain is assumed and we only need to check that
  //
  if (eq.nDmn == 1) {
    if (eq.dmn[0].phys == phys) {
      result = true;
    } 

  // Domain partition is expected
  //
  } else { 
    if (com_mod.dmnId.size() == 0) {
      throw std::runtime_error("Domain partitioning info is not provided.");
    }
    //  IF (.NOT.ALLOCATED(dmnId)) err = "Domain partitioning info "// 2      "is not provided"

    for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
      if (eq.dmn[iDmn].phys == phys) {
        if (utils::btest(com_mod.dmnId(node),eq.dmn[iDmn].Id)) {
          result = true;
          break;
        }
      }
    }
  }

  return result;
}


/// @brief Computes the JACOBIAN of an element.
//
double jacobian(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x, const Array<double>&Nxi)
{
  double Jac = 0.0;
  Array<double> xXi(nDim,nDim);

  for (int a = 0; a < eNoN; a++) {
    xXi.set_col(0, xXi.col(0) + x.col(a) * Nxi(0,a));
    xXi.set_col(1, xXi.col(1) + x.col(a) * Nxi(1,a));
    // xXi(:,1) = xXi(:,1) + (x(:,a) * Nxi(1,a))
    // xXi(:,2) = xXi(:,2) + (x(:,a) * Nxi(2,a))

    if (nDim == 3) {
      xXi.set_col(2, xXi.col(2) + x.col(a) * Nxi(2,a));
      //xXi(:,3) = xXi(:,3) + (x(:,a) * Nxi(3,a))
    }
  }

  Jac = mat_fun::mat_det(xXi,nDim);
  //Jac = MAT_DET(xXi,nDim)

  return Jac;
}


Vector<int> 
local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Vector<int>& U)
{
  Vector<int> local_vector;

  if (com_mod.ltg.size() == 0) {
    throw std::runtime_error("ltg is not set yet");
  }

  if (cm.mas(cm_mod)) {
    if (U.size() != com_mod.gtnNo) {
      throw std::runtime_error("local_is is only specified for vector with size gtnNo");
    }
  }

  if (cm.seq()) {
    local_vector.resize(com_mod.gtnNo);
    local_vector = U;
    return local_vector; 
  }

  local_vector.resize(com_mod.tnNo); 
  Vector<int> tmpU(com_mod.gtnNo);

  if (cm.mas(cm_mod)) {
    tmpU = U;
  }

  cm.bcast(cm_mod, tmpU);

  for (int a = 0; a < com_mod.tnNo; a++) {
    int Ac = com_mod.ltg[a];
    local_vector[a] = tmpU[Ac];
  }

  return local_vector;
}


Array<double> 
local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Array<double>& U)
{
  //int task_id = cm.idcm();
  //std::string msg_prefix = std::string("[local_rv:") + std::to_string(task_id) + "] ";
  //dmsg;
  //dmsg << "========== local_rv ==========";

  if (com_mod.ltg.size() == 0) {
    throw std::runtime_error("ltg is not set yet");
  }

  Array<double> local_array; 
  int m;

  if (cm.mas(cm_mod)) {
    m = U.nrows();
    if (U.ncols() != com_mod.gtnNo) {
      throw std::runtime_error("local_rv is only specified for vector with size gtnNo");
    }
  }

  if (cm.seq()) {
    local_array.resize(m, com_mod.gtnNo);
    local_array = U;
    return local_array; 
  }

  cm.bcast(cm_mod, &m);

  local_array.resize(m, com_mod.tnNo); 
  Vector<double> tmpU(m * com_mod.gtnNo);
  //ALLOCATE(LOCALRV(m,tnNo), tmpU(m*gtnNo))

  if (cm.mas(cm_mod)) {
    for (int a = 0; a < com_mod.gtnNo; a++) {
      int s = m * a;
      int e = m * (a + 1);
      for (int i = 0; i < U.nrows(); i++) {
        tmpU(i+s) = U(i, a);
      }
      // tmpU(s:e) = U(:,a)
    }
  }

  cm.bcast(cm_mod, tmpU);

  //dmsg << "Compute local_rv ";

  for (int a = 0; a < com_mod.tnNo; a++) {
    int Ac = com_mod.ltg[a];
    int s = m * Ac;
    int e = m * (Ac+1);
    //if (a < 5) dmsg << "localrv Ac: " << Ac;
    // int s  = m * (Ac-1) + 1;
    // int e  = m * Ac;
    // LOCALRV(:,a) = tmpU(s:e)
    // if (a < 5) dmsg << "localrv s and e: " << s << " " << e;
    //if (a < 5) dmsg << "localrv a: " << a;
    for (int i = 0; i < m; i++) {
      //if (a < 5) dmsg << "localrv tmpU(i+s): " << tmpU(i+s);
      local_array(i,a) = tmpU(i+s);
      //if (a < 5) dmsg << "localrv (i,a): " << local_array(i,a);
    }
  }

  return local_array;
}


Array3<double>
local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Array3<double>& U){
  if (com_mod.ltg.size() == 0) {
    throw std::runtime_error("ltg is not set yet");
  }

  Array3<double> local_array;
  int m;
  int n;
  int r;

  if (cm.mas(cm_mod)) {
    m = U.nrows(); // nsd
    r = U.ncols(); // tnNo
    n = U.nslices(); // time
    if (U.ncols() != com_mod.gtnNo) {
      throw std::runtime_error("local_rv is only specified for vector with size gtnNo");
    }
  }

  if (cm.seq()) {
    local_array.resize(m, com_mod.gtnNo, U.nslices());
    local_array = U;
    return local_array;
  }

  cm.bcast(cm_mod, &m);
  cm.bcast(cm_mod, &n);
  cm.bcast(cm_mod, &r);

  local_array.resize(m, com_mod.tnNo, n);
  Vector<double> tmpU(m * com_mod.gtnNo * n);

  if (cm.mas(cm_mod)) {
    for (int a = 0; a < com_mod.gtnNo; a++) {
      int s = m * a;
      for (int i = 0; i < n; i++) {
        int e = i * m * (com_mod.gtnNo);
        for (int j = 0; j < m; j++) {
          tmpU(j+s+e) = U(j, a, i);
        }
      }
    }
  }

  cm.bcast(cm_mod, tmpU);

  for (int a = 0; a < com_mod.tnNo; a++) {
    int Ac = com_mod.ltg[a];
    int s = m * Ac;
    for (int i = 0; i < n; i++) {
      int e = i * m * (com_mod.gtnNo);
      for (int j = 0; j < m; j++) {
        local_array(j, a, i) = tmpU(j+s+e);
      }
    }
  }

  return local_array;
}

Vector<double> 
mkc(const ComMod& com_mod, Vector<double>& U)
{
  if (U.size() != com_mod.lhs.nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Vector<double> result(com_mod.lhs.nNo);

  if (com_mod.cm.seq()) {
    result = U;
  } else {  
    for (int a = 0; a < com_mod.lhs.nNo; a++) {
      result(com_mod.lhs.map(a)) = U(a);
    }
  }

  return result;
}

Array<double> 
mkc(const ComMod& com_mod, Array<double>& U)
{
  int m = U.nrows();
  if (U.ncols() != com_mod.lhs.nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Array<double> result(m, com_mod.lhs.nNo);

  if (com_mod.cm.seq()) {
    result = U;
  } else {
    for (int a = 0; a < com_mod.lhs.nNo; a++) {
      for (int i = 0; i < m; i++) {
        result(i, com_mod.lhs.map(a)) = U(i, a);
      }
    }
  }

  return result;
}


void mkci(const ComMod& com_mod, Vector<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  if (U.size() != com_mod.lhs.nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Vector<double> tmp(com_mod.lhs.nNo);
  tmp = U;

  for (int a = 0; a < com_mod.lhs.nNo; a++) {
    U(a) = tmp(com_mod.lhs.map(a));
  }
}

void mkci(const ComMod& com_mod, Array<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  auto& lhs = com_mod.lhs;
  int nNo = lhs.nNo;
  int m = U.nrows();

  if (U.ncols() != nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Array<double> tmp(m, nNo);
  tmp = U;

  for (int a = 0; a < nNo; a++) {
    for (int i = 0; i < m; i++) {
      U(i,a) = tmp(i,lhs.map(a));
    } 
  }
}

/// @brief Set domain ID to a given number for the entire or a range of elements in a mesh.
//
void set_dmn_id(mshType& mesh, const int iDmn, const int ifirst, const int ilast)
{
  int first = 0;
  if (ifirst != consts::int_inf) {
    first = ifirst;
  }
  if ((first < 0) || (first > mesh.gnEl-1)) {
    throw std::runtime_error("Setting domain ID with a start element range " + std::to_string(first) + ".");
  }

  int last = mesh.gnEl - 1;
  if (ilast != consts::int_inf) {
    last = ilast;
  }
  if ((last < first) || (last > mesh.gnEl-1)) {
    throw std::runtime_error("Setting domain ID with a end element range " + std::to_string(last) + ".");
  }

  if (mesh.eId.size() == 0) {
    mesh.eId = Vector<int>(mesh.gnEl);
  }

  // Set the iDimn'th bit for each element ID.
  for (int e = first; e <= last; e++) {
    mesh.eId[e] |= 1UL << iDmn;
  }
}

/// @brief Computes the Skewness of an element.
//
double skewness(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x)
{
  Vector<double> coeff;

  if (nDim == 2) {
    coeff.set_values({1.0, -1.0, 1.0, -1.0});
  } else {
    coeff.set_values({1.0, 1.0, -1.0, 1.0, 1.0});
  }

  Array<double> Dmat(eNoN,nDim+2);
  Dmat = 1.0;
  Array<double> Dsub(eNoN,nDim+1);

  for (int a = 0; a < eNoN; a++) {
    auto col = x.col(a);
    Dmat(a,0) = col * col;

    for (int i = 0; i < nDim; i++) {
      Dmat(a, i+1) = x(i,a);
    }
  }

  Vector<double> detD(nDim+2);

  for (int j = 0; j < nDim+2; j++) {
    int cnt = 0;
    for (int i = 0; i < nDim+2; i++) {
      if (i == j) {
        continue;
      }
      Dsub.set_col(cnt, Dmat.col(i));
      cnt = cnt + 1;
    }
    detD(j) = coeff(j) * mat_fun::mat_det(Dsub,nDim+1);
  }

  double circumRad = 0.0; 
  for (int i = 0; i < nDim; i++) {
    circumRad += detD(i+1)*detD(i+1);
  }
  circumRad = sqrt(circumRad - 4.0*detD(0)*detD(nDim+1)) / (2.0*fabs(detD(0)));

  double integ_eq, integ_el;

  if (nDim == 2) {
    integ_eq = 0.25 * sqrt(27.0) * circumRad*circumRad;
    integ_el = 0.5 * fabs(detD(0));
  } else if (nDim == 3) {
    integ_eq = 80.0 * pow(circumRad,3.0) / sqrt(243.0);
    integ_el = fabs(detD(0)) / 6.0;
  }

  return fabs(integ_eq - integ_el) / integ_eq;
}


/// @brief Spliting "m" jobs between "n" workers. "b" contains amount of jobs
/// and "A" will store the distribution of jobs
///
/// A(nMsh, num_proc)
///
/// b(nMsh)
///
/// n = nMsh
///
/// m = num_proc
///
/// Replicates 'RECURSIVE SUBROUTINE SPLITJOBS(m, n, A, b)' in ALLFUN.f.
//
void split_jobs(int tid, int m, int n, Array<double>& A, Vector<double>& b)
{
  #define n_debug_split_jobs
  #ifdef debug_split_jobs
  DebugMsg dmsg(__func__, tid);
  dmsg.banner();
  dmsg << "m: " << m;
  dmsg << "n: " << n;
  dmsg << "b: " << b;
  #endif

  if ((m <= 0) || (n <= 0)) {
    return;
  }

  // Single worker, but multiple jobs.
  //
  if (n == 1) {
    #ifdef debug_split_jobs
    dmsg << "n == 1  " << " ";
    #endif
    int j = 0;
    for (int i = 0; i < m; i++) {
      A(i,j) = b[i];
      #ifdef debug_split_jobs
      dmsg << "n=1 A(i,0): " << A(i,0);
      #endif
    }
    return;
  }

  // Multiple workers, but a single job.
  //
  if (m == 1) {
    #ifdef debug_split_jobs
    dmsg << "m == 1  " << " ";
    #endif
    int i = 0;
    for (int j = 0; j < n; j++) {
      A(i,j) = b[i] / static_cast<double>(n); 
      #ifdef debug_split_jobs
      dmsg << "m=1 A(i,j): " << A(i,j);
      #endif
    }
    return;
  }

  // Multiple workers and multiple jobs
  // This is the initial guess for nl, nr
  int nl  = n / 2;
  int nr  = n - nl;

  // This is the total amount of work
  double sb  = b.sum();

  // The work that suppose to be done by "l"
  double sbl = sb * static_cast<double>(nl) / static_cast<double>(n); 
  #ifdef debug_split_jobs
  dmsg << "nl: " << nl;
  dmsg << "nr: " << nr;
  dmsg << "sb: " << sb;
  dmsg << "sbl: " << sbl;
  #endif

  double sl = 0.0;
  int ival = 0;
  for (int i = 0; i < m; i++) {
    if (sl + b[i] > sbl) {
      ival = i;
      break;
    }
    sl = sl + b[i];
  }

  double optsl = fabs(sl - sbl);
  int ml = ival;
  #ifdef debug_split_jobs
  dmsg;
  dmsg << "optsl: " << optsl;
  dmsg << "ml: " << ml;
  dmsg;
  dmsg << "Set j ... " << " ";
  #endif

  int j = -1;

  for (int i = ml; i < m; i++) {
    #ifdef debug_split_jobs
    dmsg << "---- i " << i;
    #endif
    if (fabs(sl + b[i] - sbl) < optsl) {
      j = i;
      optsl = fabs(sl + b[i] - sbl);
      #ifdef debug_split_jobs
      dmsg << "j: " << j;
      dmsg << "optsl: " << optsl;
      #endif
    }
  }

  if (j !=  -1) {
    ml = ml + 1;
  }

  int mr = m - ml;
  Vector<double> bl(ml); 
  Vector<double> br(mr);

  #ifdef debug_split_jobs
  dmsg;
  dmsg << "j: " << j;
  dmsg << "mr: " << mr;
  dmsg << "ml: " << ml;
  #endif

  if (j != -1) {
    for (int i = 0; i < ml; i++) {
      bl[i] = b[i];
    }
    bl[ml-1] = b[j];
    int k = 0;
    for (int i = ml-1; i < m; i++) {
      if (i == j) {
        continue; 
      }
      br[k] = b[i];
      k = k + 1;
    }
  } else { 
    for (int i = 0; i < ml; i++) {
      bl[i] = b[i];
    }
    for (int i = ml, j = 0; i < m; i++, j++) {
      br[j] = b[i];
    }
  }

  nl = round(static_cast<double>(n) * bl.sum() / sb);
  if (nl == 0) {
    nl = 1;
  }
  if (nl == n) { 
    nl = n - 1;
  }
  nr = n - nl;

  #ifdef debug_split_jobs
  dmsg << "  " << " ";
  dmsg << "nl: " << nl;
  dmsg << "nr: " << nr;
  dmsg << "  " << " ";
  dmsg << "Allocate Al: ml: " << ml;
  dmsg << "             nl: " << nl;
  dmsg << "             bl: " << bl;
  dmsg << "  " << " ";
  dmsg << "Allocate Ar: mr: " << mr;
  dmsg << "             nr: " << nr;
  dmsg << "             br: " << br;
  //dmsg << "Allocate A: " << A.nrows() << " x " << A.ncols();
  #endif

  Array<double> Al(ml,nl); 
  split_jobs(tid, ml, nl, Al, bl);

  Array<double> Ar(mr,nr);
  split_jobs(tid, mr, nr, Ar, br);


  A = 0.0;

  if (j != -1) {
    #ifdef debug_split_jobs
    dmsg << "" << " ";
    dmsg << "set A from Al ... ml: " << ml;
    #endif
    for (int i = 0; i < ml-1; i++) {
      auto Al_row = Al.row(i);
      A.set_row(i, Al_row);
      #ifdef debug_split_jobs
      for (int ii = 0; ii < n; ii++) {
        dmsg << " Al set A(i,ii): " << A(i,ii);
      }
      #endif
    }

    A.set_row(j, Al.row(ml-1));
    #ifdef debug_split_jobs
    dmsg << "----------" << " ";
    for (int ii = 0; ii < n; ii++) {
      dmsg << " set from Al:  A(j,ii): " << A(j,ii);
    }
    #endif

    int k = 0;
    for (int i = ml-1; i < m; i++) {
      if (i == j) {
        continue; 
      }
      A.set_row(i, nl, Ar.row(k));
      k = k + 1;
    }
  } else { 
    // [TODO:DaveP] another bug fix
    #ifdef debug_split_jobs
    dmsg << "----------" << " ";
    dmsg << "Set A from Al and Ar " << " ";
    dmsg << "m: " << m;
    dmsg << "ml: " << ml;
    dmsg << "mr: " << mr;
    #endif
    for (int i = 0; i < ml; i++) {
      for (int j = 0; j < nl; j++) {
        A(i, j) = Al(i, j);
       }
     }

    for (int i = 0; i < mr; i++) {
      for (int j = 0; j < nr; j++) {
        A(i+ml, j+nl) = Ar(i, j);
      }
    }

    // A(1:ml,1:nl) = Al
    // A(ml+1:m,nl+1:n) = Ar
  }

  #ifdef debug_split_jobs
  dmsg << "Returned A: " << A;
  #endif
}


};
