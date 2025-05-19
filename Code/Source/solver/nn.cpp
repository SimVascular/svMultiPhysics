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

// The functions defined here replicate the Fortran functions defined in NN.f.
//
// The functions are used to 
//
//   1) Set element properties: element type, number of Gauss integration points, ... 
//
//   2) Allocate element arrays: Gauss weights, shape functions, ...
//

#include "nn.h"

#include "Array.h"
#include "Vector.h"

#include "consts.h"
#include "mat_fun.h"
#include "utils.h"

#include "lapack_defs.h"

#include <functional>
#include <iostream> 
#include <math.h> 

namespace nn {

#define dgb_nn

using namespace consts;

// Define maps used to set element properties.
#include "nn_elem_props.h"

// Define maps used to set element Gauss integration data. 
#include "nn_elem_gip.h"

// Define maps used to set element shape function data. 
#include "nn_elem_gnn.h"

// Define maps used to get element shape function 2nd derivative data. 
#include "nn_elem_gnnxx.h"

// Define a map type used to set the bounds of element shape functions.
#include "nn_elem_nn_bnds.h"

void get_gip(const int insd, consts::ElementType eType, const int nG, Vector<double>& w, Array<double>& xi) 
{
  try {
    get_element_gauss_int_data[eType](insd, nG, w, xi);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("No support for element etype " + std::to_string(static_cast<int>(eType)) + 
        " in 'get_element_gauss_int_data'.");
  }
}

/// @brief Define Gauss integration points in local (ref) coordinates.
///
/// \todo [NOTE] There should just have a single map for mesh and face types.
//
void get_gip(mshType& mesh)
{
  try {
    set_element_gauss_int_data[mesh.eType](mesh);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("No support for mesh etype " + std::to_string(static_cast<int>(mesh.eType)) + " in 'set_element_gauss_int_data'.");
  }
}

void get_gip(Simulation* simulation, faceType& face)
{
  try {
    set_face_gauss_int_data[face.eType](face);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("No support for face type " + std::to_string(static_cast<int>(face.eType)) + " in 'set_face_gauss_int_data'.");
  }
}

/// Computes shape functions and derivatives at given natural coords.
//
void get_gnn(const int insd, consts::ElementType eType, const int eNoN, const int g, Array<double>& xi, 
    Array<double>& N, Array3<double>& Nx)
{
  try {
    get_element_shape_data[eType](insd, eNoN, g, xi, N, Nx);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[get_gnn] No support for element type " + std::to_string(static_cast<int>(eType)) + " in 'get_element_shape_data'.");
  }
}

/// @brief A big fat hack because the Fortran GETNN() operates on primitive types but
/// the C++ version does not, uses Array and Vector objects.
//
void get_gnn(const int nsd, consts::ElementType eType, const int eNoN, Vector<double>& xi, 
    Vector<double>& N, Array<double>& Nx)
{
  int size = xi.size();
  Array<double> xi_a(size,1);
  xi_a.set_col(0, xi);
  Array<double> N_a(eNoN,1);
  Array3<double> Nx_a(size,eNoN, 1);

  nn::get_gnn(nsd, eType, eNoN, 0, xi_a, N_a, Nx_a);

  xi = xi_a.col(0);
  N = N_a.col(0);
  Nx = Nx_a.slice(0);
}

void get_gnn(int gaus_pt, mshType& mesh)
{
  try {
    set_element_shape_data[mesh.eType](gaus_pt, mesh);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[get_gnn] No support for element type " + std::to_string(static_cast<int>(mesh.eType)) + " in 'set_element_shape_data'.");
  }
}

void get_gnn(Simulation* simulation, int gaus_pt, faceType& face)
{
  try {
    set_face_shape_data[face.eType](gaus_pt, face);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("No support for face type " + std::to_string(static_cast<int>(face.eType)) + " in 'set_face_shape_data'.");
  }
}

/// @brief Returns second order derivatives at given natural coords
///
/// Replicates 'SUBROUTINE GETGNNxx(insd, ind2, eType, eNoN, xi, Nxx)'.
//
void get_gn_nxx(const int insd, const int ind2, consts::ElementType eType, const int eNoN, const int gaus_pt, 
    const Array<double>& xi, Array3<double>& Nxx)
{
  using namespace consts;

  // Element types that don't have 2nd derivatives computed for them.
  static std::set<ElementType> no_derivs{ElementType::NRB, ElementType::QUD4, ElementType::HEX8, 
                                         ElementType::HEX20, ElementType::HEX27};

  if (no_derivs.count(eType) != 0) {
    return;
  }

  try {
    get_element_2nd_derivs[eType](insd, ind2, eNoN, gaus_pt, xi, Nxx);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("[get_gn_nxx] No support for element type " + std::to_string(static_cast<int>(eType)) + " in 'get_element_2nd_derivs'.");
  }
}

/// @brief Sets bounds on Gauss integration points in parametric space and
/// bounds on shape functions.
///
/// Reproduces the Fortran 'GETNNBNDS' subroutine.
//
void get_nn_bnds(const int nsd, consts::ElementType eType, const int eNoN, Array<double>& xib, Array<double>& Nb)
{
  using namespace consts;

  for (int i = 0; i < nsd; i++) {
    xib(0,i) = -1.0; 
    xib(1,i) = 1.0; 
  }

  for (int i = 0; i < eNoN; i++) {
    Nb(0,i) = 0.0; 
    Nb(1,i) = 1.0; 
  }

  switch (eType) {

    case ElementType::HEX20:
      for (int i = 0; i < 20; i++) {
        Nb(0,i) = -0.125;
      }
    break;

    case ElementType::HEX27:
      for (int i = 0; i < 20; i++) {
        Nb(0,i) = -0.125;
      }
      Nb(0,26) = 0.0;
    break;

    case ElementType::LIN2:
      Nb(0,0) = -0.125;
      Nb(0,1) = -0.125;
      Nb(0,2) = 0.0;
    break;

    case ElementType::QUD8:
      for (int i = 0; i < 8; i++) {
        Nb(0,i) = -0.125;
      }
    break;

    case ElementType::QUD9:
      for (int i = 0; i < 8; i++) {
        Nb(0,i) = -0.125;
      }
      Nb(0,8) = 0.0;
    break;

    case ElementType::TET4:
      for (int i = 0; i < nsd; i++) {
        xib(0,i) = 0.0; 
      }
    break;

    case ElementType::TET10:
      for (int i = 0; i < nsd; i++) {
        xib(0,i) = 0.0; 
      }

      for (int i = 0; i < 4; i++) {
        Nb(0,i) = -0.125;
      }
      for (int i = 4; i < 10; i++) {
        Nb(1,i) = 4.0;
      }
    break;

    case ElementType::TRI3:
      for (int i = 0; i < nsd; i++) {
        xib(0,i) = 0.0; 
      }
    break;

    case ElementType::TRI6:
      for (int i = 0; i < nsd; i++) {
        xib(0,i) = 0.0; 
      }

      for (int i = 0; i < 3; i++) {
        Nb(0,i) = -0.125;
      }

      for (int i = 3; i < 6; i++) {
        Nb(1,i) = 4.0;
      }
    break;

    case ElementType::WDG:
      xib(0,0) = 0.0;
      xib(0,1) = 0.0;
    break;

    default:
    break;
  }

  // Add a small tolerance around the bounds
  double tol = 1.0E-4;
  for (int i = 0; i < nsd; i++) {
    xib(0,i) -= tol; 
    xib(1,i) += tol; 
  }

  for (int i = 0; i < eNoN; i++) {
    Nb(0,i) -= tol; 
    Nb(1,i) += tol; 
  }
}

/// @brief Sets bounds on Gauss integration points in parametric space and
/// bounds on shape functions.
///
/// Modifies:
///   mesh % xib(2,nsd) - Bounds on Gauss integration points in parametric space
///   mesh % Nb(2,mesh % eNoN) - Bounds on shape functions
/// 
///
/// Replicates Fortran SUBROUTINE GETNNBNDS.
//
void get_nn_bnds(const ComMod& com_mod, mshType& mesh)
{
  int nsd = com_mod.nsd;
  auto eType = mesh.eType;
  int eNoN = mesh.eNoN;

  get_nn_bnds(nsd, eType, eNoN, mesh.xib, mesh.Nb);
}

void get_nnx(const int nsd, const consts::ElementType eType, const int eNoN, const Array<double>& xl, 
    const Array<double>& xib, const Array<double>& Nb, const Vector<double>& xp, Vector<double>& xi, 
    Vector<double>& N, Array<double>& Nx)
{
  #define n_debug_get_nnx 
  #ifdef debug_get_nnx 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "eType: " << eType;
  dmsg << "eNoN: " << eNoN;
  dmsg << "xl: " << xl;
  dmsg << "xib: " << xib;
  dmsg << "Nb: " << Nb;
  dmsg << "xp: " << xp;
  #endif

  bool l1;
  get_xi(nsd, eType, eNoN, xl, xp, xi, l1); 

  // Check if parameteric coordinate is within bounds
  int j = 0;

  for (int i = 0; i < nsd; i++) { 
    if (xi(i) >= xib(0,i) && xi(i) <= xib(1,i)) {
      j = j + 1;
    }
  }

  bool l2 = (j == nsd);

  get_gnn(nsd, eType, eNoN, xi, N, Nx);

  // Check if shape functions are within bounds and sum to unity
  j = 0;
  double rt = 0.0;

  for (int i = 0; i < eNoN; i++) {
    rt = rt + N(i);
    if (N(i) > Nb(0,i) && N(i) < Nb(1,i)) {
      j = j + 1;
    }
  }

  bool l3 = (j == eNoN);
  bool l4 = (rt >= 0.9999) && (rt <= 1.0001);

  l1 = (l1 && l2 && l3 && l4);

  if (!l1) {
    throw std::runtime_error("Error in computing shape functions");
  }
}

/// @brief Inverse maps {xp} to {$\xi$} in an element with coordinates {xl} using Newton's method
//
void get_xi(const int nsd, consts::ElementType eType, const int eNoN, const Array<double>& xl, const Vector<double>& xp, 
    Vector<double>& xi, bool& flag)
{
  static const int MAXITR = 5;
  static const double RTOL = 1.E-6;
  static const double ATOL = 1.E-12;

  Vector<double> xK(nsd);  
  Vector<double> rK(nsd);
  Vector<double> N(eNoN);
  Array<double> Nxi(nsd,eNoN);

  int itr = 0;
  auto xiK = xi;
  double eps = std::numeric_limits<double>::epsilon();
  bool l1, l2, l3;

  while (true) { 
     itr = itr + 1;
     nn::get_gnn(nsd, eType, eNoN, xiK, N, Nxi);
     xK = 0.0;

     for (int i = 0; i < nsd; i++) {
        for (int a = 0; a < eNoN; a++) {
           xK(i) = xK(i) + N(a)*xl(i,a);
        }
        rK(i) = xK(i) - xp(i);
     }

     double rmsA = 0.0;
     double rmsR = 0.0;

     for (int i = 0; i < nsd; i++) {
        rmsA = rmsA + pow(rK(i), 2.0);
        rmsR = rmsR + pow(rK(i) / (xK(i)+eps), 2.0);
     }
     rmsA = sqrt(rmsA/static_cast<double>(nsd));
     rmsR = sqrt(rmsR/static_cast<double>(nsd));

     l1 = itr > MAXITR;
     l2 = rmsA <= ATOL;
     l3 = rmsR <= RTOL;
     if (l1 || l2 || l3) {
       break;
     }

     Array<double> Am(nsd, nsd);

     for (int i = 0; i < nsd; i++) {
       for (int j = 0; j < nsd; j++) {
         for (int a = 0; a < eNoN; a++) {
           Am(i,j) = Am(i,j) + xl(i,a)*Nxi(j,a);
         }
       }
     }

     Am  = mat_fun::mat_inv(Am, nsd);
     rK  = mat_fun::mat_mul(Am, rK);
     xiK = xiK - rK;
  }

  // Newton's method converges
  if (l2 || l3) {
     flag = true; 

  // Newton's method failed to converge
  } else { 
    flag = false; 
  }

  xi = xiK;
}

/**
 * @brief Get shape function gradient w.r.t reference configuration coordinates. 
 * dN/dX = dN/dxi * dxi/dX, where xi is the parametric coordinates of the parent element
 * 
 * @param eNoN Number of nodes in the element
 * @param nsd Number of spatial dimensions
 * @param insd 
 * @param Nxi Shape function gradient w.r.t. parametric coordinates of the parent element
 * @param x Element node positions in reference configuration
 * @param Nx Shape function gradient w.r.t reference configuration coordinates
 * @param Jac Jacobian of element in reference configuration
 * @param ks Matrix where each component is (dxi/dX)^2
 */
void gnn(const int eNoN, const int nsd, const int insd, Array<double>& Nxi, Array<double>& x, Array<double>& Nx, 
    double& Jac, Array<double>& ks)
{
  Array<double> xXi(nsd,insd);   
  Array<double> xiX(insd,nsd);

  Jac = 0.0;
  Nx  = 0.0;
  ks  = 0.0;
  double eps = std::numeric_limits<double>::epsilon();

  if (insd == 1) {
    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < nsd; i++) {
        xXi(i,0) = xXi(i,0) + x(i,a)*Nxi(0,a);
      }
    }

    Jac = sqrt(utils::norm(xXi)) + 1.E+3*eps;
    for (int a = 0; a < eNoN; a++) {
      Nx(0,a) = Nxi(0,a) / Jac;
    }

  } else if (insd == 2) {
    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < nsd; i++) {
        xXi(i,0) = xXi(i,0) + x(i,a)*Nxi(0,a);
        xXi(i,1) = xXi(i,1) + x(i,a)*Nxi(1,a);
      }
    }

    Jac = xXi(0,0)*xXi(1,1) - xXi(0,1)*xXi(1,0);

    xiX(0,0) =  xXi(1,1) / Jac;
    xiX(0,1) = -xXi(0,1) / Jac;
    xiX(1,0) = -xXi(1,0) / Jac;
    xiX(1,1) =  xXi(0,0) / Jac;

    ks(0,0) = xiX(0,0)*xiX(0,0) + xiX(1,0)*xiX(1,0);
    ks(0,1) = xiX(0,0)*xiX(0,1) + xiX(1,0)*xiX(1,1);
    ks(1,1) = xiX(0,1)*xiX(0,1) + xiX(1,1)*xiX(1,1);
    ks(1,0) = ks(0,1);

    for (int a = 0; a < eNoN; a++) {
      Nx(0,a) = Nx(0,a)+ Nxi(0,a)*xiX(0,0) + Nxi(1,a)*xiX(1,0);
      Nx(1,a) = Nx(1,a)+ Nxi(0,a)*xiX(0,1) + Nxi(1,a)*xiX(1,1);
    }

  } else if (insd == 3) {
    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < xXi.nrows(); i++) {
        xXi(i,0) = xXi(i,0) + x(i,a)*Nxi(0,a);
        xXi(i,1) = xXi(i,1) + x(i,a)*Nxi(1,a);
        xXi(i,2) = xXi(i,2) + x(i,a)*Nxi(2,a);
      }
    }

    Jac = xXi(0,0)*xXi(1,1)*xXi(2,2) + xXi(0,1)*xXi(1,2)*xXi(2,0) + xXi(0,2)*xXi(1,0)*xXi(2,1) - 
          xXi(0,0)*xXi(1,2)*xXi(2,1) - xXi(0,1)*xXi(1,0)*xXi(2,2) - xXi(0,2)*xXi(1,1)*xXi(2,0);

    xiX(0,0) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac;
    xiX(0,1) = (xXi(2,1)*xXi(0,2) - xXi(2,2)*xXi(0,1))/Jac;
    xiX(0,2) = (xXi(0,1)*xXi(1,2) - xXi(0,2)*xXi(1,1))/Jac;
    xiX(1,0) = (xXi(1,2)*xXi(2,0) - xXi(1,0)*xXi(2,2))/Jac;
    xiX(1,1) = (xXi(2,2)*xXi(0,0) - xXi(2,0)*xXi(0,2))/Jac;
    xiX(1,2) = (xXi(0,2)*xXi(1,0) - xXi(0,0)*xXi(1,2))/Jac;
    xiX(2,0) = (xXi(1,0)*xXi(2,1) - xXi(1,1)*xXi(2,0))/Jac;
    xiX(2,1) = (xXi(2,0)*xXi(0,1) - xXi(2,1)*xXi(0,0))/Jac;
    xiX(2,2) = (xXi(0,0)*xXi(1,1) - xXi(0,1)*xXi(1,0))/Jac;

    ks(0,0) = xiX(0,0)*xiX(0,0)+xiX(1,0)*xiX(1,0)+xiX(2,0)*xiX(2,0);
    ks(0,1) = xiX(0,1)*xiX(0,0)+xiX(1,1)*xiX(1,0)+xiX(2,1)*xiX(2,0);
    ks(0,2) = xiX(0,2)*xiX(0,0)+xiX(1,2)*xiX(1,0)+xiX(2,2)*xiX(2,0);
    ks(1,1) = xiX(0,1)*xiX(0,1)+xiX(1,1)*xiX(1,1)+xiX(2,1)*xiX(2,1);
    ks(1,2) = xiX(0,1)*xiX(0,2)+xiX(1,1)*xiX(1,2)+xiX(2,1)*xiX(2,2);
    ks(2,2) = xiX(0,2)*xiX(0,2)+xiX(1,2)*xiX(1,2)+xiX(2,2)*xiX(2,2);
    ks(1,0) = ks(0,1);
    ks(2,0) = ks(0,2);
    ks(2,1) = ks(1,2);

    for (int a = 0; a < eNoN; a++) {
      Nx(0,a) = Nx(0,a) + Nxi(0,a)*xiX(0,0) + Nxi(1,a)*xiX(1,0) + Nxi(2,a)*xiX(2,0);
      Nx(1,a) = Nx(1,a) + Nxi(0,a)*xiX(0,1) + Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1);
      Nx(2,a) = Nx(2,a) + Nxi(0,a)*xiX(0,2) + Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2);
    }
  }
}

/// @brief This routine returns a surface normal vector at element "e" and Gauss point
/// 'g' of face 'lFa' that is the normal weighted by Jac, i.e.
/// Jac = SQRT(NORM(n)), the Jacobian of the mapping from parent surface element to
/// reference/old/new configuration.
///
/// cfg denotes which configuration (reference/timestep 0, old/timestep n, or new/timestep n+1). Default reference
///
/// Reproduce Fortran 'GNNB'.
//
void gnnb(const ComMod& com_mod, const faceType& lFa, const int e, const int g, const int nsd, const int insd, 
    const int eNoNb, const Array<double>& Nx, Vector<double>& n, MechanicalConfigurationType cfg)
{
  auto& cm = com_mod.cm;

  #define n_debug_gnnb 
  #ifdef debug_gnnb 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "e: " << e+1;
  dmsg << "g: " << g+1;
  dmsg << "nsd: " << nsd;
  dmsg << "insd: " << insd;
  dmsg << "eNoNb: " << eNoNb;
  dmsg << "cfg: " << cfg;
  #endif

  int iM = lFa.iM;
  int Ec = lFa.gE(e);
  auto& msh = com_mod.msh[iM];
  int eNoN = msh.eNoN;

  #ifdef debug_gnnb 
  dmsg << "iM: " << iM+1;
  dmsg << "Ec: " << Ec+1;
  dmsg << "eNoN: " << eNoN;
  dmsg << "msh.IEN.nrows: " << msh.IEN.nrows();
  dmsg << "msh.IEN.ncols: " << msh.IEN.ncols();
  #endif

  Array<double> lX(nsd,eNoN); 
  Vector<int> ptr(eNoN); 
  std::vector<bool> setIt(eNoN);

  // Creating a ptr list that contains pointer to the nodes of elements
  // that are at the face at the beginning of the list and the rest at
  // the end
  //
  std::fill(setIt.begin(), setIt.end(), true);

  for (int a = 0; a < eNoNb; a++) {
    int Ac = lFa.IEN(a,e);
    int b = 0;
    bool found_node = false;

    for (int ib = 0; ib < eNoN; ib++) {
      b = ib;
      if (setIt[ib]) {
        int Bc = msh.IEN(ib,Ec);
        if (Bc == Ac) {
          found_node = true;
          break;
        }
      }
    }

    if (!found_node) {
      throw std::runtime_error("[svMultiPhysics::gnnb] ERROR: The '" + lFa.name + "' face node " + std::to_string(Ac) + 
          " could not be matched to a node in the '" + msh.name + "' volume mesh.");
    }

    ptr(a) = b;
    setIt[b] = false;
  }

  int a = eNoNb;

  for (int b = 0; b < eNoN; b++) {
    if (setIt[b]) {
      ptr(a) = b;
      a = a + 1;
    }
  }

  // Correcting the geometry if mesh is moving or if we want the 
  // area-weighted normal in a different configuration
  //
  for (int a = 0; a < eNoN; a++) {
    int Ac = msh.IEN(a,Ec);
    for (int i = 0; i < lX.nrows(); i++) {
      // Get position vector
      lX(i,a) = com_mod.x(i,Ac);
    }
    if (com_mod.mvMsh) {
      for (int i = 0; i < lX.nrows(); i++) {
        // Add mesh displacement
        lX(i,a) = lX(i,a) + com_mod.Do(i+nsd+1,Ac);
      }
    }
    else {
      switch (cfg) {
        case MechanicalConfigurationType::reference:
          // Do nothing
          break;
        case MechanicalConfigurationType::old_timestep:
          for (int i = 0; i < lX.nrows(); i++) {
            // Add displacement at timestep n
            lX(i,a) = lX(i,a) + com_mod.Do(i,Ac);
          }
          break;
        case MechanicalConfigurationType::new_timestep:
          for (int i = 0; i < lX.nrows(); i++) {
            // Add displacement at timestep n+1
            lX(i,a) = lX(i,a) + com_mod.Dn(i,Ac);
          }
          break;
        default:
          throw std::runtime_error("gnnb: invalid MechanicalConfigurationType provided");
      }
    }
  }

  // Calculating surface deflation
  if (msh.lShl) {
    // Since the face has only one parametric coordinate (edge), find
    // its normal from cross product of mesh normal and interior edge

    // Update shape functions if NURBS
    if (msh.eType == ElementType::NRB) {
     //  CALL NRBNNX(msh(iM), Ec)
    }

    // Compute adjoining mesh element normal
    //
    Array<double> xXi(nsd,nsd-1);

    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < insd; i++) {
        for (int j = 0; j < nsd; j++) {
          xXi(j,i) = xXi(j,i) + lX(j,a)*msh.Nx(i,a,g);
        }
      }
    }

    auto v = utils::cross(xXi);
    for (int i = 0; i < nsd; i++) {
      v(i) = v(i) / sqrt(utils::norm(v));
    }

    // Face element surface deflation
    xXi.resize(nsd,1);
    for (int a = 0; a < eNoNb; a++) {
      int b = ptr(a);
      for (int i = 0; i < nsd; i++) {
        xXi(i,0) = xXi(i,0) + lFa.Nx(0,a,g)*lX(i,b);
      }
    }

    // Face normal
    n(0) = v(1)*xXi(2,0) - v(2)*xXi(1,0);
    n(1) = v(2)*xXi(0,0) - v(0)*xXi(2,0);
    n(2) = v(0)*xXi(1,0) - v(1)*xXi(0,0);

    // I choose Gauss point of the mesh element for calculating
    // interior edge
    v = 0.0;
    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < nsd; i++) {
        v(i) = v(i) + lX(i,a)*msh.N(a,g);
      }
    }

    int a = ptr(0);
    for (int i = 0; i < nsd; i++) {
      v(i) = lX(i,a) - v(i);
    }

    if (utils::norm(n,v) < 0.0) {
      n = -n;
    }

    return;

  } else {

    Array<double> xXi(nsd,insd);

    for (int a = 0; a < eNoNb; a++) {
      int b = ptr(a);
      for (int i = 0; i < insd; i++) {
        for (int j = 0; j < nsd; j++) {
          xXi(j,i) = xXi(j,i) + Nx(i,a)*lX(j,b);
        }
      }
    }

    n = utils::cross(xXi);
  }

  // Changing the sign if neccessary. 'a' locates on the face and 'b'
  // in the interior of the element. v points outward along ba
  //
  a = ptr(0);
  int b = ptr(lFa.eNoN);
  Vector<double> v(nsd);

  for (int i = 0; i < nsd; i++) {
    v(i) = lX(i,a) - lX(i,b);
  }

  if (utils::norm(n,v) < 0.0) {
    n = -n;
  }
}

/// @brief Compute shell kinematics: normal vector, covariant & contravariant basis vectors
///
/// Replicates 'SUBROUTINE GNNS(eNoN, Nxi, xl, nV, gCov, gCnv)' defined in NN.f.
//
void gnns(const int nsd, const int eNoN, const Array<double>& Nxi, Array<double>& xl, Vector<double>& nV, 
    Array<double>& gCov, Array<double>& gCnv) 
{
  int insd = nsd - 1;

  Array<double> xXi(nsd,insd); 
  Array<double> Gmat(insd,insd);

  // Calculating surface deflation
  //
  for (int a = 0; a < eNoN; a++) {
    for (int i = 0; i < insd; i++) {
      for (int j = 0; j < nsd; j++) {
        xXi(j,i) = xXi(j,i) + xl(j,a)*Nxi(i,a);
      }
    }
  }

  nV = utils::cross(xXi);

  // Covariant basis
  //
  gCov = xXi;

  // Metric tensor g_i . g_j
  for (int i = 0; i < insd; i++) {
    for (int j = 0; j < insd; j++) {
      for (int a = 0; a < eNoN; a++) {
        Gmat(i,j) = Gmat(i,j) + gCov(a,i)*gCov(a,j);
      }
    }
  }

  // Contravariant basis
  //
  Gmat = mat_fun::mat_inv(Gmat, insd);

  for (int i = 0; i < insd; i++) {
    for (int j = 0; j < insd; j++) {
      for (int k = 0; k < nsd; k++) {
        gCnv(k,i) += Gmat(i,j)*gCov(k,j);
      }
    }
  }
}

/// @brief Compute second order derivative on parent element
//
void gn_nxx(const int l, const int eNoN, const int nsd, const int insd, Array<double>& Nxi, Array<double>& Nxi2, Array<double>& lx,
    Array<double>& Nx, Array<double>& Nxx)
{
  Array<double> xXi(nsd,insd); 
  Array<double> xXi2(nsd,l); 
  Array<double> K(l,l); 
  Array<double> B(l,eNoN);

  double t = 2.0;

  if (insd == 2) {
    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < nsd; i++) {
        xXi(i,0) = xXi(i,0) + lx(i,a)*Nxi(0,a);
        xXi(i,1) = xXi(i,1) + lx(i,a)*Nxi(1,a);
        xXi2(i,0) = xXi2(i,0) + lx(i,a)*Nxi2(0,a);
        xXi2(i,1) = xXi2(i,1) + lx(i,a)*Nxi2(1,a);
        xXi2(i,2) = xXi2(i,2) + lx(i,a)*Nxi2(2,a);
      }
    }

    K.set_row(0, {xXi(0,0)*xXi(0,0), xXi(1,0)*xXi(1,0), t*xXi(0,0)*xXi(1,0)});
    K.set_row(1, {xXi(0,1)*xXi(0,1), xXi(1,1)*xXi(1,1), t*xXi(0,1)*xXi(1,1)});
    K.set_row(2, {xXi(0,0)*xXi(0,1), xXi(1,0)*xXi(1,1), xXi(0,0)*xXi(1,1) + xXi(0,1)*xXi(1,0)});

    for (int a = 0; a < eNoN; a++) {
      B(0,a) = Nxi2(0,a) - Nx(0,a)*xXi2(0,0) - Nx(1,a)*xXi2(1,0);
      B(1,a) = Nxi2(1,a) - Nx(0,a)*xXi2(0,1) - Nx(1,a)*xXi2(1,1);
      B(2,a) = Nxi2(2,a) - Nx(0,a)*xXi2(0,2) - Nx(1,a)*xXi2(1,2);
    }

    // Compute the solution to the linear equations K * X = B.
    //
    Vector<int> IPIV(l);
    int INFO;

    dgesv_(&l, &eNoN, K.data(), &l, IPIV.data(), B.data(), &l, &INFO);

    if (INFO != 0) {
      throw std::runtime_error("[gn_nxx] Error in Lapack");
    }

    Nxx = B;

  } else if (insd == 3) {

    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < nsd; i++) {
        xXi(i,0) = xXi(i,0) + lx(i,a)*Nxi(0,a);
        xXi(i,1) = xXi(i,1) + lx(i,a)*Nxi(1,a);
        xXi(i,2) = xXi(i,2) + lx(i,a)*Nxi(2,a);

        xXi2(i,0) = xXi2(i,0) + lx(i,a)*Nxi2(0,a);
        xXi2(i,1) = xXi2(i,1) + lx(i,a)*Nxi2(1,a);
        xXi2(i,2) = xXi2(i,2) + lx(i,a)*Nxi2(2,a);
        xXi2(i,3) = xXi2(i,3) + lx(i,a)*Nxi2(3,a);
        xXi2(i,4) = xXi2(i,4) + lx(i,a)*Nxi2(4,a);
        xXi2(i,5) = xXi2(i,5) + lx(i,a)*Nxi2(5,a);
      }
    }

    for (int i = 0; i < 3; i++) { 
      K.set_row(i, { xXi(0,i)*xXi(0,i), xXi(1,i)*xXi(1,i), xXi(2,i)*xXi(2,i), 
                     t*xXi(0,i)*xXi(1,i), t*xXi(1,i)*xXi(2,i), t*xXi(0,i)*xXi(2,i) 
                   } );
    }

    int i = 0;
    int j = 1;

    K.set_row(3, { xXi(0,i)*xXi(0,j), xXi(1,i)*xXi(1,j),
                   xXi(2,i)*xXi(2,j),
                   xXi(0,i)*xXi(1,j) + xXi(0,j)*xXi(1,i),
                   xXi(1,i)*xXi(2,j) + xXi(1,j)*xXi(2,i),
                   xXi(0,i)*xXi(2,j) + xXi(0,j)*xXi(2,i) 
                 } );

     i = 1;
     j = 2;
     K.set_row(4, { xXi(0,i)*xXi(0,j), xXi(1,i)*xXi(1,j),
                    xXi(2,i)*xXi(2,j),
                    xXi(0,i)*xXi(1,j) + xXi(0,j)*xXi(1,i),
                    xXi(1,i)*xXi(2,j) + xXi(1,j)*xXi(2,i),
                    xXi(0,i)*xXi(2,j) + xXi(0,j)*xXi(2,i) 
                  } );

     i = 0;
     j = 2;
     K.set_row(5, { xXi(0,i)*xXi(0,j), xXi(1,i)*xXi(1,j),
                    xXi(2,i)*xXi(2,j),
                    xXi(0,i)*xXi(1,j) + xXi(0,j)*xXi(1,i),
                    xXi(1,i)*xXi(2,j) + xXi(1,j)*xXi(2,i),
                    xXi(0,i)*xXi(2,j) + xXi(0,j)*xXi(2,i) 
                  } );


    for (int a = 0; a < eNoN; a++) {
      for (int i = 0; i < 6; i++) {
        B(i,a) = Nxi2(i,a) - Nx(0,a)*xXi2(0,i) - Nx(1,a)*xXi2(1,i) - Nx(2,a)*xXi2(2,i);
      }
    }

    // Compute the solution to the linear equations K * X = B.
    //
    Vector<int> IPIV(l);
    int INFO;

    dgesv_(&l, &eNoN, K.data(), &l, IPIV.data(), B.data(), &l, &INFO);

    if (INFO != 0) {
      throw std::runtime_error("[gn_nxx] Error in Lapack");
    }

    Nxx = B;
  }
}

/// @brief Set mesh properties for the input element type. 
///
/// Mesh data set 
///   mesh % eType - element type (e.g. eType_TET4)
///   mesh % nG - number of element gauss points 
///   mesh % vtkType - element VTK type (e.g. 10 for tet4)
///   mesh % nEf - number of element faces
///   mesh % lShpF - if the basis function is linear
///     
/// Mesh arrays allocated 
///   mesh % w(mesh % nG) - Gauss weights
///   mesh % xi(insd,mesh % nG) - Gauss integration points in parametric space
///   mesh % N(mesh % eNoN,mesh % nG) - Parent shape function
///   mesh % Nx(insd, mesh % eNoN, mesh % nG) - Parent shape functions gradient
///   mesh % xib(2,nsd) - Bounds on Gauss integration points in parametric space
///   mesh % Nb(2,mesh % eNoN) - Bounds on shape functions
//
void select_ele(const ComMod& com_mod, mshType& mesh)
{
  // Set integration dimension: shell, fiber or solid (2d or 3d)
  int insd;
  if (mesh.lShl) {
    insd = com_mod.nsd - 1;
  } else if (mesh.lFib) {
    insd = 1;
  } else {
    insd = com_mod.nsd;
  }

  // Set element properties based on integration dimension 
  // and number of element nodes.
  //
  try {
    if (insd == 3) { 
      set_3d_element_props[mesh.eNoN](insd, mesh);
    } else if (insd == 2) { 
      set_2d_element_props[mesh.eNoN](insd, mesh);
    } else if (insd == 1) { 
      set_1d_element_props[mesh.eNoN](insd, mesh);
    }
  } catch (const std::bad_function_call& exception) {
      throw std::runtime_error("[select_ele] No support for " + std::to_string(mesh.eNoN) + " noded " + 
          std::to_string(insd) + "D elements.");
  }

  // Set mesh 'w' and 'xi' arrays used for Gauss integration.
  mesh.w = Vector<double>(mesh.nG); 
  mesh.xi = Array<double>(insd, mesh.nG); 
  get_gip(mesh);

  // Create mesh 'N' and 'Nx' shape function arrays.
  mesh.N = Array<double>(mesh.eNoN, mesh.nG); 
  mesh.Nx = Array3<double>(insd, mesh.eNoN, mesh.nG); 
  for (int g = 0; g < mesh.nG; g++) {
    get_gnn(g, mesh);
  }

  // Create bounds on Gauss integration points and shape functions.
  mesh.xib = Array<double>(2, com_mod.nsd); 
  mesh.Nb = Array<double>(2, mesh.eNoN); 
  get_nn_bnds(com_mod, mesh);
}

/// @brief Set face properties for the input element type. 
///
/// Face data set 
///   face % eType - element type (e.g. eType_TET4)
//
void select_eleb(Simulation* simulation, mshType& mesh, faceType& face)
{
  // Get the object storing global variables.
  auto& com_mod = simulation->get_com_mod();

  int insd = com_mod.nsd - 1;
  if (mesh.lShl) {
    insd = insd - 1;
   }
  if (mesh.lFib) {
    insd = 0;
  }

  // [NOTE] Not implemented. 
  if (mesh.eType == ElementType::NRB) { 
    face.eType = ElementType::NRB; 
    if (insd == 1) { 
      //ALLOCATE(lFa%Nxx(1,lFa%eNoN,lFa%nG))
    } else if (insd == 2) {
      //ALLOCATE(lFa%Nxx(3,lFa%eNoN,lFa%nG))
    }

  return; 
  }

  // Set element properties based on integration dimension and number of element nodes.
  //
  try {
    set_face_element_props[face.eNoN](insd, face);
  } catch (const std::bad_function_call& exception) {
    throw std::runtime_error("No support for " + std::to_string(face.eNoN) + " noded " +
      std::to_string(insd) + "D elements in 'set_face_element_props'.");
  }

  // Set face 'w' and 'xi' arrays used for Gauss integration.
  face.w = Vector<double>(face.nG);
  face.xi = Array<double>(insd, face.nG);
  get_gip(simulation, face);

  // Create mesh 'N' and 'Nx' shape function arrays.
  face.N = Array<double>(face.eNoN, face.nG); 
  face.Nx = Array3<double>(insd, face.eNoN, face.nG); 
  for (int g = 0; g < face.nG; g++) {
    get_gnn(simulation, g, face);
  }
}

};

