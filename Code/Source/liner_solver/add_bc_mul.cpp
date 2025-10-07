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

#include "add_bc_mul.h"

#include "dot.h"

namespace add_bc_mul {

/// @brief The contribution of coupled BCs is added to the matrix-vector
/// product operation. Depending on the type of operation (adding the
/// contribution or computing the PC contribution) different
/// coefficients are used.
///
/// For reference, see 
/// Moghadam et al. 2013 eq. 27 (https://doi.org/10.1016/j.jcp.2012.07.035) and
/// Moghadam et al. 2013b (https://doi.org/10.1007/s00466-013-0868-1).
///
/// Reproduces code in ADDBCMUL.f.
/// @param lhs The left-hand side of the linear system. 0D resistance is stored in the face(i).res field.
/// @param op_Type The type of operation (addition or PC contribution)
/// @param dof The number of degrees of freedom.
/// @param X The input vector.
/// @param Y The current matrix-vector product (Y = K*X), to which we add K^BC * X = res * v * v^T * X.
/// The expression is slightly different if preconditioning.
void add_bc_mul(FSILS_lhsType& lhs, const BcopType op_Type, const int dof, const Array<double>& X, Array<double>& Y)
{
  Vector<double> coef(lhs.nFaces); 
  Array<double> v(dof,lhs.nNo);

  //Setting coef depending on adding resistance to stiffness or
  //computing preconditioner

  if (op_Type == BcopType::BCOP_TYPE_ADD) {
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = lhs.face[i].res;
    }
  } else if (op_Type == BcopType::BCOP_TYPE_PRE) {
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = -lhs.face[i].res / (1.0 + (lhs.face[i].res*lhs.face[i].nS));
    }
  } else { 
    //PRINT *, "FSILS: op_Type is not defined"
    //STOP "FSILS: FATAL ERROR"
  }

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    // Cap faces do not contribute to the tangent matrix
    if (face.isCap) {
      continue;
    }

    // In the following calculations, we are computing the product of the
    // coupled BC tangent contribution with the vector X (refer to Moghadam
    // et al. 2013 eq. 27). This is computed by first constructing the vector
    // v, which one of the integrals found in the expression, int{N_A * n_i} dGamma.
    // Then, v is dotted with X to yield a quantity S. Then S is multiplied by
    // by v again, and also multiplied by the appropriate coefficients in
    // the expression.
    // The calculations are complicated somewhat if there is a capping surface,
    // but these complications are explained below.


    // Calculating S, which is the inner product of the right integral (v) and
    // the vector to be multiplied (X).

    int nsd = std::min(face.dof, dof);

    if (face.coupledFlag) {
      // If face is shared across procs
      if (face.sharedFlag) {
        v = 0.0;
        // Setting vector v = int{N_A n_i} dGamma
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }
        // Computing S = coef * v^T * X
        double S = coef(faIn) * dot::fsils_dot_v(dof, lhs.mynNo, lhs.commu, v, X);
        
        // Add cap surface contribution to S
        // Cap surfaces contribute to flow rate but not pressure
        if (face.isCapped) {
          int faInCap = face.faInCap;
          auto& faceCap = lhs.face[faInCap];

          Array<double> vcap(dof,lhs.nNo);
          vcap = 0.0;
          for (int a = 0; a < faceCap.nNo; a++) {
            int Ac = faceCap.glob(a);
            for (int i = 0; i < nsd; i++) {
              vcap(i,Ac) = faceCap.valM(i,a);  
            }
          }
          // Add capping surface contribution to S
          S = S + coef(faIn) * dot::fsils_dot_v(dof, lhs.mynNo, lhs.commu, vcap, X);
        }
        
        // Computing Y = Y + v * S
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            Y(i,Ac) = Y(i,Ac) + v(i,Ac)*S;
          }
        }
      } 
      // If face is not shared across procs
      else  {
        // Computing S = coef * v^T * X
        double S = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            S = S + face.valM(i,a)*X(i,Ac);
          }
        }
        
        // If capping surface is present add its contribution to S
        if (face.isCapped) {
          int faInCap = face.faInCap;
          auto& faceCap = lhs.face[faInCap];

          if (!faceCap.coupledFlag) {
            std::cerr << "ADDBCMUL(): Cap face is not coupled. Probably cap face has zero resistance." << std::endl;
            throw std::runtime_error("FSILS: FATAL ERROR");
          }

          for (int a = 0; a < faceCap.nNo; a++) {
            int Ac = faceCap.glob(a);
            for (int i = 0; i < nsd; i++) {
              S = S + faceCap.valM(i,a)*X(i,Ac);
            }
          }
        }

        // Multiply S by the resistance or related quantity if
        // preconditioning
        S = coef(faIn) * S;

        // Computing Y = Y + v * S
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            Y(i,Ac) = Y(i,Ac) + face.valM(i,a)*S;
          }
        }
      }
    }
  }

}

};
