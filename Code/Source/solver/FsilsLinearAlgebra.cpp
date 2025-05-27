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

#include "FsilsLinearAlgebra.h"
#include "fsils_api.hpp"
#include "lhsa.h"
#include <iostream>

/////////////////////////////////////////////////////////////////
//             F s i l s L i n e a r A l g e b r a             //
/////////////////////////////////////////////////////////////////
// The following methods implement the FSILS LinearAlgebra interface.

std::set<consts::LinearAlgebraType> FsilsLinearAlgebra::valid_assemblers = {
  consts::LinearAlgebraType::none,
  consts::LinearAlgebraType::fsils,
};

FsilsLinearAlgebra::FsilsLinearAlgebra()
{
  interface_type = consts::LinearAlgebraType::fsils; 
  assembly_type = consts::LinearAlgebraType::fsils; 
  preconditioner_type = consts::PreconditionerType::PREC_FSILS;
}

/// @brief Allocate data arrays.
void FsilsLinearAlgebra::alloc(ComMod& com_mod, eqType& lEq)
{
  #define n_debug_alloc
  #ifdef debug_alloc
  std::cout << "[FsilsLinearAlgebra::alloc] ---------- alloc ---------- " << std::endl;
  #endif
  int dof = com_mod.dof;

  com_mod.Val.resize(dof*dof, com_mod.lhs.nnz);
}

/// @brief Assemble local element arrays.
void FsilsLinearAlgebra::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR)
{
  #define n_debug_assemble
  #ifdef debug_assemble
  std::cout << "[FsilsLinearAlgebra::assemble] ---------- assemble ---------- " << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] num_elem_nodes: " << num_elem_nodes << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] eqN.size(): " << eqN.size() << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] lK.size(): " << lK.size() << std::endl;
  std::cout << "[FsilsLinearAlgebra::assemble] lR.size(): " << lR.size() << std::endl;
  #endif

  lhsa_ns::do_assem(com_mod, num_elem_nodes, eqN, lK, lR);
}

/// @brief Check the validity of the preconditioner and assembly types options. 
void FsilsLinearAlgebra::check_options(const consts::PreconditionerType prec_cond_type, 
  const consts::LinearAlgebraType assembly_type)
{
  using namespace consts;
  auto prec_cond_type_name = consts::preconditioner_type_to_name.at(prec_cond_type);
  auto assembly_type_name = LinearAlgebra::type_to_name.at(assembly_type);
  std::string error_msg;

  if (valid_assemblers.count(assembly_type) == 0) {
    error_msg = "fsils linear algebra can't use '" + assembly_type_name + "' for assembly.";
  }

  if (fsils_preconditioners.count(prec_cond_type) == 0) { 
    error_msg = "fsils linear algebra can't use '" + prec_cond_type_name + "' for a preconditioner.";
  }

  if (error_msg != "") { 
    throw std::runtime_error("[svMultiPhysics] ERROR: " + error_msg);
  }
}

/// @brief Initialize framework.
void FsilsLinearAlgebra::initialize(ComMod& com_mod, eqType& lEq)
{
  // Nothing is needed to initialize FSILS.
}

/// @brief Set the linear algebra package for assmbly.
void FsilsLinearAlgebra::set_assembly(consts::LinearAlgebraType atype)
{
  if (atype == consts::LinearAlgebraType::none) {
    return;
  }

  if (valid_assemblers.count(atype) == 0) {
    auto str_type = LinearAlgebra::type_to_name.at(atype);
    throw std::runtime_error("[FsilsLinearAlgebra] ERROR: Can't set fsils linear algebra to use '" +
      str_type + "' for assembly.");
  }

  assembly_type = atype;
}

/// @brief Set the preconditioner.
void FsilsLinearAlgebra::set_preconditioner(consts::PreconditionerType prec_type)
{
  if (consts::fsils_preconditioners.count(prec_type) == 0) {
    auto prec_cond_type_name = consts::preconditioner_type_to_name.at(prec_type);
    throw std::runtime_error("[FsilsLinearAlgebra] ERROR: fsils linear algebra can't use '" + 
        prec_cond_type_name + "' for a preconditioner.");
    return;
  }

  preconditioner_type = prec_type;
}

/// @brief Solve a system of linear equations.
void FsilsLinearAlgebra::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  auto& lhs = com_mod.lhs;
  int dof = com_mod.dof;
  auto& R = com_mod.R;      
  auto& Val = com_mod.Val;
  auto preconditioner = lEq.linear_algebra_preconditioner;

  fsi_linear_solver::fsils_solve(lhs, lEq.FSILS, dof, R, Val, preconditioner, incL, res);
}

