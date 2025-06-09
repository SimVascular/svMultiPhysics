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

#include "TrilinosLinearAlgebra.h"
#include <iostream>

// Include Trilinos-dependent data structures and functions. 
//
#ifdef WITH_TRILINOS
#include "trilinos_impl.cpp"

// If Trilinos is not used then define TrilinosImpl with noop methods.
//
#else
class TrilinosLinearAlgebra::TrilinosImpl {
  public:
    TrilinosImpl(){};
    void alloc(ComMod& com_mod, eqType& lEq){};
    void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR){};
    void initialize(ComMod& com_mod) {};
    void set_preconditioner(consts::PreconditionerType prec_type) {};
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
    void solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
};
#endif

/////////////////////////////////////////////////////////////////
//          T r i l i n o s L i n e a r A l g e b r a          //
/////////////////////////////////////////////////////////////////
// The following methods implement the Trilinos LinearAlgebra interface.

/// Trilinos can use fsils or trilinos for assembly.
std::set<consts::LinearAlgebraType> TrilinosLinearAlgebra::valid_assemblers = {
  consts::LinearAlgebraType::fsils,
  consts::LinearAlgebraType::trilinos
};
 
TrilinosLinearAlgebra::TrilinosLinearAlgebra()
{
  #ifndef WITH_TRILINOS
  throw std::runtime_error("[TrilinosLinearAlgebra] svMultiPhysics has not been built with Trilinos.");
  #else
  impl = new TrilinosLinearAlgebra::TrilinosImpl();
  interface_type = consts::LinearAlgebraType::trilinos; 
  assembly_type = consts::LinearAlgebraType::trilinos;
  preconditioner_type = consts::PreconditionerType::PREC_TRILINOS_DIAGONAL;
  #endif
}

TrilinosLinearAlgebra::~TrilinosLinearAlgebra()
{
  if (fsils_solver != nullptr) { 
    delete fsils_solver;
  }
}

/// @brief Allocate data arrays.
void TrilinosLinearAlgebra::alloc(ComMod& com_mod, eqType& lEq)
{
  if (use_fsils_assembly) {
    int dof = com_mod.dof;
    com_mod.Val.resize(dof*dof, com_mod.lhs.nnz);
    initialize_fsils(com_mod, lEq);
  }

  impl->alloc(com_mod, lEq);
}

/// @brief Assemble local element arrays.
///
/// Assembly can be performed using fsils or trilinos.
///
void TrilinosLinearAlgebra::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR)
{
  if (use_fsils_assembly) {
    fsils_solver->assemble(com_mod, num_elem_nodes, eqN, lK, lR);
  } else {
    impl->assemble(com_mod, num_elem_nodes, eqN, lK, lR);
  }
}

/// @brief Check the validity of the precondition and assembly options. 
/// 
/// Trilinos can use fsils or trilinos for assembly.
//
void TrilinosLinearAlgebra::check_options(const consts::PreconditionerType prec_cond_type, 
    const consts::LinearAlgebraType assembly_type)
{
  using namespace consts;
  std::string error_msg;

  if (trilinos_preconditioners.count(prec_cond_type) == 0) {
    auto prec_cond_type_name = consts::preconditioner_type_to_name.at(prec_cond_type);
    error_msg = "trilinos linear algebra can't use '" + prec_cond_type_name + "' for a preconditioner.";
  }

  if ((assembly_type != consts::LinearAlgebraType::none) && 
      (valid_assemblers.count(assembly_type) == 0)) {
    auto assembly_type_name = LinearAlgebra::type_to_name.at(assembly_type);
    error_msg = "trilinos linear algebra can't use '" + assembly_type_name + "' for assembly.";
  }

  if (error_msg != "") {
    throw std::runtime_error("[svMultiPhysics] ERROR: " + error_msg);
  }
}

/// @brief Initialize Trilinos framework.
void TrilinosLinearAlgebra::initialize(ComMod& com_mod, eqType& lEq)
{
  impl->initialize(com_mod);
}

/// @brief Create an fsils linear algebra interface for assembly.
void TrilinosLinearAlgebra::initialize_fsils(ComMod& com_mod, eqType& lEq)
{
  #define n_debug_initialize_fsils
  #ifdef debug_initialize_fsils
  std::cout << "[TrilinosLinearAlgebra::initialize_fsils] ---------- initialize_fsils ---------- " << std::endl;
  std::cout << "[TrilinosLinearAlgebra::initialize_fsils] preconditioner_type: " << preconditioner_type << std::endl;
  #endif

  if (fsils_solver != nullptr) { 
    delete fsils_solver;
  }

  fsils_solver = LinearAlgebraFactory::create_interface(consts::LinearAlgebraType::fsils);
  fsils_solver->initialize(com_mod, lEq);
  fsils_solver->alloc(com_mod, lEq);
  fsils_solver->set_assembly(consts::LinearAlgebraType::fsils);
}

/// @brief Set the linear algebra package for assmbly.
void TrilinosLinearAlgebra::set_assembly(consts::LinearAlgebraType atype)
{
  if (assembly_type == consts::LinearAlgebraType::none) {
    return;
  }

  if (valid_assemblers.count(atype) == 0) {
    auto str_type = LinearAlgebra::type_to_name.at(atype);
    throw std::runtime_error("[TrilinosLinearAlgebra] ERROR: Can't set Trilinos linear algebra to use '" + 
      str_type + "' for assembly." + " Trilinos can only use 'trilinos' for assembly.");
  }

  assembly_type = atype;

  if (assembly_type == consts::LinearAlgebraType::fsils) {
    use_fsils_assembly = true;
  }
}
 
/// @brief Set the proconditioner.
void TrilinosLinearAlgebra::set_preconditioner(consts::PreconditionerType prec_type)
{
  if (consts::trilinos_preconditioners.count(prec_type) == 0) {
    auto str_type = consts::preconditioner_type_to_name.at(prec_type);
    throw std::runtime_error("[TrilinosLinearAlgebra] ERROR: trilinos linear algebra can't use '" +
        str_type + "' for a preconditioner.");
  }

  preconditioner_type = prec_type;
  impl->set_preconditioner(prec_type);
}

/// @brief Solve a system of linear equations.
void TrilinosLinearAlgebra::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  //std::cout << "[TrilinosLinearAlgebra::solve] solve " << std::endl;
  if (assembly_type == consts::LinearAlgebraType::trilinos) {
    impl->solve_assembled(com_mod, lEq, incL, res);
  } else {
    impl->solve(com_mod, lEq, incL, res);
  }
}

