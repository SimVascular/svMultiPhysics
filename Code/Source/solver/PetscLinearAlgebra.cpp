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

#include "PetscLinearAlgebra.h"

#include <iostream>

// Include PETSc-dependent data structures and functions. 
//
#ifdef USE_PETSC
#include "petsc_impl.cpp"

// If PETSc is not used then define PetscImpl with noop methods.
//
#else
class PetscLinearAlgebra::PetscImpl {
  public:
    PetscImpl(){};
    void initialize(ComMod& com_mod, eqType& lEq) {};
    void set_preconditioner(consts::PreconditionerType prec_type) {};
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) {};
};
#endif

/////////////////////////////////////////////////////////////////
//             P e t s c L i n e a r A l g e b r a             //
/////////////////////////////////////////////////////////////////
// The following methods implement the LinearAlgebra interface.

std::set<consts::LinearAlgebraType> PetscLinearAlgebra::valid_assemblers = {
  consts::LinearAlgebraType::none,
  consts::LinearAlgebraType::fsils
};

PetscLinearAlgebra::PetscLinearAlgebra()
{
  #ifndef USE_PETSC
  throw std::runtime_error("[PetscLinearAlgebra] svMultiPhysics hase not been built with the PETSc package.");
  #else
  impl = new PetscLinearAlgebra::PetscImpl();
  interface_type = consts::LinearAlgebraType::petsc; 
  assembly_type = consts::LinearAlgebraType::none;
  preconditioner_type = consts::PreconditionerType::PREC_PETSC_JACOBI;
  #endif
}

PetscLinearAlgebra::~PetscLinearAlgebra()
{
}

/// @brief Allocate data arrays.
void PetscLinearAlgebra::alloc(ComMod& com_mod, eqType& lEq)
{
  initialize_fsils(com_mod, lEq);
}

/// @brief Assemble local element arrays.
void PetscLinearAlgebra::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN, 
    const Array3<double>& lK, const Array<double>& lR)
{
  fsils_solver->assemble(com_mod, num_elem_nodes, eqN, lK, lR);
}

/// @brief Check the validity of the precondition and assembly types options. 
void PetscLinearAlgebra::check_options(const consts::PreconditionerType prec_cond_type, 
    const consts::LinearAlgebraType assembly_type)
{
  using namespace consts;
  auto prec_cond_type_name = consts::preconditioner_type_to_name.at(prec_cond_type);
  auto assembly_type_name = LinearAlgebra::type_to_name.at(assembly_type);
  std::string error_msg;

  if (valid_assemblers.count(assembly_type) == 0) {
    auto assembly_type_name = LinearAlgebra::type_to_name.at(assembly_type);
    error_msg = "petsc linear algebra can't use '" + assembly_type_name + "' for assembly.";
  }

  if (petsc_preconditioners.count(prec_cond_type) == 0) { 
    error_msg = "petsc linear algebra can't use '" + prec_cond_type_name + "' for preconditioning.";
  }

  if (error_msg != "") {
    throw std::runtime_error("[svMultiPhysics] ERROR: " + error_msg);
  }
}

/// @brief Initialize the PETSc framework.
void PetscLinearAlgebra::initialize(ComMod& com_mod, eqType& lEq)
{
  impl->initialize(com_mod, lEq);
}

/// @brief Initialize an FsilsLinearAlgebra object used for assembly.
void PetscLinearAlgebra::initialize_fsils(ComMod& com_mod, eqType& lEq)
{
  fsils_solver = LinearAlgebraFactory::create_interface(consts::LinearAlgebraType::fsils);
  fsils_solver->initialize(com_mod, lEq);
  fsils_solver->alloc(com_mod, lEq);
}

/// @brief Set the linear algebra package for assembly.
void PetscLinearAlgebra::set_assembly(consts::LinearAlgebraType atype)
{
  if (atype == consts::LinearAlgebraType::none) {
    return;
  }

  if (valid_assemblers.count(atype) == 0) {
    auto str_type = LinearAlgebra::type_to_name.at(atype);
    throw std::runtime_error("[PetscLinearAlgebra] ERROR: Can't set petsc linear algebra to use '" +
      str_type + "' for assembly." + " petsc can only use 'fsils' for assembly.");
  }

  assembly_type = atype;

  if (assembly_type == consts::LinearAlgebraType::fsils) {
    use_fsils_assembly = true;
  }
}

/// @brief Set the proconditioner.
void PetscLinearAlgebra::set_preconditioner(consts::PreconditionerType prec_type)
{
  if (consts::petsc_preconditioners.count(prec_type) == 0) {
    auto str_type = consts::preconditioner_type_to_name.at(prec_type);
    throw std::runtime_error("[PetscLinearAlgebra] ERROR: petsc linear algebra can't use '" +
        str_type + "' for a preconditioner.");
  }

  preconditioner_type = prec_type;
  impl->set_preconditioner(prec_type);
}

/// @brief Solve a system of linear equations.
void PetscLinearAlgebra::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  //std::cout << "[PetscLinearAlgebra] solve" << std::endl;
  impl->solve(com_mod, lEq, incL, res);
}

