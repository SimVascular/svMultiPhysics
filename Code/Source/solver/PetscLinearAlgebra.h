// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PETSC_LINEAR_ALGEBRA_H 
#define PETSC_LINEAR_ALGEBRA_H 

#include "LinearAlgebra.h"

/// @brief The PetscLinearAlgebra class implements the LinearAlgebra 
/// interface for the PETSc numerical linear algebra package.
///
class PetscLinearAlgebra : public virtual LinearAlgebra {

  public:
    PetscLinearAlgebra();
    ~PetscLinearAlgebra();
    virtual void alloc(ComMod& com_mod, eqType& lEq);
    virtual void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN, 
        const Array3<double>& lK, const Array<double>& lR);
    virtual void check_options(const consts::PreconditionerType prec_cond_type, const consts::LinearAlgebraType assembly_type);
    virtual void initialize(ComMod& com_mod, eqType& lEq);
    virtual void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    virtual void set_assembly(consts::LinearAlgebraType assembly_type);
    virtual void set_preconditioner(consts::PreconditionerType prec_type);
    virtual void finalize();

  private:
    static std::set<consts::LinearAlgebraType> valid_assemblers;
    void initialize_fsils(ComMod& com_mod, eqType& lEq);
    /// @brief The FsilsLinearAlgebra object used to assemble local element matrices.
    LinearAlgebra* fsils_solver = nullptr;
    // Private class used to hide PETSc implementation details.
    class PetscImpl;
    PetscImpl* impl = nullptr;
    bool use_fsils_assembly = false;
};

#endif

