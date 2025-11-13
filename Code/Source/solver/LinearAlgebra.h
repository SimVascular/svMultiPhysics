// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LINEAR_ALGEBRA_H 
#define LINEAR_ALGEBRA_H 

#include "ComMod.h"
#include "consts.h"

/// @brief The LinearAlgebra class provides an abstract interface to linear algebra 
/// frameworks: FSILS, Trilinos, PETSc, etc.
//
class LinearAlgebra {
  public:
    static const std::map<std::string, consts::LinearAlgebraType> name_to_type;
    static const std::map<consts::LinearAlgebraType, std::string> type_to_name;
    static void check_equation_compatibility(const consts::EquationType eq_phys, 
        const consts::LinearAlgebraType lin_alg_type, const consts::LinearAlgebraType assembly_type);

    LinearAlgebra();
    virtual ~LinearAlgebra() { };
    virtual void alloc(ComMod& com_mod, eqType& lEq) = 0;
    virtual void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN, 
        const Array3<double>& lK, const Array<double>& lR) = 0;
    virtual void check_options(const consts::PreconditionerType prec_cond_type, const consts::LinearAlgebraType assembly_type) = 0;
    virtual void initialize(ComMod& com_mod, eqType& lEq) = 0;
    virtual void set_assembly(consts::LinearAlgebraType assembly_type) = 0;
    virtual void set_preconditioner(consts::PreconditionerType prec_type) = 0;
    virtual void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res) = 0;
    virtual void finalize() = 0;
    
    virtual consts::LinearAlgebraType get_interface_type() { return interface_type; }

    consts::LinearAlgebraType interface_type = consts::LinearAlgebraType::none;
    consts::LinearAlgebraType assembly_type = consts::LinearAlgebraType::none;
    consts::PreconditionerType preconditioner_type = consts::PreconditionerType::PREC_NONE;
};

/// @brief The LinearAlgebraFactory class provides a factory used to create objects derived from LinearAlgebra. 
class LinearAlgebraFactory {
  public:
    static LinearAlgebra* create_interface(consts::LinearAlgebraType interface_type);
};


#endif

