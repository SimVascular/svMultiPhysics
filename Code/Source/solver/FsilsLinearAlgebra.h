// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef FSILS_LINEAR_ALGEBRA_H 
#define FSILS_LINEAR_ALGEBRA_H 

#include "LinearAlgebra.h"

/// @brief The FsilsLinearAlgebra class implements the LinearAlgebra 
/// interface for the FSILS numerical linear algebra included in svFSIplus.
///
class FsilsLinearAlgebra : public virtual LinearAlgebra {

  public:
    FsilsLinearAlgebra();
    ~FsilsLinearAlgebra() { };

    virtual void alloc(ComMod& com_mod, eqType& lEq);
    virtual void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR);
    virtual void check_options(const consts::PreconditionerType prec_cond_type, const consts::LinearAlgebraType assembly_type);
    virtual void initialize(ComMod& com_mod, eqType& lEq);
    virtual void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    virtual void set_assembly(consts::LinearAlgebraType atype);
    virtual void set_preconditioner(consts::PreconditionerType prec_type);
    virtual void finalize();

  private:
    /// @brief A list of linear algebra interfaces that can be used for assembly.
    static std::set<consts::LinearAlgebraType> valid_assemblers; 
};

#endif

