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

