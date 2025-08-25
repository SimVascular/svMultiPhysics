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

#ifndef CMM_H 
#define CMM_H 

#include "ComMod.h"
#include "Simulation.h"

/// @brief These subroutines implement the Coupled Momentum Method (CMM).
namespace cmm {

void cmm_3d(ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, const Array<double>& Nx,  
    const Array<double>& al, const Array<double>& yl, const Array<double>& bfl, const Array<double>& Kxi,     
    Array<double>& lR, Array3<double>& lK);

void cmm_b(ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const int e, const Array<double>& al, const Array<double>& dl, 
    const Array<double>& xl, const Array<double>& bfl, const Vector<double>& pS0l, const Vector<double>& vwp, 
    const Vector<int>& ptr);

void bcmmi(ComMod& com_mod, const int eNoN, const int idof, const double w, const Vector<double>& N, const Array<double>& Nxi,
    const Array<double>& xl, const Array<double>& tfl, Array<double>& lR);

void cmmi(ComMod& com_mod, const mshType& lM, const Array<double>& al, const Array<double>& dl, const Array<double>& xl,
    const Array<double>& bfl, const Array<double>& pS0l, const Vector<double>& vwp, const Vector<int>& ptr);

void cmm_mass(ComMod& com_mod, const double w, const Vector<double>& N, const Array<double>& al, 
    const Array<double>& bfl, const Vector<double>& vwp, Array<double>& lR, Array3<double>& lK);

void cmm_stiffness(ComMod& com_mod, const Array<double>& Nxi, const Array<double>& xl, const Array<double>& dl,                        
    const Vector<double>& pS0l, const Vector<double>& vwp, Vector<double>& pSl, Array<double>& lR, Array3<double>& lK);

void construct_cmm(ComMod& com_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

