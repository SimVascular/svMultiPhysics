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

#ifndef RIS_H 
#define RIS_H

#include "ComMod.h"

namespace ris {

void ris_meanq(ComMod& com_mod, CmMod& cm_mod);
void ris_resbc(ComMod& com_mod, const Array<double>& Yg, const Array<double>& Dg);
void setbc_ris(ComMod& com_mod, const bcType& lBc, const mshType& lM, const faceType& lFa, 
    const Array<double>& Yg, const Array<double>& Dg);

void ris_updater(ComMod& com_mod, CmMod& cm_mod);
void ris_status(ComMod& com_mod, CmMod& cm_mod);

void doassem_ris(ComMod& com_mod, const int d, const Vector<int>& eqN, 
    const Array3<double>& lK, const Array<double>& lR); 

void doassem_velris(ComMod& com_mod, const int d, const Array<int>& eqN, 
    const Array3<double>& lK, const Array<double>& lR);

void clean_r_ris(ComMod& com_mod);
void setbcdir_ris(ComMod& com_mod, Array<double>& lA, Array<double>& lY, Array<double>& lD);

// TODO: RIS 0D code
void ris0d_bc(ComMod& com_mod, CmMod& cm_mod, const Array<double>& Yg, const Array<double>& Dg);
void ris0d_status(ComMod& com_mod, CmMod& cm_mod); //, const Array<double>& Yg, const Array<double>& Dg);

};

#endif

