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

#ifndef URIS_H 
#define URIS_H
 
#include "ComMod.h"
#include "Simulation.h"

namespace uris {

void uris_meanp(ComMod& com_mod, CmMod& cm_mod, const int iUris); // done

void uris_meanv(ComMod& com_mod, CmMod& cm_mod, const int iUris); // done

void uris_update_disp(ComMod& com_mod, CmMod& cm_mod);

void uris_find_tetra(ComMod& com_mod, CmMod& cm_mod, const int iUris);

void inside_tet(ComMod& com_mod, int& eNoN, Vector<double>& xp, 
                Array<double>& xl, int& flag, bool ext); // done

void uris_read_msh(Simulation* simulation); // done

void uris_write_vtus(ComMod& com_mod); // done

void uris_calc_sdf(ComMod& com_mod); // done

void uris_read_sv(Simulation* simulation, mshType& mesh, const URISFaceParameters* mesh_param); //done

int in_poly(Vector<double>& P, Array<double>& P1, bool ext); // done

int same_side(Vector<double>& v1, Vector<double>& v2, Vector<double>& v3,
              Vector<double>& v4, Vector<double>& p, bool ext); //done

}

#endif

