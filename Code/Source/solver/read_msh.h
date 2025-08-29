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

#ifndef READ_MSH_H 
#define READ_MSH_H 

#include "ComMod.h"
#include "Simulation.h"
#include "Vector.h"

#include "utils.h"

#include <string>

namespace read_msh_ns {

  class blkType
  {
    public:
      int n = 0;
      Vector<int> gN;
  };

  void calc_elem_ar(ComMod& com_mod, const CmMod& cm_mod, mshType& lM, bool& rflag);
  void calc_elem_jac(ComMod& com_mod, const CmMod& cm_mod, mshType& lM, bool& rflag);
  void calc_elem_skew(ComMod& com_mod, const CmMod& cm_mod, mshType& lM, bool& rflag);

  void calc_mesh_props(ComMod& com_mod, const CmMod& cm_mod, const int nMesh, std::vector<mshType>& mesh);

  void calc_nbc(mshType& mesh, faceType& face);

  void check_ien(Simulation* simulation, mshType& mesh);
  void check_line_conn(mshType& mesh);
  void check_hex8_conn(mshType& mesh);
  void check_hex20_conn(mshType& mesh);
  void check_hex27_conn(mshType& mesh);
  void check_quad4_conn(mshType& mesh);
  void check_tet_conn(mshType& mesh);
  void check_tri3_conn(mshType& mesh);
  void check_tri6_conn(mshType& mesh);
  void check_wedge_conn(mshType& mesh);

  int find_blk(const int nsd, const int nBkd, const std::vector<bool>& nFlt, const Vector<double>&xMin, const Vector<double>&dx, const Vector<double>& x);

  void load_var_ini(Simulation* simulation, const ComMod& com_mod);

  void match_faces(const ComMod& com_mod, const faceType& face1, const faceType& face2, const double tol, utils::stackType& lPrj);
  void match_nodes(const ComMod& com_mod, const faceType& lFa, const faceType& pFa, 
                   const double ptol, const int nNds, Array<int>& map);

  void read_fib_nff(Simulation* simulation, mshType& mesh, const std::string& fName, const std::string& kwrd, const int idx);
  void read_msh(Simulation* simulation);

  void set_dmn_id_ff(Simulation* simulation, mshType& mesh, const std::string& file_name);
  void set_dmn_id_vtk(Simulation* simulation, mshType& mesh, const std::string& file_name, const std::string& kwrd);
  void set_projector(Simulation* simulation, utils::stackType& avNds);
  void set_ris_projector(Simulation* simulation);
  void set_uris_meshes(Simulation* simulation);
  
};

#endif

