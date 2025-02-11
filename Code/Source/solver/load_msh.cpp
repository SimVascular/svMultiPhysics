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

// The functions defined here replicate the Fortran functions defined in LOADMSH.f.

#include "load_msh.h"

#include "consts.h"
#include "nn.h"
#include "read_msh.h"
#include "vtk_xml.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <string>
#include <iomanip>
#include <algorithm>
#include <math.h>

namespace load_msh {

    /// @brief Generate hash maps for mesh nodes and elements

class MeshHashMaps {
public:
    MeshHashMaps() {}

    /// @brief Create a hash map for the global mesh nodes
    ///
    /// \code {.cpp}
    /// mesh
    /// nsd
    /// \endcode
    ///
    std::unordered_map<std::string, int> createNodeHashMap(const mshType& mesh, const int& nsd) {
        std::unordered_map<std::string, int> mesh_node_map;
        for (int i = 0; i < mesh.gnNo; i++) {
            std::ostringstream key;
            key << std::scientific << std::setprecision(16);
            for (int j = 0; j < nsd; j++) {
                key << mesh.x(j,i) <<",";
            }
            mesh_node_map[key.str()] = i;
        }
        if (mesh_node_map.size() != mesh.gnNo) {
            throw std::runtime_error("There may be a duplicate nodes within the mesh.");
        }
        return mesh_node_map;
    }

    /// @brief Create a hash map for the global mesh elements
    ///
    /// \code {.cpp}
    /// mesh
    /// \endcode
    ///
    std::unordered_map<std::string, int> createElementHashMap(const mshType& mesh) {
        std::unordered_map<std::string, int> mesh_element_set;
        for (int i = 0; i < mesh.gnEl; i++) {
            for (unsigned int j = 0; j < mesh.ordering.size(); j++) {
                std::vector<int> element_nodes;
                for (unsigned int k = 0; k < mesh.ordering[j].size(); k++) {
                    element_nodes.push_back(mesh.gIEN(mesh.ordering[j][k], i));
                }
                std::sort(element_nodes.begin(), element_nodes.end());
                std::string key = "";
                for (int node: element_nodes) {
                    key += std::to_string(node) + ",";
                }
                mesh_element_set[key] = i;
            }
        }
        return mesh_element_set;
    }
};


#define ndbg_load_msh

/// @brief Read mesh position coordinates and connectivity.
///
/// [NOTE] Not imiplemented.
//
void read_ccne(Simulation* simulation, mshType& mesh, const MeshParameters* mesh_param)
{
  auto mesh_path = mesh_param->mesh_file_path();
  auto mesh_name = mesh_param->name();
  throw std::runtime_error("[read_ccne] read_ccne() is not implemented.");
}

/// @brief Read list of end nodes from a file into the face data structure.
///
/// This sets the face data: 
/// \code {.cpp}
///   face.nNo 
///   face.gN 
///   face.nEl 
///   face.eNoN 
/// \endcode
//
void read_ndnlff(const std::string& file_name, faceType& face)
{
  std::ifstream end_nodes_file;
  end_nodes_file.open (file_name);
  if (!end_nodes_file.is_open()) {
    throw std::runtime_error("Failed to open the end nodes face file '" + file_name + "'.");
  }

  // Read the integer node IDs.
  //
  int node_id;
  std::string str_value;
  std::vector<int> node_ids;

  while (end_nodes_file >> str_value) {
    std::istringstream str_stream(str_value);
    if (!(str_stream >> node_id)) {
      throw std::runtime_error("Incorrect integer value '" + str_value + "' found in end nodes face file '" + file_name + "'."); 
    }
    node_ids.push_back(node_id);
  }

  if (node_ids.size() == 0) { 
    throw std::runtime_error("Failed to read the end nodes face file '" + file_name + "'.");
  }

  // Set the face node IDs.
  //
  face.nNo = node_ids.size();
  face.gN = Vector<int>(face.nNo);
  for (int i = 0; i < face.nNo; i++) { 
    face.gN[i] = node_ids[i] - 1; 
    //std::cout << "[read_ndnlff] a Ac: " << i+1 << " " << node_ids[i] << std::endl;
  }

  // Set the face element properties.
  face.nEl = 0;
  face.eNoN = 1;
}

/// @brief Create data for a mesh.
///
/// Replicates Fortran READSV subroutine defined in LOADMSH.f.
///
///   SUBROUTINE READSV(list, lM)
//
void read_sv(Simulation* simulation, mshType& mesh, const MeshParameters* mesh_param) {
        auto mesh_path = mesh_param->mesh_file_path();
        auto mesh_name = mesh_param->get_name();
#define n_dbg_read_sv
#ifdef dbg_read_sv
        DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
        dmsg.banner();
        dmsg << "Mesh name: " << mesh_name;
        dmsg << "Mesh path: " << mesh_path;
        dmsg << "mesh.lShl: " << mesh.lShl;
#endif
        // Read in volume mesh.
        vtk_xml::read_vtu(mesh_path, mesh);

        // Check that the input number of spatial dimensions is consistent 
        // with the types of elements defined for the simulation mesh.
        //
        int nsd = simulation->com_mod.nsd;
        int elem_dim = consts::element_dimension.at(mesh.eType);
        auto elem_type = consts::element_type_to_string.at(mesh.eType);

        if (mesh.lShl) { 
          if (nsd == 1) {
            throw std::runtime_error("The number of spatial dimensions (" + std::to_string(nsd) + 
                ") is not consistent with the mesh '" + mesh.name + "' which contains shell elements.");
          }

        } else if (!mesh.lFib) { 
          if (elem_dim != nsd) {
            throw std::runtime_error("The number of spatial dimensions (" + std::to_string(nsd) + 
              ") is not consistent with the mesh '" + mesh.name + "' which contains " + elem_type + " elements.");
          }
        }

        // Set mesh element properites for the input element type.
        nn::select_ele(simulation->com_mod, mesh);

        // Check the mesh element node ordering.
        //
        // Note: This may change element node ordering.
        //
        auto &com_mod = simulation->get_com_mod();
        if (com_mod.usePrecomp) {
            vtk_xml::read_precomputed_solution_vtu(com_mod.precompFileName, com_mod.precompFieldName, mesh);
        }
        if (com_mod.ichckIEN) {
            read_msh_ns::check_ien(simulation, mesh);
        }
        // Read face meshes.
        //
        // Creates nodal coordinates and element connectivity data.
        //
        mesh.nFa = mesh_param->face_parameters.size();
        mesh.fa.resize(mesh.nFa);

        if (mesh.lFib && (mesh.nFa > 1)) {
            throw std::runtime_error("There are " + std::to_string(mesh.nFa) + " faces defined for the '" +
                                     mesh.name + "' mesh. Only one face is allowed for a 1D fiber-based mesh.");
        }

        for (int i = 0; i < mesh.nFa; i++) {
            auto face_param = mesh_param->face_parameters[i];
            auto &face = mesh.fa[i];
            face.name = face_param->name();
#ifdef dbg_read_sv
            dmsg << "face.name: " << face.name;
#endif

            face.qmTRI3 = face_param->quadrature_modifier_TRI3();
            if (face.qmTRI3 < (1.0 / 3.0) || face.qmTRI3 > 1.0) {
                throw std::runtime_error("Quadrature_modifier_TRI3 must be in the range [1/3, 1].");
            }
#ifdef dbg_read_sv
            dmsg << "face.qmTRI3: " << face.qmTRI3;
#endif

            if (mesh.lFib) {
                auto face_path = face_param->end_nodes_face_file_path();
#ifdef dbg_read_sv
                dmsg << "Read end nodes face file ... " << " ";
                dmsg << "face_path: " << face_path;
#endif
                if (face_path == "") {
                    throw std::runtime_error("No end nodes face file path provided.");
                }

                read_ndnlff(face_path, face);

            } else {
                auto face_path = face_param->face_file_path();
                vtk_xml::read_vtp(face_path, face);

                // If node IDs were not read then create them.
                if (face.gN.size() == 0) {
                    read_msh_ns::calc_nbc(mesh, face);
                    // Reset the connectivity with the new node IDs?
                    for (int e = 0; e < face.nEl; e++) {
                        for (int a = 0; a < face.eNoN; a++) {
                            int Ac = face.IEN(a, e);
                            Ac = face.gN(Ac);
                            face.IEN(a, e) = Ac;
                        }
                    }
                }
            }
        }
        if (!mesh.lFib) {
          // Create a hash map for nodes and elements.
          MeshHashMaps mesh_hash_maps;
          auto mesh_node_map = mesh_hash_maps.createNodeHashMap(mesh, com_mod.nsd);
          auto mesh_element_set = mesh_hash_maps.createElementHashMap(mesh);

          for (int i = 0; i < mesh.nFa; i++) {
            auto &face = mesh.fa[i];
            if (!mesh.lFib) {
              // Set face global node IDs
              for (int a = 0; a < face.nNo; a++) {
                std::ostringstream key;
                key << std::scientific << std::setprecision(16);
                for (int j = 0; j < com_mod.nsd; j++) {
                  key << face.x(j, a) << ",";
                }
                face.gN(a) = mesh_node_map[key.str()];
              }

              // Set face element IDs
              for (int e = 0; e < face.nEl; e++) {
                std::vector<int> element_nodes;
                for (int a = 0; a < face.eNoN; a++) {
                  int Ac = face.IEN(a, e);
                  Ac = face.gN(Ac);
                  face.IEN(a, e) = Ac;
                  element_nodes.push_back(Ac);
                }

                std::string key = "";
                std::sort(element_nodes.begin(), element_nodes.end());
                for (int a = 0; a < face.eNoN; a++) {
                  key += std::to_string(element_nodes[a]) + ",";
                }

                face.gE(e) = mesh_element_set[key];
              }
            }
          }
        }

        for (int i = 0; i<mesh.nFa; i++){
            auto &face = mesh.fa[i];
            nn::select_eleb(simulation, mesh, face);
        }
    }
};

