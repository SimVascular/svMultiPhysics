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

// The functions defined here replicate the Fortran functions defined in VTKXML.f.

#include "vtk_xml.h"
#include "vtk_xml_parser.h"
#include "VtkData.h"

#include "all_fun.h"
#include "consts.h"
#include "post.h"

#include <iomanip>
#include <sstream>
#include <stdio.h>

#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>

namespace vtk_xml {

#define dbg_vtk_xml
#define n_dbg_read_vtu_pdata 

void do_test()
{

  std::string file_name_1 = "x_remesh_restart_0__cm.bin";
  std::string file_name_2 = "x_remesh_restart_1__cm.bin";

  int num_1 = 769;
  Array<double> coords_1(3,num_1);
  coords_1.read(file_name_1);

  int num_2 = 767;
  Array<double> coords_2(3,num_2);
  coords_2.read(file_name_2);

  int num_dupe = 0;
  double tol = 1e-4;

  for (int i = 0; i < num_1; i++) {
    double x1 = coords_1(0,i);
    double y1 = coords_1(1,i);
    double z1 = coords_1(2,i);

    for (int j = 0; j < num_2; j++) {
      double x2 = coords_2(0,j);
      double y2 = coords_2(1,j);
      double z2 = coords_2(2,j);
      double d = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

      if (d < tol) {
        num_dupe += 1;
        break;
      }
    }
  }

  std::cout << "Num dupe: " << num_dupe << std::endl;

  vtkSmartPointer<vtkUnstructuredGrid> mesh;

  std::string fileName = "mesh-complete/mesh-complete.mesh.vtu";
  auto reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName(fileName.c_str());
  reader->Update();
  mesh = reader->GetOutput();

  vtkIdType m_NumPoints = mesh->GetNumberOfPoints();
  vtkIdType m_NumCells = mesh->GetNumberOfCells();
  std::cout << "  Number of points: " << m_NumPoints << std::endl;
  std::cout << "  Number of cells: " << m_NumCells << std::endl;

  auto numPoints = mesh->GetNumberOfPoints();
  auto points = mesh->GetPoints();

  double pt[3];
  int num_found_1 = 0;
  int num_found_2 = 0;

  for (int i = 0; i < numPoints; i++) {
    points->GetPoint(i,pt);
    double x = pt[0];
    double y = pt[1];
    double z = pt[2];

    for (int i = 0; i < num_1; i++) {
      double x1 = coords_1(0,i);
      double y1 = coords_1(1,i);
      double z1 = coords_1(2,i);
      double d = sqrt((x1-x)*(x1-x) + (y1-y)*(y1-y) + (z1-z)*(z1-z));

      if (d < tol) {
        num_found_1 += 1;
        break;
      }
    }

    for (int j = 0; j < num_2; j++) {
      double x2 = coords_2(0,j);
      double y2 = coords_2(1,j);
      double z2 = coords_2(2,j);
      double d = sqrt((x2-x)*(x2-x) + (y2-y)*(y2-y) + (z2-z)*(z2-z));
      if (d < tol) {
        num_found_2 += 1;
        break;
      }
    }
  }

  std::cout << "num_found_1: " << num_found_1 << std::endl;
  std::cout << "num_found_2: " << num_found_2 << std::endl;

  //exit(0);
}

/// @brief This routine prepares data array of a regular mesh
///
/// Parameters:
///
///   d - Stores data for writing to VTK files.
///   outDof - Data dof (i.e. nsd).
///   nOute - Number of element data outputs.
///
/// Modifies:
///
///  d.nNo 
///  d.eNoN 
///  d.IEN 
///  d.xe - Element based variables to be written
///  d.gx - variables after transformation to global format ? 
///  d.x - clears this but may contain nodal coords?
///
/// Reproduces 'SUBROUTINE INTMSHDATA(lM, d, outDof, nOute)'
//
void int_msh_data(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, dataType& d, const int outDof, const int nOute)
{
  auto& cm = com_mod.cm;

  #define n_debug_int_msh_data
  #ifdef debug_int_msh_data
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "nOute: " << nOute;
  dmsg << "outDof: " << outDof;
  dmsg << "d.x.size(): " << d.x.size();
  dmsg << "com_mod.savedOnce: " << com_mod.savedOnce;
  #endif

  d.eNoN = lM.eNoN;
  d.vtkType = lM.vtkType;

  if (cm.mas(cm_mod)) {
    d.nNo  = lM.gnNo;
    d.nEl  = lM.gnEl;
    d.eNoN = lM.eNoN;
  } else { 
    d.nNo  = 0;
    d.nEl  = 0;
  }

  #ifdef debug_int_msh_data
  dmsg << "d.nNo: " << d.nNo;
  dmsg << "d.nEl: " << d.nEl;
  dmsg << "d.eNoN: " << d.eNoN;
  #endif

  d.IEN.resize(d.eNoN,d.nEl); 
  d.gx.resize(outDof,d.nNo);

  d.gx = all_fun::global(com_mod, cm_mod, lM, d.x);
  #ifdef debug_int_msh_data
  dmsg << "d.gx.size(): " << d.gx.size();
  #endif

  d.x.clear();

  if (cm.mas(cm_mod)) {
    for (int e = 0; e < lM.gnEl; e++) {
      int Ec = lM.otnIEN(e);
      d.IEN.set_col(e, lM.gIEN.col(Ec));
    }
  }

  // Element variables
  //
  // 'm' appears to count the number of element output variables.
  //
  int m = nOute;

  if (!com_mod.savedOnce || com_mod.nMsh > 1) {
    if (com_mod.savedOnce) {
      m = m + 1;
    } else { 
      m = m + 2;
    }
  }

  // Add Element based variables to lDe.
  //
  // Note: d.xe may be allocated with size 0.
  //
  Array<double> lDe;

  if (d.xe.size() != 0) {

    if (nOute != 0) {
      lDe.resize(nOute, lM.nEl);
    }

    for (int e = 0; e < lM.nEl; e++) {
      for (int i = 0; i < nOute; i++) {
        lDe(i,e) = d.xe(i,e);
      }
    }
    d.xe.clear();
  }

  // Check if ghost cells are defined.
  //
  if (lM.iGC.size() != 0) {
    m = m + 1;
  }

  if (m != 0) {
    d.xe.resize(d.nEl,m);
  }

  Vector<int> sCount;
  Vector<int> tmpI;
  Vector<int> disps;
  Array<double> gDe;

  m = 0;

  // If files have not been written or there are multiple meshes.
  // 
  if (!com_mod.savedOnce || com_mod.nMsh > 1) {
    #ifdef debug_int_msh_data
    dmsg << "!com_mod.savedOnce || com_mod.nMsh > 1 ";
    dmsg << "com_mod.dmnId.size(): " << com_mod.dmnId.size();
    #endif
    m = 1;

    // If more than one domain.
    //
    // lM.eDist is Element distribution between processors.
    //
    if (com_mod.dmnId.size() != 0) {
      #ifdef debug_int_msh_data
      dmsg << "com_mod.dmnId.size() != 0 ";
      #endif
      sCount.resize(cm.np()); 
      tmpI.resize(d.nEl);
      for (int i = 0; i < cm.np(); i++) {
        sCount(i) = lM.eDist(i+1) - lM.eDist(i);
        #ifdef debug_int_msh_data
        dmsg << "sCount(i): " << sCount(i);
        #endif
      }

      #ifdef debug_int_msh_data
      dmsg << "Send lM.eId.data to tmpI.data ... ";
      dmsg << "  lM.nEl: " << lM.nEl;
      dmsg << "  lM.eId.size(): " << lM.eId.size();
      dmsg << "  tmpI.size(): " << tmpI.size();
      #endif

      MPI_Gatherv(lM.eId.data(), lM.nEl, cm_mod::mpint, tmpI.data(), sCount.data(), lM.eDist.data(), cm_mod::mpint, 
          cm_mod.master, cm.com());

      if (cm.mas(cm_mod)) {
        for (int e = 0; e < lM.gnEl; e++) {
          int Ec = lM.otnIEN(e);
          d.xe(e,m-1) = static_cast<double>(tmpI(Ec));
        }
      }
    } else { 
      for (int j = 0; j < d.xe.nrows(); j++) {
        d.xe(j,m-1) = 1.0;
      }
    }

    if (!com_mod.savedOnce) {
      m = m + 1;
      tmpI.resize(d.nEl);

      if (cm.mas(cm_mod)) {
        if (!cm.seq()) {
          int i = 0;

          for (int e = 0; e < lM.gnEl; e++) {
            if (e+1 > lM.eDist(i)) {
              while (true) { 
                i = i + 1;
                if (e+1 <= lM.eDist(i)) {
                  break;
                }
              }
            }
            tmpI(e) = i;
          }

          sCount.resize(lM.gnEl);
          sCount = tmpI;
          for (int e = 0; e < lM.gnEl; e++) {
            int Ec = lM.otnIEN(e);
            d.xe(e,m-1) = static_cast<double>(sCount(Ec));
          }
        }
      }
    }
  }

  // [NOTE] Skip this code, problems with 0-sized arrays.

  if (0) {
    if (cm.mas(cm_mod)) {
      disps.resize(cm.np()); 
      sCount.resize(cm.np()); 
      gDe.resize(nOute,d.nEl);

      for (int i = 0; i < cm.np(); i++) {
        disps(i) = lM.eDist(i)*nOute;
        sCount(i) = lM.eDist(i+1)*nOute - disps(i);
        #ifdef debug_int_msh_data
        dmsg << "disps(i): " << disps(i);
        dmsg << "sCount(i): " << sCount(i);
        #endif
      }
    } else { 
      disps.clear(); 
      sCount.clear(); 
      gDe.clear();
    }

    #ifdef debug_int_msh_data
    dmsg << "Send lDe to gDe ...";
    dmsg << "  lDe.size(): " << lDe.size();
    dmsg << "  nOute*lM.nEl: " << nOute*lM.nEl;
    dmsg << "  gDe.size(): " << gDe.size();
    #endif
    MPI_Gatherv(lDe.data(), nOute*lM.nEl, cm_mod::mpreal, gDe.data(), sCount.data(), disps.data(), cm_mod::mpreal, cm_mod.master, cm.com());

    if (cm.mas(cm_mod)) {
      for (int e = 0; e < lM.gnEl; e++) {
        int Ec = lM.otnIEN(e);
        for (int i = 0; i < nOute; i++) {
          d.xe(e,m+i) = gDe(i,Ec);
        }
      }
    }
    m = m + nOute;
  }

  // [NOTE] Skip this code, problems with 0-sized arrays.
  //
  if (0) {
  //if (lM.iGC.size() != 0) {
    #ifdef debug_int_msh_data
    dmsg;
    dmsg << "lM.iGC.size() != 0 ";
    #endif
    m = m + 1;
    sCount.resize(cm.np()); 
    tmpI.resize(d.nEl);
    for (int i = 0; i < cm.np(); i++) {
      sCount(i) = lM.eDist(i+1) - lM.eDist(i);
    }

    #ifdef debug_int_msh_data
    dmsg << "Send iGC to tmpI .. ";
    dmsg << "  lM.nEl: " << lM.nEl;
    dmsg << "  lM.iGC.size(): " << lM.iGC.size();
    dmsg << "  tmpI.size(): " << tmpI.size();
    #endif
    MPI_Gatherv(lM.iGC.data(), lM.nEl, cm_mod::mpint, tmpI.data(), sCount.data(), lM.eDist.data(), cm_mod::mpint, cm_mod.master, cm.com());

    if (cm.mas(cm_mod)) {
      #ifdef debug_int_msh_data
      dmsg << "2 m: " << m;
      dmsg << "lM.gnEl: " << lM.gnEl;
      dmsg << "5 Set d.xe m " << m;
      #endif
      for (int e = 0; e < lM.gnEl; e++) {
        int Ec = lM.otnIEN(e);
	d.xe(e,m-1) = static_cast<double>(tmpI(Ec));
      }
    }
  }
}

//----------
// read_vtp
//----------
// Read a face from a SimVascular .vtp file.
//
// Face variables set
//   face.eNoN - number of noders per element
//   face.gebc - EBC array (gE + gIEN) 
//   face.gnEl - globel number of elements
//   face.nEl - number of elements
//   face.nNo - number of nodes
//   face.x - node coordinates
//
// Replicates Fortran READVTP subroutine defined in VTKXML.f.
//
void read_vtp(const std::string& file_name, faceType& face)
{
  using namespace vtk_xml_parser;

  if (FILE *file = fopen(file_name.c_str(), "r")) {
    fclose(file);
  } else {
    throw std::runtime_error("The VTK face file '" + file_name + "' can't be read.");
  }

  // Check if the face mesh needs to be read from a VTP or VTU file.
  //
  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);

  if (file_ext == "vtp") {
    vtk_xml_parser::load_vtp(file_name, face);

  } else if (file_ext == "vtu") {
    vtk_xml_parser::load_vtu(file_name, face);
  }

  if (face.gN.size() == 0) {
    std::cout << "[WARNING] No node IDs found in the '" << file_name << "' face file.";
  }
  //else {
  //  for (int e = 0; e < face.nEl; e++) {
  //    for (int a = 0; a < face.eNoN; a++) {
        //int Ac = face.IEN(a,e);
        //Ac = face.gN(Ac);
        //face.IEN(a,e) = Ac;
  //    }
  //  }
  //}

  // Create essential BC array.
  //
  if (face.gE.size() == 0) {
    std::cout << "[WARNING] No element IDs found in the '" << file_name << "' face file.";
  } else {
    face.gnEl = face.nEl;
    face.gebc = Array<int>(face.eNoN+1, face.gnEl);

    for (int i = 0; i < face.gE.size(); i++) {
      face.gebc(0,i) = face.gE(i);
      //std::cout << "[read_vtp] i: " << i+1 << "  gebc: " << face.gebc(0,i) << " "; 
      for (int j = 0; j < face.eNoN; j++) {
        face.gebc(j+1,i) = face.IEN(j,i);
        //std::cout << face.gebc(j+1,i) << " "; 
      }
      //std::cout << std::endl; 
    }
  }
}

//----------------
// read_vtp_pdata
//----------------
// Read prestress data from a vtp file.
//
// Face variables set:
//
//   face.x - 
//
// Arguments:
//   fName - The name of the vtu file to read.
//   kwrd - The name of the vtk data array to read.
//   m - The number of data components (?) 
//
// Reproduces Fortran 'SUBROUTINE READVTPPDATA(lFa, fName, kwrd, m, idx)'.
//
void read_vtp_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, 
    const int idx, faceType& face)
{
  if (FILE *file = fopen(fName.c_str(), "r")) {
      fclose(file);
  } else {
    throw std::runtime_error("The VTK VTP data file '" + fName + "' can't be read.");
  }

  // Read the vtk file.
  auto vtk_data = VtkData::create_reader(fName);
  int num_elems = vtk_data->num_elems();
  int num_points = vtk_data->num_points();

  if (num_points != face.nNo) {
    throw std::runtime_error("The number of nodes (" + std::to_string(num_points) +
        ") in the prestress VTK file '" + fName + "' is not equal to the number of nodes ("
        + std::to_string(face.nNo) + ") for the face named '" + face.name + "'.");
  }

  // Check that the vtk file has prestress data.
  if (!vtk_data->has_point_data(kwrd)) { 
    throw std::runtime_error("No PointData DataArray named '" + kwrd + 
        "' found in the prestress VTK file '" + fName + "' for the '" + face.name + "' face.");
  }

  if (m == nsd) {
    // Set the stress data.
    Array<double> tmpR(consts::maxNSD, face.nNo);
    vtk_data->copy_point_data(kwrd, tmpR);

    for (int a = 0; a < face.nNo; a++) {
      for (int i = idx*nsd; i < (idx+1)*nsd; i++) {
        face.x(i, a) = tmpR(i, a);
      }
    }

  // Copy data directly into face.x.
  //
  } else {
    vtk_data->copy_point_data(kwrd, face.x);
  }
}

//----------
// read_vtu
//----------
// Read a mesh from a SimVascular .vtu or .vtp file.
// 
// Mesh variables set
//   mesh.gnNo - number of nodes
//   mesh.gnEl - number of elements
//   mesh.eNoN - number of noders per element
//   mesh.x - node coordinates
//   mesh.gIEN - element connectivity
//   mesh.eType - element type
//
// Replicates Fortran READVTU subroutine defined in VTKXML.f.
//
//   SUBROUTINE READVTU(lM, fName) 
//
void read_vtu(const std::string& file_name, mshType& mesh)
{
  using namespace vtk_xml_parser;

  if (FILE *file = fopen(file_name.c_str(), "r")) {
    fclose(file);
  } else {
    throw std::runtime_error("The VTU mesh file '" + file_name + "' can't be read.");
  }

  // Read data from a VTK file.
  //
  #define n_read_vtu_use_VtkData 
  #ifdef read_vtu_use_VtkData 
  auto vtk_data = VtkData::create_reader(file_name);
  int num_elems = vtk_data->num_elems(); 
  int np_elem = vtk_data->np_elem(); 
  int elem_type = vtk_data->elem_type(); 

  // Set mesh data.
  mesh.nEl = num_elems;
  mesh.eNoN = np_elem;
  mesh.IEN = vtk_data->get_connectivity();
  mesh.x = vtk_data->get_points();

  delete vtk_data;

  #else

  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);
  if (file_ext == "vtp") {
    vtk_xml_parser::load_vtp(file_name, mesh);
  } else if (file_ext == "vtu") {
    vtk_xml_parser::load_vtu(file_name, mesh);
  }
  #endif
}

//----------
// read_precomputed_solution_vtu
//----------

/// @brief Read a mesh file with state-variable solution fields from a .vtu or .vtp file.
///
///  Mesh variables set
///   mesh.Ys = precomputed state-variable solutions (e.g. velocity, pressure, etc.)

void read_precomputed_solution_vtu(const std::string& file_name, const std::string& field_name, mshType& mesh)
{
  using namespace vtk_xml_parser;

  if (FILE *file = fopen(file_name.c_str(), "r")) {
    fclose(file);
  } else {
    throw std::runtime_error("The VTU mesh file '" + file_name + "' can't be read.");
  }

  // Read data from a VTK file.
  //
  #define n_read_vtu_use_VtkData
  #ifdef read_vtu_use_VtkData
  auto vtk_data = VtkData::create_reader(file_name);
  int num_elems = vtk_data->num_elems();
  int np_elem = vtk_data->np_elem();

  // Set mesh data.
  mesh.nEl = num_elems;
  mesh.eNoN = np_elem;
  mesh.IEN = vtk_data->get_connectivity();
  mesh.x = vtk_data->get_points();

  delete vtk_data;

  #else

  auto file_ext = file_name.substr(file_name.find_last_of(".") + 1);
  if (file_ext == "vtp") {
    vtk_xml_parser::load_time_varying_field_vtu(file_name, field_name, mesh);
  } else if (file_ext == "vtu") {
    vtk_xml_parser::load_time_varying_field_vtu(file_name, field_name, mesh);
  }
  #endif
}

//----------------
// read_vtu_pdata
//----------------
// Read prestress data from a vtu file.
//
// Mesh variables set:
//
//   mesh.x - 
//
// Arguments:
//   fName - The name of the vtu file to read.
//   kwrd - The name of the vtk data array to read.
//   nsd - The number of space dimensions.
//   m - The number of data components (?) 
//
void read_vtu_pdata(const std::string& fName, const std::string& kwrd, const int nsd, const int m, const int idx, mshType& mesh)
{
  if (FILE *file = fopen(fName.c_str(), "r")) {
      fclose(file);
  } else {
    throw std::runtime_error("The VTK VTU pressure data file '" + fName + "' can't be read.");
  }

  // Read the vtk file.
  //
  auto vtk_data = VtkData::create_reader(fName);
  int num_elems = vtk_data->num_elems();
  int num_points = vtk_data->num_points();

  if (num_points != mesh.gnNo) {
    throw std::runtime_error("The number of nodes (" + std::to_string(num_points) +
        ") in the prestress VTK file '" + fName + "' is not equal to the number of nodes ("
        + std::to_string(mesh.gnNo) + ") for the mesh named '" + mesh.name + "'.");
  }

  // Check that the vtk file has prestress data.
  if (!vtk_data->has_point_data(kwrd)) { 
    throw std::runtime_error("No PointData DataArray named '" + kwrd + "' found in the prestress VTK file '" + fName + 
        "' for the '" + mesh.name + "' mesh.");
  }

  if (m == nsd) {
    // Set the stress data.
    Array<double> tmpR(consts::maxNSD, mesh.gnNo);
    vtk_data->copy_point_data(kwrd, tmpR);

    for (int a = 0; a < mesh.gnNo; a++) {
      for (int i = idx*nsd; i < (idx+1)*nsd; i++) {
        mesh.x(i, a) = tmpR(i, a);
      }
    }

  // Copy data directly into mesh.x.
  //
  } else {
    vtk_data->copy_point_data(kwrd, mesh.x);
  }
}

//-----------
// read_vtus
//-----------
// Reproduces 'SUBROUTINE READVTUS(lA, lY, lD, fName)'
//
void read_vtus(Simulation* simulation, Array<double>& lA, Array<double>& lY, Array<double>& lD, const std::string& fName)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& cep_mod = simulation->cep_mod;
  const int gtnNo = com_mod.gtnNo;

  if (FILE *file = fopen(fName.c_str(), "r")) {
    fclose(file);
  } else {
    throw std::runtime_error("The VTK VTU data file '" + fName + "' can't be read.");
  }

  // Read the vtk file.
  auto vtk_data = VtkData::create_reader(fName);
  int nNo = vtk_data->num_points();

  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];

    for (int iOut = 0; iOut < eq.nOutput; iOut++) {
      auto& output = eq.output[iOut];

      if (!output.options.spatial) {
        continue;
      }

      int l = output.l;
      int s = eq.s + output.o;
      int e = s + l - 1;
      auto varName = output.name;
      auto oGrp = output.grp;

      Array<double> tmpGS;

      switch (oGrp) {
        case OutputNameType::outGrp_A: 
        case OutputNameType::outGrp_Y:
        case OutputNameType::outGrp_D:
          if (l > 1) {
            tmpGS.resize(maxNSD,nNo);
          } else {
            tmpGS.resize(1,nNo);
          }

          vtk_data->copy_point_data(varName, tmpGS);
          Array<double> tmpS;

          if (nNo != gtnNo) {
            if (l > 1) {
              tmpS.resize(maxNSD,gtnNo);
            } else {
              tmpS.resize(1,gtnNo);
            }

           int b = 0;

           for (int iM = 0; iM < com_mod.nMsh; iM++) {
             auto& msh = com_mod.msh[iM];

             for (int a = 0; a < msh.gnNo; a++) {
               int Ac = msh.gpN(a);
               for (int i = 0; i < maxNSD; i++) {
                 tmpS(i,Ac) = tmpGS(i,b);
               }
               b = b + 1;
             }
           }

           tmpGS = tmpS;
         }
        break;
      }

      switch (oGrp) {
        case OutputNameType::outGrp_A:
          for (int i = 0; i < l; i++) {
            for (int j = 0; j < gtnNo; j++) {
              lA(i+s,j) = tmpGS(i,j);
            }
          }
        break;

        case OutputNameType::outGrp_Y:
          for (int i = 0; i < l; i++) {
            for (int j = 0; j < gtnNo; j++) {
              lY(i+s,j) = tmpGS(i,j);
            }
          }
        break;

        case OutputNameType::outGrp_D:
          for (int i = 0; i < l; i++) {
            for (int j = 0; j < gtnNo; j++) {
              lD(i+s,j) = tmpGS(i,j);
            }
          }
        break;
      }
    }
  }
}

//-----------
// write_vtp
//-----------
// Reproduces Fortran 'SUBROUTINE WRITEVTP(lFa, fName)'
//
void write_vtp(ComMod& com_mod, faceType& lFa, const std::string& fName)
{
  const int nsd = com_mod.nsd;
  //std::cout << "[write_vtp] ========== write_vtp ==========" << std::endl;
  //std::cout << "[write_vtp] lFa.x.size(): " << lFa.x.size() << std::endl;
  //std::cout << "[write_vtp] lFa.IEN.size(): " << lFa.IEN.size() << std::endl;

  auto vtk_writer = VtkData::create_writer(fName);
  vtk_writer->set_points(lFa.x);
  vtk_writer->set_connectivity(nsd, lFa.IEN);

  if (lFa.gN.size() != 0) {
    vtk_writer->set_point_data("GlobalNodeID", lFa.gN);
  }

  if (lFa.gE.size() != 0) {
    vtk_writer->set_point_data("GlobalElementID", lFa.gE);
  }

  vtk_writer->write();
  delete vtk_writer;
}

//-----------
// write_vtu
//-----------
// Reproduces Fortran 'SUBROUTINE WRITEVTU(lM, fName)'
//
void write_vtu(ComMod& com_mod, mshType& lM, const std::string& fName)
{
  const int nsd = com_mod.nsd;
  //std::cout << "[write_vtu] ========== write_vtu ==========" << std::endl;
  //std::cout << "[write_vtu] fName: " << fName << std::endl;

  auto vtk_writer = VtkData::create_writer(fName);

  vtk_writer->set_points(lM.x);
  vtk_writer->set_connectivity(nsd, lM.gIEN);
  vtk_writer->write();

  delete vtk_writer;
}

//-----------------
// write_vtu_debug
//-----------------
// This function can be called withing the code for distributed meshes for debugging.
//
void write_vtu_debug(ComMod& com_mod, mshType& lM, const std::string& fName)
{
  const int nsd = com_mod.nsd;
  //std::cout << "[write_vtu_debug] ========== write_vtu_debug ==========" << std::endl;
  //std::cout << "[write_vtu_debug] fName: " << fName << std::endl;
  //std::cout << "[write_vtu_debug] nsd: " << nsd << std::endl;
  //std::cout << "[write_vtu_debug] lM.nNo: " << lM.nNo << std::endl;
  //std::cout << "[write_vtu_debug] lM.gIEN nrows: " << lM.gIEN.nrows() << std::endl;
  //std::cout << "[write_vtu_debug]         ncols: " << lM.gIEN.ncols() << std::endl;
  //std::cout << "[write_vtu_debug] lM.IEN nrows: " << lM.IEN.nrows() << std::endl;
  //std::cout << "[write_vtu_debug]        ncols: " << lM.IEN.ncols() << std::endl;
  //std::cout << "[write_vtu_debug] lM.gN size: " << lM.gN.size() << std::endl;

  Array<double> x(nsd, lM.nNo);

  for (int a = 0; a < lM.nNo; a++) {
    int Ac = lM.gN(a);
    for (int i = 0; i < nsd; i++) {
      x(i,a) = com_mod.x(i,Ac);
    }
  }

  auto vtk_writer = VtkData::create_writer(fName);

  vtk_writer->set_points(x);

  vtk_writer->set_connectivity(nsd, lM.IEN);

  vtk_writer->write();

  delete vtk_writer;
}

//------------
// write_vtus
//------------
// Reproduces 'SUBROUTINE WRITEVTUS(lA, lY, lD, lAve)'
//
void write_vtus(Simulation* simulation, const Array<double>& lA, const Array<double>& lY, const Array<double>& lD, const bool lAve)
{
  #define n_debug_write_vtus
  #ifdef debug_write_vtus 
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm = com_mod.cm;
  auto& cm_mod = simulation->cm_mod;
  auto& cep_mod = simulation->cep_mod;

  bool lIbl = false;
  bool lD0  = false;

  int a = com_mod.iblank.sum();
  a = cm.reduce(cm_mod, a);

  if (a > 0) {
    lIbl = true;
  }

  const int nsd = com_mod.nsd;
  const int nEq = com_mod.nEq;
  const auto& eqs = com_mod.eq;
  const auto& Dinit = com_mod.Dinit;

  for (int iEq = 0; iEq < nEq; iEq++) {
    if (eqs[iEq].phys == EquationType::phys_CMM && Dinit.size() != 0) {
      lD0 = true;
      break; 
     }
  }

  const int nMsh = com_mod.nMsh;
  const auto& meshes = com_mod.msh;

  int nOut = 1;
  int nOute  = 0;
  int outDof = nsd;
  int nFn = 0;

  for (int iEq = 0; iEq < nEq; iEq++) {
    auto& eq = eqs[iEq];

    for (int iOut = 0; iOut < eq.nOutput; iOut++) {
      if (!eq.output[iOut].options.spatial) {
        continue; 
      }

      auto oGrp = eq.output[iOut].grp;
      #ifdef debug_write_vtus 
      dmsg << "oGrp: " << oGrp;
      #endif

      if (oGrp == OutputNameType::outGrp_fN) {
        nFn = 1;
        for (int iM = 0; iM < nMsh; iM++) {
          nFn = std::max(nFn, meshes[iM].nFn);
        }
        nOut = nOut + nFn;
        outDof = outDof + eq.output[iOut].l * nFn;
      } else { 
        nOut = nOut + 1;
        outDof = outDof + eq.output[iOut].l;
      }

      if (oGrp == OutputNameType::outGrp_J || oGrp == OutputNameType::outGrp_mises ||
          oGrp == OutputNameType::outGrp_I1) { 
        nOute = nOute + 1;
      }
    }
  }

  // iblank array for immersed bodies
  if (lIbl) {
    nOut = nOut + 1;
    outDof = outDof + 1;
  }

  // Initial displacements for CMM equation
  if (lD0) {
    nOut = nOut + 1;
    outDof = outDof + nsd;
  }

  std::vector<std::string> outNames(nOut); 
  std::vector<int> outS(nOut+1); 
  std::vector<std::string>outNamesE(nOute);

  // Prepare all solultions in to dataType d
  //
  std::vector<dataType> d(nMsh);

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = meshes[iM];
    d[iM].x.resize(outDof, msh.nNo);
    d[iM].xe.resize(nOute,msh.nEl);
    Array<double> tmpV(consts::maxNSD,msh.nNo);

    int cOut = 0;
    outS[cOut] = 0;
    outS[cOut+1] = nsd;
    outNames[cOut] = "";

    for (int a = 0; a < msh.nNo; a++) {
      int Ac = msh.gN(a);
      for (int i = 0; i < nsd; i++) {
        d[iM].x(i,a) = com_mod.x(i,Ac) / msh.scF;
      }
    }

    if (lD0) {
      cOut = cOut + 1;
      int is = outS[cOut];
      int ie = is + nsd - 1;
      outS[cOut+1] = ie + 1;
      outNames[cOut] = "Initial_displacement";

      for (int a = 0; a < msh.nNo; a++) {
        int Ac = msh.gN(a);
        for (int i = 0; i < nsd; i++) {
          d[iM].x(i+is,a) = Dinit(i,Ac);
        }
      }
    }

    nOute = 0;
    std::fill(outNamesE.begin(), outNamesE.end(), "");

    for (int iEq = 0; iEq < nEq; iEq++) {
      auto& eq = eqs[iEq];

      for (int iOut = 0; iOut < eq.nOutput; iOut++) {
        if (!eq.output[iOut].options.spatial) {
          continue;
        }

        int l = eq.output[iOut].l;
        int s = eq.s + eq.output[iOut].o;
        int e = s + l - 1;

        cOut = cOut + 1;
        int is = outS[cOut];
        int ie = is + l - 1;
        outS[cOut+1] = ie + 1;

        auto oGrp = eq.output[iOut].grp;
        outNames[cOut] = eq.output[iOut].name;
        Vector<double> tmpVe;

        switch (oGrp) {
          case OutputNameType::outGrp_NA:
            throw std::runtime_error("Undefined output grp in VTK");
          break;

          case OutputNameType::outGrp_A:
            for (int a = 0; a < msh.nNo; a++) {
              int Ac = msh.gN(a);
              for (int i = 0; i < l; i++) {
                 d[iM].x(i+is,a) = lA(i+s,Ac);
              }
            }
          break;

          case OutputNameType::outGrp_Y:
            if (eq.phys != EquationType::phys_heatF) {
                for (int a = 0; a < msh.nNo; a++) {
                    int Ac = msh.gN(a);
                    for (int i = 0; i < l; i++) {
                        d[iM].x(i + is, a) = lY(i + s, Ac);
                    }
                }
            } else {
                if (eq.output[iOut].name == "Velocity") {
                    for (int a = 0; a < msh.nNo; a++) {
                        int Ac = msh.gN(a);
                        for (int i = 0; i < l; i++) {
                            d[iM].x(i + is, a) = lY(i, Ac);
                        }
                    }
                } else {
                    for (int a = 0; a < msh.nNo; a++) {
                        int Ac = msh.gN(a);
                        for (int i = 0; i < l; i++) {
                            d[iM].x(i + is, a) = lY(i + s, Ac);
                        }
                    }
                }
            }
          break;

          case OutputNameType::outGrp_D:
            #ifdef debug_write_vtus 
            dmsg << "case " << " outGrp_D ";
            dmsg << "is: " << is;
            dmsg << "ie: " << ie;
            dmsg << "s: " << s;
            dmsg << "e: " << e;
            #endif
            for (int a = 0; a < msh.nNo; a++) {
              int Ac = msh.gN(a);
              // Divide by mesh scale factor
              for (int i = 0; i < l; i++) {
                 d[iM].x(i+is,a) = lD(i+s,Ac) / msh.scF;
              }
            }
          break;

          case OutputNameType::outGrp_WSS:
          case OutputNameType::outGrp_trac:
            post::bpost(simulation, msh, tmpV, lY, lD, oGrp);

            for (int a = 0; a < msh.nNo; a++) {
              for (int i = 0; i < l; i++) {
                 d[iM].x(i+is,a) = tmpV(i,a); 
              }
            }
          break;

          case OutputNameType::outGrp_vort: 
          case OutputNameType::outGrp_eFlx: 
          case OutputNameType::outGrp_hFlx: 
          case OutputNameType::outGrp_stInv: 
          case OutputNameType::outGrp_vortex: 
          case OutputNameType::outGrp_Visc: 
            post::post(simulation, msh, tmpV, lY, lD, oGrp, iEq);
            for (int a = 0; a < msh.nNo; a++) {
              int Ac = msh.gN(a);
              for (int i = 0; i < l; i++) {
                 d[iM].x(i+is,a) = tmpV(i,a); 
              }
            }
          break;

          case OutputNameType::outGrp_absV: 
            for (int a = 0; a < msh.nNo; a++) {
              int Ac = msh.gN(a);
              for (int i = 0; i < l; i++) {
                d[iM].x(i+is,a) = lY(i,Ac) - lY(i+nsd+1,Ac);
              }
            }
          break;

          case OutputNameType::outGrp_fN:
            cOut = cOut - 1;
            tmpV.resize(nFn*nsd,msh.nNo);
            if (msh.nFn != 0) {
              post::fib_dir_post(simulation, msh, nFn, tmpV, lD, iEq);
            }
            for (int iFn = 0; iFn < nFn; iFn++) {
              is = outS[cOut];
              ie = is + l - 1;
              outS[cOut+1] = ie + 1;
              outNames[cOut] = eq.output[iOut].name + std::to_string(iFn);
              cOut = cOut + 1;

              for (int a = 0; a < msh.nNo; a++) {
                for (int i = 0; i < l; i++) {
                  d[iM].x(i+is,a) = tmpV(i+iFn*nsd,a);
                }
              }
            }
            tmpV.resize(consts::maxNSD, msh.nNo);
          break;

          case OutputNameType::outGrp_fA:
            tmpV.resize(1,msh.nNo);
            if (msh.nFn == 2) {
              post::fib_algn_post(simulation, msh, tmpV, lD, iEq);
            }
            for (int a = 0; a < msh.nNo; a++) {
              d[iM].x(is,a) = tmpV(0,a);
            }
            tmpV.resize(consts::maxNSD, msh.nNo);
          break;

          case OutputNameType::outGrp_stress:
          case OutputNameType::outGrp_cauchy:
          case OutputNameType::outGrp_mises:
            #ifdef debug_write_vtus 
            dmsg << "case " << " outGrp_stress";
            #endif
            tmpV.resize(l,msh.nNo); 
            tmpVe.resize(msh.nEl);

            if (com_mod.pstEq) {
              for (int a = 0; a < msh.nNo; a++) {
                int Ac = msh.gN(a);
                for (int i = 0; i < l; i++) {
                  tmpV(i,a) = com_mod.pS0(i,Ac);
                }
              }
            }

            if (msh.lShl) {
              post::shl_post(simulation, msh, l, tmpV, tmpVe, lD, iEq, oGrp);
              //CALL SHLPOST(msh(iM), l, tmpV, tmpVe, lD, iEq,oGrp)
            } else { 
              if (!com_mod.cmmInit) {
                post::tpost(simulation, msh, l, tmpV, tmpVe, lD, lY, iEq, oGrp);
              }
            }

            for (int a = 0; a < msh.nNo; a++) {
              for (int i = 0; i < l; i++) {
                d[iM].x(i+is,a) = tmpV(i,a);
              }
            }

            if (oGrp == OutputNameType::outGrp_mises) {
              outNamesE[nOute] = "E_VonMises";
              for (int a = 0; a < msh.nEl; a++) {
                d[iM].xe(nOute,a) = tmpVe(a);
              }
              nOute = nOute + 1;
            }

            tmpV.resize(consts::maxNSD,msh.nNo);
          break;

          case OutputNameType::outGrp_J:
          case OutputNameType::outGrp_F:
          case OutputNameType::outGrp_strain:
          case OutputNameType::outGrp_fS:
          case OutputNameType::outGrp_I1:
            #ifdef debug_write_vtus 
            dmsg << "case " << " outGrp_J";
            #endif
            tmpV.resize(l,msh.nNo); 
            tmpVe.resize(msh.nEl);

            if (msh.lShl) {
              post::shl_post(simulation, msh, l, tmpV, tmpVe, lD, iEq, oGrp);
              //CALL SHLPOST(msh(iM), l, tmpV, tmpVe, lD, iEq,oGrp)
            } else {
              post::tpost(simulation, msh, l, tmpV, tmpVe, lD, lY, iEq, oGrp);
            }

            for (int a = 0; a < msh.nNo; a++) {
              for (int i = 0; i < l; i++) {
                d[iM].x(i+is,a) = tmpV(i,a);
              }
            }

            if (oGrp == OutputNameType::outGrp_J) {
              outNamesE[nOute] = "E_Jacobian";
              for (int a = 0; a < msh.nEl; a++) {
                d[iM].xe(nOute,a) = tmpVe(a);
              }
              nOute = nOute + 1;
            }

            if (oGrp == OutputNameType::outGrp_I1) {
              outNamesE[nOute] = "E_CG_I1";
              nOute = nOute + 1;
              for (int a = 0; a < msh.nEl; a++) {
                d[iM].xe(nOute,a) = tmpVe(a);
              }
            }

            tmpV.resize(consts::maxNSD,msh.nNo);
          break;

          case OutputNameType::outGrp_divV:
            tmpV.resize(l,msh.nNo); 
            post::div_post(simulation, msh, tmpV, lY, lD, iEq);
            for (int a = 0; a < msh.nNo; a++) {
              d[iM].x(is,a) = tmpV(0,a);
            }
            tmpV.resize(consts::maxNSD,msh.nNo);
          break;

          default:
            throw std::runtime_error("Undefined output");
          break;

        } // switch 

      } // iOut for loop 

    } // iEq for loop 

    if (lIbl) {
      int is = outS[cOut];
      int ie = is;
      outS[cOut+1] = ie + 1;
      outNames[cOut] = "IBLANK";
      cOut = cOut + 1;

      for (int a = 0; a < msh.nNo; a++) {
        int Ac = msh.gN(a);
        for (int i = 0; i < nsd; i++) {
          d[iM].x(i+is,a) = static_cast<double>(com_mod.iblank(Ac));
        }
      } 
    } 

  } // iM for loop 


  // Integrate data from all processors
  //
  int nNo = 0;
  int nEl = 0;

  for (int iM = 0; iM < nMsh; iM++) {
    if (meshes[iM].eType == ElementType::NRB) {
      // CALL INTNRBDATA(msh(iM), d(iM), outDof)
    } else { 
      int_msh_data(com_mod, cm_mod, com_mod.msh[iM], d[iM], outDof, nOute);
    }

    nNo = nNo + d[iM].nNo;
    nEl = nEl + d[iM].nEl;
  }

  if (cm.slv(cm_mod)) {
    com_mod.savedOnce = true;
    return; 
  }

  Array<double> tmpV(consts::maxNSD, nNo);

  // Writing to vtu file (master only)
  //
  std::string fName;

  if (com_mod.cTS > 1000 || lAve) {
    fName = std::to_string(com_mod.cTS);
  } else { 
    std::ostringstream ss;
    ss << std::setw(3) << std::setfill('0') << com_mod.cTS;
    fName = ss.str();
  }

  fName = com_mod.saveName + "_" + fName + ".vtu";
  auto vtk_writer = VtkData::create_writer(fName);

  // Writing the position data
  //
  int iOut = 0;
  int s = outS[iOut];
  int e = outS[iOut+1] - 1;
  int nSh  = 0;
  tmpV = 0.0;

  for (int iM = 0; iM < nMsh; iM++) {
    for (int a = 0; a < d[iM].nNo; a++) {
      for (int i = 0; i < nsd; i++) {
        tmpV(i,a+nSh) = d[iM].gx(i+s,a);
      }
    }
    nSh = nSh + d[iM].nNo;
  }

  vtk_writer->set_points(tmpV);

  // Writing the connectivity data
  //
  nSh = 0;     // For 0-based IDs.
  int num_elems = 0;

  for (int iM = 0; iM < nMsh; iM++) {
    int eNoN = d[iM].eNoN;
    int nEl = d[iM].nEl;
    Array<int> tmpI(eNoN, nEl);

    for (int e = 0; e < nEl; e++) {
      for (int i = 0; i < eNoN; i++) {
        tmpI(i,e) = d[iM].IEN(i,e) + nSh;
      }
    }

    vtk_writer->set_connectivity(nsd, tmpI);
    nSh = nSh + d[iM].nNo;
  }

  // Writing all solutions
  //
  for (int iOut = 1; iOut < nOut; iOut++) {
    int s = outS[iOut];
    int e = outS[iOut+1] - 1;
    int l = e - s + 1;

    Array<double> tmpV(l, nNo);
    int nSh = 0;

    for (int iM = 0; iM < nMsh; iM++) {
      for (int a = 0; a < d[iM].nNo; a++) {
        for (int i = 0; i < l; i++) {
          tmpV(i,a+nSh) = d[iM].gx(i+s,a);
        }
      }

      nSh = nSh + d[iM].nNo;
    }

    vtk_writer->set_point_data(outNames[iOut], tmpV);
  }

  // Write element-based variables
  //
  int ne = -1;

  if (!com_mod.savedOnce || nMsh > 1) {
    Array<int> tmpI(1,nEl);

    // Write the domain ID
    //
    ne = 0;

    if (com_mod.dmnId.size() != 0) {
      int Ec = 0;

      for (int iM = 0; iM < nMsh; iM++) {
        for (int e = 0; e < d[iM].nEl; e++) {
          tmpI(0,Ec) = static_cast<int>(d[iM].xe(e,ne));
          Ec = Ec + 1;
        }
      }
      vtk_writer->set_element_data("Domain_ID", tmpI);
    }

    if (!com_mod.savedOnce) {
      com_mod.savedOnce = true;
      ne = ne + 1;

      // Write partition data
      if (!cm.seq()) {
        int Ec = 0;
        for (int iM = 0; iM < nMsh; iM++) {
          for (int e = 0; e < d[iM].nEl; e++) {
            tmpI(0,Ec) = static_cast<int>(d[iM].xe(e,ne));
            Ec = Ec + 1;
          }
        }
        vtk_writer->set_element_data("Proc_ID", tmpI);
      }
    }

    // Write the mesh ID
    //
    if (nMsh > 1) {
      int Ec = 0;
      for (int iM = 0; iM < nMsh; iM++) {
        for (int e = 0; e < d[iM].nEl; e++) {
          tmpI(0,Ec) = iM;
          Ec = Ec + 1;
        }
      }
      vtk_writer->set_element_data("Mesh_ID", tmpI);
    }
  }  // if (com_mod.savedOnce || nMsh > 1)

  // Write element Jacobian and von Mises stress if necessary
  //
  for (int l = 0; l < nOute; l++) {
    ne = ne + 1;
    Array<double> tmpVe(1,nEl);
    int Ec = 0;

    for (int iM = 0; iM < nMsh; iM++) {
      if (d[iM].xe.size() != 0) {
        for (int e = 0; e < d[iM].nEl; e++) {
          tmpVe(0,Ec) = d[iM].xe(e,ne);
          Ec = Ec + 1;
        }
      }
    }
    vtk_writer->set_element_data(outNamesE[l], tmpVe);
  }

  // Write element ghost cells if necessary
  if (lIbl) {
    ne = ne + 1;
    Array<int> tmpI(1,nEl);
    int Ec = 0;
    for (int iM = 0; iM < nMsh; iM++) {
      if (d[iM].xe.size() != 0) {
        for (int e = 0; e < d[iM].nEl; e++) {
          tmpI(0,Ec) = static_cast<int>(d[iM].xe(e,ne));
          Ec = Ec + 1;
         }
       }
     }
     vtk_writer->set_element_data("EGHOST", tmpI);
  }

  vtk_writer->write();
  delete vtk_writer;
}

};

