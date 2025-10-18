// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "ComMod.h"

#include <iostream>

//--------
// ComMod 
//--------
//
ComMod::ComMod() 
{
  mvMsh        = false;

  stFileFlag   = false;
  stFileRepl   = false;

  bin2VTK      = false;
  saveAve      = false;
  sepOutput    = false;
  saveATS      = 1;
  saveIncr     = 10;
  nITs         = 0;
  startTS      = 0;
  stFileName   = "stFile";
  iniFilePath  = "";
  stopTrigName = "STOP_SIM";
  ichckIEN     = true; 
  zeroAve      = false;
  cmmInit      = false;
  cmmVarWall   = false;
  shlEq        = false;
  pstEq        = false;
  sstEq        = false;
  ibFlag       = false;
  risFlag      = false;

}

//---------
// ~ComMod 
//---------
//
ComMod::~ComMod() 
{
}

///////////////////////
//   a d j T y p e   //
///////////////////////

void adjType::destroy()
{
  nnz = 0;

  pcol.clear();

  prow.clear();
}

///////////////////////
// c p l B C T y p e //
///////////////////////

cplBCType::cplBCType()
{
  schm = consts::CplBCType::cplBC_NA;
}

///////////////////
// d m n T y p e //
///////////////////

dmnType::dmnType()
{
}

dmnType::~dmnType()
{
}

/////////////////
// e q T y p e //
/////////////////

eqType::eqType()
{
  roInf = 0.2;
}

eqType::~eqType()
{
}

////////////////////////////
// s t m o d e l  T y p e //
////////////////////////////

stModelType::stModelType()
{
}

//////////////////////
// f a c e  T y p e //
//////////////////////

faceType::faceType()
{
}

faceType::~faceType()
{
}

/// @brief Free memory and reset some data members.
///
/// This replicates the Fortran 'SUBROUTINE DESTROYFACE(lFa)' 
/// implemented in ALLFUN.f. 
void faceType::destroy()
{
  gE.clear();     
  gN.clear();    
  lN.clear();   
  IEN.clear();  
  gebc.clear();
  w.clear();  
  x.clear(); 
  xi.clear(); 
  N.clear();      
  nV.clear();    
  Nx.clear();   
  Nxx.clear(); 

  nAdj.destroy();
  eAdj.destroy();

  for (int i = 0; i < nFs; i++) {
    fs[i].destroy();
  } 

  eType = consts::ElementType::NA;
  nEl = 0;
  nNo = 0;
  gnEl= 0;
}

///////////////////
// f s T y p e   //
///////////////////

fsType::fsType()
{
}

/// @brief SUBROUTINE DESTROYFS(fs)
void fsType::destroy()
{
  eType = consts::ElementType::NA;

  w.clear();
  xi.clear();
  xib.clear();
  N.clear();
  Nb.clear();
  Nx.clear();
  Nxx.clear();
}

///////////////////
// m s h T y p e //
///////////////////

mshType::mshType()
{
  //std::cout << "+ + + + +  mshType ctor + + + + + " << std::endl;
  eType = consts::ElementType::NA;
}

/////////////////////
// r m s h T y p e //
/////////////////////

rmshType::rmshType()
{
  isReqd  = false;
}


/////////////////////////////////////////////////////////
// s v Z e r o D S o l v e r I n t e r f a c e T y p e //
/////////////////////////////////////////////////////////

void svZeroDSolverInterfaceType::set_data(const svZeroDSolverInterfaceParameters& params)
{
  if (!params.defined()) {
    return;
  }

  configuration_file= params.configuration_file();

  coupling_type = params.coupling_type();

  solver_library = params.shared_library();

  if (params.initial_flows.defined()) { 
    have_initial_flows = true;
    initial_flows = params.initial_flows();
  }

  if (params.initial_pressures.defined()) { 
    have_initial_pressures = true;
    initial_pressures = params.initial_pressures();
  }

  has_data = true;
}

//----------------
// add_block_face
//----------------
// Add a block name / face name pair representing the coupling of a 
// 0D block with a 3D face.
//
void svZeroDSolverInterfaceType::add_block_face(const std::string& block_name, const std::string& face_name)
{
  block_surface_map[block_name] = face_name;
}


