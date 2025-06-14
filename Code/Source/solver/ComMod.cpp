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


