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

///////////////////////////////////////////////////////////
//      S e t   M a t e r i a l   P r o p e r t i e s    //
///////////////////////////////////////////////////////////
//
// The 'set_material_props' map defined here sets material
// properties from values read in from a file.

using SeMaterialPropertiesMapType = std::map<consts::ConstitutiveModelType, 
    std::function<void(DomainParameters*, double, double, double, dmnType&)>>;

SeMaterialPropertiesMapType set_material_props = {

//---------------------------//
//       stIso_nHook         //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_nHook, [](DomainParameters* domain_params, double mu, double kap, double lam, 
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_nHook;
  lDmn.stM.C10 = 0.5 * mu;
} },

//---------------------------//
//       stIso_lin           //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_lin, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_lin;
  lDmn.stM.C10 = mu;
} },

//---------------------------//
//       stIso_StVK          //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_StVK, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_StVK;
  lDmn.stM.C10 = lam;
  lDmn.stM.C01 = mu;
  lDmn.stM.Kpen = kap;
} },

//---------------------------//
//       stIso_mStVK         //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_mStVK, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_mStVK;
  lDmn.stM.C10 = kap;
  lDmn.stM.C01 = mu;
  lDmn.stM.Kpen = kap;
} },

//---------------------------//
//       stIso_MR            //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_MR, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_MR;
  lDmn.stM.C10 = domain_params->constitutive_model.mooney_rivlin.c1.value();
  lDmn.stM.C01 = domain_params->constitutive_model.mooney_rivlin.c2.value();
} },

//---------------------------//
//       stIso_HGO           //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_HGO, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_HGO;
  auto& params = domain_params->constitutive_model.holzapfel_gasser_ogden;

  lDmn.stM.C10 = mu * 0.5;
  lDmn.stM.aff = params.a4.value();
  lDmn.stM.bff = params.b4.value();
  lDmn.stM.ass = params.a6.value();
  lDmn.stM.bss = params.b6.value();
  lDmn.stM.kap = params.kappa.value();
} },

//---------------------------//
//       stIso_Gucci         //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_Gucci, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_Gucci;
  auto& params = domain_params->constitutive_model.guccione;

  lDmn.stM.C10 = params.c.value();
  lDmn.stM.bff = params.bf.value();
  lDmn.stM.bss = params.bt.value();
  lDmn.stM.bfs =params.bfs.value(); 
} },

//---------------------------//
//       stIso_HO            //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_HO, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_HO;
  auto& params = domain_params->constitutive_model.holzapfel;

  lDmn.stM.a = params.a.value();
  lDmn.stM.b = params.b.value();
  lDmn.stM.aff = params.a4f.value(); 
  lDmn.stM.bff = params.b4f.value(); 
  lDmn.stM.ass = params.a4s.value();
  lDmn.stM.bss = params.b4s.value();
  lDmn.stM.afs = params.afs.value();
  lDmn.stM.bfs = params.bfs.value();
  lDmn.stM.khs = params.k.value();
} },

//---------------------------//
//       stIso_HO_ma          //
//---------------------------//
//
{consts::ConstitutiveModelType::stIso_HO_ma, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_HO_ma;
  auto& params = domain_params->constitutive_model.holzapfel;
  lDmn.stM.a = params.a.value();
  lDmn.stM.b = params.b.value();
  lDmn.stM.aff = params.a4f.value(); 
  lDmn.stM.bff = params.b4f.value(); 
  lDmn.stM.ass = params.a4s.value();
  lDmn.stM.bss = params.b4s.value();
  lDmn.stM.afs = params.afs.value();
  lDmn.stM.bfs = params.bfs.value();
  lDmn.stM.khs = params.k.value();
} },

//--------------------------------------------//
//       stArtificialNeuralNet    //
//--------------------------------------------//
//
{consts::ConstitutiveModelType::stArtificialNeuralNet, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stArtificialNeuralNet;
  auto& params = domain_params->constitutive_model.cann;

  lDmn.stM.paramTable.num_rows = params.rows.size();

  lDmn.stM.paramTable.invariant_indices.resize(lDmn.stM.paramTable.num_rows);
  lDmn.stM.paramTable.activation_functions.resize(lDmn.stM.paramTable.num_rows,3);
  lDmn.stM.paramTable.weights.resize(lDmn.stM.paramTable.num_rows,3);

  // Populate components of the table in stM
  for (size_t i = 0; i < lDmn.stM.paramTable.num_rows; i++)
  {
    // Store invariant index
    lDmn.stM.paramTable.invariant_indices[i] = params.rows[i]->row.invariant_index.value_; 

    // Store activation function and weight values
    for (size_t j = 0; j < 3; j++)
    {
      lDmn.stM.paramTable.activation_functions(i,j) = params.rows[i]->row.activation_functions.value_[j];
      lDmn.stM.paramTable.weights(i,j) = params.rows[i]->row.weights.value_[j];
    }

  }
  
} },

//---------------------------//
//       stIso_LS            //
//---------------------------//
// Lee-Sacks material model.
//
{consts::ConstitutiveModelType::stIso_LS, [](DomainParameters* domain_params, double mu, double kap, double lam,
    dmnType& lDmn) -> void
{
  lDmn.stM.isoType = consts::ConstitutiveModelType::stIso_LS;
  auto& params = domain_params->constitutive_model.lee_sacks;

  lDmn.stM.a = params.a.value();
  lDmn.stM.a0 = params.a0.value();
  lDmn.stM.b1 = params.b1.value();
  lDmn.stM.b2 = params.b2.value();
  lDmn.stM.mu0 = params.mu0.value();

} },

};
