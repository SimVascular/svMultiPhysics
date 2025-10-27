// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "LinearAlgebra.h"
#include "PetscLinearAlgebra.h"
#include "FsilsLinearAlgebra.h"
#include "TrilinosLinearAlgebra.h"

const std::map<std::string, consts::LinearAlgebraType> LinearAlgebra::name_to_type = {
  {"none", consts::LinearAlgebraType::none},
  {"fsils", consts::LinearAlgebraType::fsils},
  {"petsc", consts::LinearAlgebraType::petsc},
  {"trilinos", consts::LinearAlgebraType::trilinos}
};

const std::map<consts::LinearAlgebraType, std::string> LinearAlgebra::type_to_name = {
  {consts::LinearAlgebraType::none, "none"},
  {consts::LinearAlgebraType::fsils, "fsils"},
  {consts::LinearAlgebraType::petsc, "petsc"},
  {consts::LinearAlgebraType::trilinos, "trilinos"}
};

/// @brief Check that equation physics is compatible with LinearAlgebra type.
//
void LinearAlgebra::check_equation_compatibility(const consts::EquationType eq_physics,
    const consts::LinearAlgebraType lin_alg_type, const consts::LinearAlgebraType assembly_type)
{
  using namespace consts;

  // ustruct physics requires fsils assembly. 
  //
  if (eq_physics == EquationType::phys_ustruct) {
    if ((lin_alg_type == LinearAlgebraType::trilinos) &&
        (assembly_type != LinearAlgebraType::fsils)) {
      throw std::runtime_error("[svMultiPhysics] Equations with ustruct physics must use fsils for assembly.");
    }
  }
}

LinearAlgebra::LinearAlgebra()
{
}

/// @brief Create objects derived from LinearAlgebra. 
LinearAlgebra* LinearAlgebraFactory::create_interface(consts::LinearAlgebraType interface_type)
{
  LinearAlgebra* interface = nullptr;

  switch (interface_type) {
    case consts::LinearAlgebraType::fsils:
      interface = new FsilsLinearAlgebra();
    break;

    case consts::LinearAlgebraType::petsc:
      interface = new PetscLinearAlgebra();
    break;

    case consts::LinearAlgebraType::trilinos:
      interface = new TrilinosLinearAlgebra();
    break;
  }

  return interface;
}

