// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// This file contains PETSc-dependent data structures. 

#ifndef PETSC_INTERFACE_H
#define PETSC_INTERFACE_H

#include <petscksp.h>
#include <petscao.h>
#include <stdbool.h>

#include <utility>
#include <vector>

#include "consts.h"

struct PetscFaceTangent {
    PetscReal resistance;
    std::vector<std::pair<PetscInt, PetscReal>> owned_rows;
    std::vector<PetscInt> columns;
    std::vector<PetscReal> flux;
};

//--------
// LHSCtx
//---------
// PETSc lhs context. 
//
typedef struct {
    PetscBool created;  /* Whether petsc lhs is created */

    PetscInt  nNo;      /* local number of vertices */
    PetscInt  mynNo;    /* number of owned vertices */

    PetscInt *map;      /* local to local mapping, map[O2] = O1 */
    PetscInt *ltg;      /* local to global in PETSc ordering */
    PetscInt *ghostltg; /* local to global in PETSc ordering */
    PetscInt *rowPtr;   /* row pointer for adjacency info */
    PetscInt *colPtr;   /* column pointer for adjacency info */
} LHSCtx;

//-------
// LSCtx
//-------
// PETSc linear solver context. 
//
typedef struct {
    PetscBool created;  /* Whether mat and vec is created */
    const char *pre;      /* option prefix for different equations */

    PetscInt  DirPts;   /* number of dofs with Dirichlet BC */
    PetscInt *DirBC;    /* PETSc index for dofs with Dirichlet BC */

    Vec       b;        /* rhs/solution vector of owned vertices */
    Mat       A;        /* stiffness matrix */
    KSP       ksp;      /* linear solver context */

    PetscBool rcs;      /* whether rcs preconditioner is activated */
    Vec       Dr;       /* diagonal matrix from row maxabs */
    Vec       Dc;       /* diagonal matrix from col maxabs */
} LSCtx;

void petsc_initialize(const PetscInt nNo, const PetscInt mynNo, 
    const PetscInt nnz, const PetscInt nEq, const PetscInt *svFSI_ltg, 
    const PetscInt *svFSI_map, const PetscInt *svFSI_rowPtr, 
    const PetscInt *svFSI_colPtr, char *inp);

/// @brief Create the vector and the matrix of the linear system of an
/// equation.
///
/// Does nothing if they are already created. The nonzero pattern of the
/// matrix is fixed here. It holds the entries of the mesh adjacency and, for
/// the lumped parameter boundary conditions, one entry per pressure-load row and flow-derivative
/// column on the same coupled face.
///
/// The O2 index of a degree of freedom is the O2 index of its node times dof,
/// plus its component.
///
/// @param[in] dof Number of degrees of freedom per node.
/// @param[in] iEq Index of the equation.
/// @param[in] nEq Number of equations.
/// @param[in] svFSI_DirBC For each degree of freedom on this process, by O2
///   index, zero if a Dirichlet boundary condition constrains it and
///   nonzero otherwise.
/// @param[in] faces Owned rows and global columns for each coupled face.
void petsc_create_linearsystem(const PetscInt dof, const PetscInt iEq,
                               const PetscInt nEq, const PetscReal *svFSI_DirBC,
                               const std::vector<PetscFaceTangent>& faces);

void petsc_create_linearsolver(const consts::PreconditionerType lsType, const consts::PreconditionerType pcType, 
    const PetscInt kSpace, const PetscInt maxIter, const PetscReal relTol, 
    const PetscReal absTol, const consts::EquationType phys, const PetscInt dof, 
    const PetscInt iEq, const PetscInt nEq);

/// @brief Fill the vector and the matrix of the linear system of an equation.
///
/// Each coupled face adds its signed pressure-load/flow-derivative outer
/// product, including columns owned by other processes. The Dirichlet boundary
/// conditions are then applied, and the system is scaled if the row and
/// column scaling preconditioner is active.
///
/// @param[in] dof Number of degrees of freedom per node.
/// @param[in] iEq Index of the equation.
/// @param[in] R Right-hand side, dof values for each node on this process in
///   O1 ordering.
/// @param[in] Val Matrix, dof * dof values for each entry of the mesh
///   adjacency.
/// @param[in] svFSI_DirBC For each degree of freedom on this process, by O2
///   index, zero if a Dirichlet boundary condition constrains it and
///   nonzero otherwise.
/// @param[in] faces Signed tangent and its two vectors for each coupled face.
void petsc_set_values(const PetscInt dof, const PetscInt iEq,
                      const PetscReal *R, const PetscReal *Val,
                      const PetscReal *svFSI_DirBC,
                      const std::vector<PetscFaceTangent>& faces);

void petsc_solve(PetscReal *resNorm,  PetscReal *initNorm,  PetscReal *dB, 
    PetscReal *execTime, bool *converged, PetscInt *numIter, 
    PetscReal *R, const PetscInt maxIter, const PetscInt dof, 
    const PetscInt iEq);

void petsc_destroy_all(const PetscInt);

PetscErrorCode petsc_create_lhs(const PetscInt, const PetscInt, const PetscInt,  
                                const PetscInt *, const PetscInt *, 
                                const PetscInt *, const PetscInt *);

/// @brief Store the Dirichlet degrees of freedom of an equation.
///
/// @param[in] dof Number of degrees of freedom per node.
/// @param[in] cEq Index of the equation.
/// @param[in] svFSI_DirBC For each degree of freedom on this process, by O2
///   index, zero if a Dirichlet boundary condition constrains it and
///   nonzero otherwise.
///
/// @return A PETSc error code.
PetscErrorCode petsc_create_bc(const PetscInt dof, const PetscInt cEq,
                               const PetscReal *svFSI_DirBC);

PetscErrorCode petsc_create_vecmat(const PetscInt, const PetscInt, const PetscInt,
                                  const std::vector<PetscFaceTangent>&);

PetscErrorCode petsc_set_vec(const PetscInt, const PetscInt, const PetscReal *);

PetscErrorCode petsc_set_mat(const PetscInt, const PetscInt, const PetscReal *);

/// @brief Add the lumped parameter boundary conditions to the matrix of an
/// equation, then apply its Dirichlet boundary conditions to the matrix and
/// the right-hand side.
///
/// @param[in] cEq Index of the equation.
/// @param[in] DirBC For each degree of freedom on this process, by O2 index,
///   zero if a Dirichlet boundary condition constrains it and nonzero
///   otherwise.
/// @param[in] faces Per-face tangents, as in petsc_set_values().
///
/// @return A PETSc error code.
PetscErrorCode petsc_set_bc(const PetscInt cEq, const PetscReal *DirBC,
                            const std::vector<PetscFaceTangent>& faces);

PetscErrorCode petsc_set_pcfieldsplit(const PetscInt, const PetscInt);

PetscErrorCode petsc_pc_rcs(const PetscInt, const PetscInt);


PetscErrorCode petsc_debug_save_vec(const char *, Vec);
PetscErrorCode petsc_debug_save_mat(const char *, Mat);

#endif
