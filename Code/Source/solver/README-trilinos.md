# Trilinos Linear Solver Implementation

Developer guide for the Trilinos-based parallel linear solver used in svMultiPhysics.

## Overview

This implementation uses the **Trilinos** framework (specifically **Tpetra**, **Belos**, **Ifpack2**, and **MueLu** packages) to solve large-scale sparse linear systems arising from finite element discretizations in parallel MPI environments. For detailed instructions on the Trilinos built used for this project, please refer to the svMultiPhysics/Docker/ubuntu/dockerfile. The instructions found in the dockerfile can be easily adapted for a Linux--based installation.

## User guide
The Trilinos interface of svMultiPhysics allows users to select linear solvers and preconditioners from the Trilinos linear algebra package for solving large-scale sparse linear systems. This section describes the solvers and preconditioners available through the interface, using the exact names that can be specified in the input .xml file. The linear solvers available are:
- `GMRES`: Generalized Minimal Residual method, a robust and widely used Krylov--subspace solver suitable for general non-symmetric sparse systems. It is typically used for CFD and FSI applications due to its stability and flexibility.
- `CG`: Conjugate Gradient method, designed for symmetric positive-definite systems. Within svMultiPhysics, this method should be used for solving the mesh motion equation in the FSI ALE framework, where the system satisfies these mathematical properties.
- `BiCGS`: Bi--Conjugate Gradient Stabilized, an iterative solver for general non--symmetric problems. It is often used when GMRES becomes too expensive (e.g., due to memory growth with restart size).
The preconditioners available are:
- `trilinos-diagonal`: the simplest option, performing a row and column scaling using the inverse square root of the diagonal entries. This preconditioner is inexpensive and can be useful for mildly ill-conditioned systems.
- `trilinos-blockjacobi`: applies block-Jacobi scaling using the diagonal block submatrices. This often provides better conditioning than simple diagonal scaling, especially for block-structured systems.
- `trilinos-ilu`: computes an Incomplete LU (ILU) factorization within an additive Schwarz framework. The ILU has no additional fill beyond the original sparsity pattern, making it relatively efficient while still improving solver robustness.
- `trilinos-ilut`: similar to ILU but with a controlled threshold (set to 1.0e-2) that determines the amount of fill. A higher threshold generates a more accurate (but more expensive) preconditioner. Implemented using the additive Schwarz method.
- `trilinos-riluk0`: a reduced ILU(0) factorization with no fill. Implemented with additive Schwarz. It is fast to compute and is typically an excellent choice for FSI problems using velocity--based structural formulations (`ustruct` equation for solid).
- `trilinos-riluk1`: similar to riluk0 but includes one level of fill, improving accuracy at the cost of additional computation.
- `trilinos-ml`: an algebraic multigrid (AMG) preconditioner from the MueLu package. This is generally the preferred choice for displacement--based structural problems in FSI (`struct` equation), and it has also proven effective for large--scale CFD systems. AMG preconditioners can be expensive to build but scale very well for large problems, making them ideal for high-resolution simulations. 
### Notes on the trilinos-ml
AMG can be extremely powerful, but only when configured appropriately for the specific PDE, mesh size, and physical model. Because of the large number of tunable parameters and because optimal settings differ between CFD, FSI, and structural mechanics, exposing everything through the standard .xml input file would be cumbersome and hard--coding problem specific AMG parameters in the source code requires solver rebuilding everytime a change is implemented.
Instead, svMultiPhysics allows a dedicated AMG parameters file, letting advanced users:
- override defaults (hard--coded in the solver already) for smoother and coarsening choices
- tailor AMG behavior to their problem class
- experiment with performance tuning without modifying the main simulation input files
This approach keeps the main input format clean while still providing expert--level control for those who need it. In order to use the AMG parameters from file, a file named exactly `mueluOptions.xml` must be in the same folder where the simulation input file is. An example `mueluOptions.xml` file is included in the svMultiPhysics/Code/Source/solver directory.

**Files:**
- `trilinos_impl.h` — type definitions, function declarations, and key data structures
- `trilinos_impl.cpp` — implementation of assembly, matrix construction, and solve routines

---

## Key Data Structures

### 1. `Trilinos` Struct

Central container holding all Tpetra/Trilinos objects for a linear system solve:

```cpp
struct Trilinos {
  Teuchos::RCP<const Tpetra_Map> Map;           // DOF ownership map (owned nodes only)
  Teuchos::RCP<const Tpetra_Map> ghostMap;      // DOF map including ghost nodes
  Teuchos::RCP<Tpetra_MultiVector> F;           // RHS vector (owned DOFs)
  Teuchos::RCP<Tpetra_MultiVector> ghostF;      // RHS vector (owned + ghost DOFs)
  Teuchos::RCP<Tpetra_CrsMatrix> K;             // Global stiffness matrix
  Teuchos::RCP<Tpetra_Vector> X;                // Solution vector (owned DOFs)
  Teuchos::RCP<Tpetra_Vector> ghostX;           // Solution vector (owned + ghost)
  Teuchos::RCP<Tpetra_Import> Importer;         // Import object (ghost communication)
  std::vector<Teuchos::RCP<Tpetra_MultiVector>> bdryVec_list;  // Coupled Neumann BCs
  Teuchos::RCP<const Teuchos::Comm<int>> comm;  // MPI communicator
  Teuchos::RCP<Tpetra_CrsGraph> K_graph;        // Sparse graph structure

  Teuchos::RCP<Tpetra_Operator> MueluPrec;      // MueLu (algebraic multigrid) preconditioner
  Teuchos::RCP<Ifpack2_Preconditioner> ifpackPrec;  // Ifpack2 preconditioner
};
```

**Ownership model:**
- `Map`: defines distribution of DOFs across MPI ranks (owned DOFs only, non-overlapping).
- `ghostMap`: extends `Map` to include ghost DOFs (nodes shared with neighboring ranks).
- Vectors (`F`, `X`) use `Map`; ghost vectors (`ghostF`, `ghostX`) use `ghostMap`.

---

### 2. Tpetra Type Aliases

Defined in `trilinos_impl.h` for clarity and portability:

| Alias | Trilinos Type | Description |
|-------|---------------|-------------|
| `Tpetra_Map` | `Tpetra::Map<LO, GO, Node>` | Parallel distribution map |
| `Tpetra_CrsMatrix` | `Tpetra::CrsMatrix<Scalar_d, LO, GO, Node>` | Sparse matrix (CSR format) |
| `Tpetra_CrsGraph` | `Tpetra::CrsGraph<LO, GO, Node>` | Sparse graph (topology) |
| `Tpetra_MultiVector` | `Tpetra::MultiVector<Scalar_d, LO, GO, Node>` | Multi-vector (RHS, etc.) |
| `Tpetra_Vector` | `Tpetra::Vector<Scalar_d, LO, GO, Node>` | Single vector (solution) |
| `Tpetra_Operator` | `Tpetra::Operator<Scalar_d, LO, GO, Node>` | Abstract linear operator |
| `Belos_LinearProblem` | `Belos::LinearProblem<...>` | Belos problem wrapper |
| `Belos_SolverFactory` | `Belos::TpetraSolverFactory<...>` | Creates Belos solvers |

**Ordinals:**
- `LO` (local ordinal): `int` — indices within a rank
- `GO` (global ordinal): `int` — global node/DOF indices across all ranks
- `Node`: Kokkos node type (typically default)

---

### 3. `TrilinosMatVec` Operator

Custom `Tpetra::Operator` that applies the linear operator:

$$
A x = K x + \sum_i v_i (v_i^T x)
$$

where:
- $K$ is the global stiffness matrix
- $v_i$ are coupled Neumann boundary vectors (rank-1 updates for resistance BCs)

**Purpose:** enables efficient matrix-free application of boundary terms without explicitly forming outer products $v_i v_i^T$.

**Key method:**
```cpp
void TrilinosMatVec::apply(const Tpetra_MultiVector& x, Tpetra_MultiVector& y, ...)
```
- Computes `y = K*x` via `trilinos_->K->apply(...)`
- Adds coupled Neumann contributions: `y += v_i * (v_i^T * x)` for each boundary vector

---

## Workflow

### Phase 1: Initialization and Graph Construction

**Function:** `trilinos_lhs_create(...)`

**Inputs:**
- `numGlobalNodes`: total mesh nodes (across all ranks)
- `numLocalNodes`: nodes owned by this rank
- `numGhostAndLocalNodes`: owned + ghost nodes
- `nnz`: number of nonzero entries in CSR node-adjacency structure
- `ltgSorted`: local-to-global map (sorted/reordered for better partitioning)
- `ltgUnsorted`: local-to-global map (original CSR ordering)
- `rowPtr`, `colInd`: CSR row pointers and column indices (node-level connectivity)
- `Dof`: degrees of freedom per node (e.g., 3 for 3D velocity, 4 for velocity+pressure)

**Steps:**

1. **Create DOF-based maps:**
   - Expand node GIDs to DOF GIDs: `dofGID = nodeGID * dof + d`
   - Build `Map` (owned DOFs) and `ghostMap` (owned + ghost DOFs)

2. **Build sparse graph (`K_graph`):**
   - Compute `nnzPerDofRow`: for each DOF row (using Map ordering), determine number of nonzeros by mapping node GID → unsorted index → `rowPtr` to get node neighbor count, then multiply by `dof`.
   - Construct `Tpetra_CrsGraph` with pre-allocated `nnzPerDofRow`.
   - Insert global column indices via `insertGlobalIndices(rowGID, rowCols)` for each DOF row.
   - Call `fillComplete()` to finalize graph (communication, optimization).

3. **Create matrix and vectors:**
   - Matrix `K` built from finalized graph.
   - Vectors `F`, `ghostF`, `X`, `ghostX` created from respective maps.
   - `Importer` created for ghost node communication (`Map` → `ghostMap`).

**Key detail:** DOF expansion — node-level CSR connectivity is expanded to DOF-level by iterating over `dof × dof` blocks per node-node connection.

---

### Phase 2: Assembly

Two assembly modes are supported:

#### A. Element-by-Element Assembly (Trilinos-native)

**Function:** `trilinos_doassem_(...)`

**Inputs:**
- `numNodesPerElement`: nodes in current finite element
- `eqN`: element connectivity (local node indices → proc-local indices)
- `lK`: element stiffness matrix (dense, size `dof*dof × numNodesPerElement × numNodesPerElement`)
- `lR`: element force vector (size `dof × numNodesPerElement`)

**Steps:**
1. Map element nodes to global node IDs via `localToGlobalUnsorted`.
2. For each element node pair `(a, b)`:
   - Compute global row/col DOF indices: `rowGID = nodeGID_a * dof + i`, `colGID = nodeGID_b * dof + j`.
   - Extract block `lK[a,b]` (size `dof × dof`).
   - Call `K->sumIntoGlobalValues(rowGID, dof, vals, cols)` to accumulate into global matrix.
3. Sum force vector entries into `ghostF` via `sumIntoGlobalValue(...)`.

**Note:** assembly uses `ghostF` (includes ghost nodes); final communication done later.

#### B. Global Assembly (FSILS-assembled)

**Function:** `trilinos_global_solve_(...)`

**Inputs:**
- `Val`: preassembled matrix values (CSR format, already assembled by FSILS)
- `RHS`: preassembled RHS vector

**Steps:**
1. Loop over all nodes (owned + ghost), extract CSR rows from `Val` and `RHS`.
2. Insert/sum values into `K` and `ghostF` using global indices.
3. Proceed to solve (next phase).

**Use case:** when matrix/RHS are assembled externally (e.g., legacy FSILS assembly), this path avoids redundant element loops.

---

### Phase 3: Solve

**Function:** `trilinos_solve_(...)`

**Inputs:**
- `dirW`: Dirichlet boundary condition weights (1 for free DOFs, 0 for constrained)
- Solver parameters: `lsType` (GMRES, BiCGStab, CG), `relTol`, `maxIters`, `kspace` (Krylov space size)
- `precondType`: preconditioner selection (diagonal, block Jacobi, ILU, ILUT, RILUk, ML)

**Steps:**

1. **Finalize matrix:**
   - `K->fillComplete()` — finalizes parallel assembly (ghost communication, CRS optimization).

2. **Export RHS:**
   - `F->doExport(*ghostF, exporter, Tpetra::ADD)` — sum contributions from ghost nodes to owned nodes (if using element assembly).
   - Or `REPLACE` mode if RHS already assembled correctly.

3. **Construct Jacobi scaling (diagonal preconditioning + Dirichlet BCs):**
   - Extract diagonal of `K`.
   - Modify diagonal: set zero entries to 1, compute `D^{-1/2}`.
   - Apply Dirichlet weights from `dirW`.
   - Left-scale and right-scale `K` by diagonal: `K ← D^{-1/2} K D^{-1/2}`.
   - Scale `F` and boundary vectors similarly.

4. **Setup Belos linear problem:**
   - Operator: `TrilinosMatVec` (applies `K` + coupled Neumann terms due to resistance or RCR BCs).
   - Solution: `X`, RHS: `F`.
   - Preconditioner: set via `setPreconditioner(...)` (creates Ifpack2 or MueLu preconditioner).

5. **Create and configure Belos solver:**
   - Solver type: Block GMRES, BiCGStab, or Pseudoblock CG.
   - Parameters: convergence tolerance, max iterations, Krylov space dimension, verbosity.
   - Factory pattern: `Belos::TpetraSolverFactory` creates solver manager.

6. **Solve:**
   - `solverManager->solve()` — iterative solve.
   - Returns convergence status, iteration count, residual norms.

7. **Post-process solution:**
   - Unscale `X` by multiplying with diagonal.
   - Import solution to ghost nodes: `ghostX->doImport(*X, *Importer, Tpetra::INSERT)`.
   - Copy `ghostX` to output array `x`.

8. **Cleanup:**
   - Zero out `F`, `ghostF`, `X`, boundary vectors.
   - Nullify preconditioners and matrix.

---

## Preconditioners

Implemented in `setPreconditioner(...)`:

| Preconditioner | Type | Description |
|----------------|------|-------------|
| `NO_PRECONDITIONER` | None | No preconditioning (may converge slowly) |
| `TRILINOS_DIAGONAL_PRECONDITIONER` | Diagonal | Jacobi scaling only (handled separately) |
| `TRILINOS_BLOCK_JACOBI_PRECONDITIONER` | Ifpack2 Jacobi | Block Jacobi relaxation (1 sweep) |
| `TRILINOS_ILU_PRECONDITIONER` | Ifpack2 Schwarz+ILU(0) | Overlapping domain decomposition with ILU(0) |
| `TRILINOS_ILUT_PRECONDITIONER` | Ifpack2 Schwarz+ILUT | ILU with threshold dropping |
| `TRILINOS_RILUK0_PRECONDITIONER` | Ifpack2 RILUK(0) | Relaxed ILU (level 0) |
| `TRILINOS_RILUK1_PRECONDITIONER` | Ifpack2 RILUK(1) | Relaxed ILU (level 1) |
| `TRILINOS_ML_PRECONDITIONER` | MueLu | Algebraic multigrid (smoothed aggregation) |

**MueLu (ML) preconditioner:**
- Configured in `setMueLuPreconditioner(...)`.
- Parameters tuned for large FSI problems: 6 levels, V-cycle, Gauss-Seidel smoother, KLU coarse solver.
- Can load custom parameters from `mueluOptions.xml` if present.

---

## Coupled Neumann Boundary Conditions

**Motivation:** resistance boundary conditions add low-rank updates to the stiffness matrix:

$$
A = K + \sum_i v_i v_i^T
$$

where $v_i$ are boundary vectors (normal vectors scaled by resistance coefficients).

**Implementation:**
- Vectors stored in `trilinos_->bdryVec_list` (one per coupled BC face).
- Populated in `trilinos_bc_create_(...)` from user-provided BC data.
- Applied during matrix-vector multiply via `TrilinosMatVec::apply(...)`:
  - Compute dot product: $\alpha_i = v_i^T x$
  - Add contribution: $y \leftarrow y + \alpha_i v_i$

**Advantages:**
- Avoids explicitly forming outer product matrices (memory-intensive).
- Efficient vectorized operations.

---

## Graph and Matrix Construction Details

### DOF Ordering

- **Node-based CSR:** input adjacency (`rowPtr`, `colInd`) is node-to-node connectivity.
- **DOF expansion:** each node has `dof` degrees of freedom; global DOF index = `nodeGID * dof + d`.
- **Map ordering vs CSR ordering:**
  - `ltgSorted`: node ordering used to build Tpetra `Map` (may differ from CSR order for better partitioning).
  - `ltgUnsorted`: original CSR ordering.
  - Graph construction uses a hash map (`gidToUnsortedIndex`) to map node GID → unsorted index → `rowPtr` for computing neighbor counts.

### Allocation Strategy

- Pre-allocate `nnzPerDofRow` for each DOF row (reduces memory overhead).
- For each DOF `d` of node `n`:
  - Find node's neighbor count: `rowPtr[unsortedIdx+1] - rowPtr[unsortedIdx]`.
  - Allocate `neighborCount * dof` entries (because each node-node connection expands to `dof × dof` block).

---

## Debugging and Diagnostics

### Matrix/Vector Printing

Functions available for debugging (write ASCII files):
- `printMatrixToFile(trilinos_)` → writes `K.txt` (global row/col/value triples).
- `printRHSToFile(trilinos_)` → writes `F.txt`.
- `printSolutionToFile(trilinos_)` → writes `X.txt`.

---

## Solver Types

Defined in `trilinos_impl.h`:

| Macro | Belos Solver | Use Case |
|-------|--------------|----------|
| `TRILINOS_GMRES_SOLVER` | Block GMRES | General nonsymmetric systems (default) |
| `TRILINOS_BICGSTAB_SOLVER` | BiCGStab | Faster for some nonsymmetric problems |
| `TRILINOS_CG_SOLVER` | Pseudoblock CG | Symmetric positive definite systems |

**GMRES parameters:**
- `kspace`: Krylov subspace dimension (typically 50-300).
- `maxRestarts`: computed as `maxIters / kspace + 1`.
- Orthogonalization: DGKS (recommended for robustness).

---

## Memory Management

- **Reference-counted pointers:** all Trilinos objects use `Teuchos::RCP` (similar to `std::shared_ptr`).
- **Cleanup:** at end of solve, set `K`, preconditioners to `Teuchos::null` to free memory.
- **Kokkos:** initialized once (`Kokkos::initialize()`) and finalized at program exit.

---

## Integration Points

### Called from svMultiPhyisics

Functions bridge C++ Trilinos to svMultiPhysics implementation:
- `trilinos_lhs_create(...)` — graph/matrix setup (once per timestep/Newton iteration).
- `trilinos_doassem_(...)` — default element assembly (called per element).
- `trilinos_global_solve_(...)` — solve with pre-assembled data.
- `trilinos_bc_create_(...)` — setup coupled Neumann BCs.

### TrilinosLinearAlgebra Class

Higher-level C++ wrapper (`TrilinosLinearAlgebra::TrilinosImpl`) provides:
- `alloc(...)` — allocate data structures.
- `assemble(...)` — element-by-element assembly.
- `solve(...)` / `solve_assembled(...)` — solve linear system.
- `init_dir_and_coup_neu(...)` — setup boundary conditions.

---

## Tips for new development

1. **Adding a new preconditioner:**
   - Define a new macro in `trilinos_impl.h` (e.g., `#define TRILINOS_CUSTOM_PREC 709`).
   - Add case in `setPreconditioner(...)` to configure Ifpack2/MueLu parameters.
   - If the preconditioner is user--defined, then it should be implemented as a new Tpetra::Operator class and added to the Belos problem.

2. **Changing graph structure:**
   - Graph is fixed after `fillComplete()`.
   - To modify: destroy graph (`K_graph = Teuchos::null`), rebuild in `trilinos_lhs_create(...)`.

3. **Debugging convergence issues:**
   - Enable Belos verbosity: comment out `#define NOOUTPUT` in `trilinos_impl.cpp`.
   - Print matrix/RHS: call `printMatrixToFile()`, `printRHSToFile()`.
   - Check conditioning: ensure diagonal has no zeros (handled by `checkDiagonalIsZero()`).
   - To check the sparse pattern of the matrix, consider using hdf5 library. This will allow to save large files that can be read by using a simple Python script

4. **Performance profiling:**
   - Use Teuchos::Time for timing sections (example: `Belos Solve Timer` in code).
   - Trilinos has built-in profiling via `Teuchos::TimeMonitor` (enable with CMake flag).

---

## References

- [Trilinos Documentation](https://trilinos.github.io/)
- [Tpetra User Guide](https://docs.trilinos.org/dev/packages/tpetra/doc/html/index.html)
- [Belos User Guide](https://docs.trilinos.org/dev/packages/belos/doc/html/index.html)
- [MueLu User Guide](https://trilinos.github.io/pdfs/mueluguide.pdf)
- [Ifpack2 Documentation](https://docs.trilinos.org/dev/packages/ifpack2/doc/html/index.html)

---

**Last updated:** November 2025
