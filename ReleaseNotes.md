# svMultiPhysics 1.0.0.3

This release contains changes made since the June 11, 2026 svMultiPhysics package release.

## Highlights

- Correct the tangent used by 3D-0D coupled structural simulations, substantially improving nonlinear convergence for deforming coupling surfaces. (#802)
- Add Darcy flow physics. (#620)
- Add 3D-1D coupling. (#563)
- Add the Land-Niederer and Regazzoni active-stress models. (#652, #596)
- Improve electromechanics active-stress and force-length coupling. (#650)
- Support multiple and geometrically filtered cellular-electrophysiology stimuli. (#586, #579)
- Add scaffold-surface support and velocity handling for URIS valves. (#569, #565)
- Add experimental line search to the nonlinear integrator. (#701)
- Add tools for initializing coupled 0D RCR boundary conditions. (#613, #637)
- Add validated, package-based container images for the default and Trilinos configurations. (#788)

## Coupled Simulations

- Add the missing shape-derivative contribution to the tangent for 3D-0D coupled structural simulations. (#802)
- Apply the tangent correction to both physical walls and coupling caps. (#802)
- Evaluate the unit pressure load at the intermediate configuration, matching its treatment in the residual. (#802)
- Enforce outward cap-normal orientation when loading cap surfaces. (#802)
- Cache cap connectivity and node mappings instead of repeatedly reconstructing and communicating them. (#802)
- Refactor coupled-boundary-condition support across the FSILS, Trilinos, and PETSc backends. (#802)
- Allow the 0D solver to couple to equations other than the first equation. (#605)
- Add utilities for initializing coupled 0D RCR boundary conditions. (#613, #637)

## Physics and Modeling

- Initialize the 3D pressure field from a scalar value. (#598)
- Add prestress contributions to stress-tensor output. (#675)
- Refactor fiber handling and active-stress implementations. (#725, #578)
- Skip fiber stretch and rate calculations when they are not required by the selected active-stress model. (#610)
- Correct two-dimensional structural simulations. (#562)
- Correct concentration output and its reference tests. (#642)
- Correct the timing of Neumann boundary-condition evaluation to use \(t_{n+\alpha}\). (#777)

## Linear Solvers and Numerical Methods

- Use modified Gram-Schmidt orthogonalization in GMRES. (#600)
- Correct GMRES handling when the initial guess is already converged. (#574)
- Correct GMRES early-termination results. (#630)
- Correct FSILS RCS diagonal cancellation and general block scaling. (#647)
- Prevent the Navier-Stokes solver from breaking down on dependent search directions. (#674)
- Set the MueLu number of equations from the degrees of freedom. (#617)
- Correct weak Dirichlet-condition assembly indices for three-dimensional fluid simulations. (#585)
- Correct the ParMETIS coarse-graph hash reset. (#691)
- Replace variable-length arrays in the 0D coupling implementation. (#713)
- Move repeated matrix multiplication outside loops and use fixed-size operations where appropriate. (#608)

## Output and Reliability

- Make `Vector`, `Array`, and `Array3` behavior consistent when resized to zero or copied from empty objects. (#792)
- Improve VTU output. (#631)
- Correct history output. (#624)
- Preserve spaces in string parameters and stop when `genBC` fails. (#703)
- Correct mesh-data copying. (#732)
- Remove the final restart file only on the master MPI rank. (#721)
- Correct CANN Eigen matrix resizing in Debug builds. (#672)
- Make material convergence tests independent of execution order and random samples. (#699)
- Correct FitzHugh-Nagumo Jacobian handling. (#700)
- Add standalone reference tests for cellular ODE models. (#627)

## Finite-Element Infrastructure

- Introduce the object-oriented finite-element basis implementation. (#561)
- Add bounded Gaussian line quadrature. (#645)
- Add centroid, Duffy, and symmetric triangle quadrature rules. (#735, #745)
- Introduce a common solver math layer. (#727)
- Modularize the solver Core and finite-element components with per-directory CMake targets. (#622)
- Reorganize physics sources and public headers. (#670)
- Use fixed-size Eigen operations for `struct` and `ustruct` calculations and reuse the viscous response for linear elements. (#640)
- Remove unnecessary locking from the solver basis cache. (#686)

## Build and Platform Support

- Add Intel OneAPI compiler and Intel MPI compatibility fixes. (#584)
- Add CMake 4 compatibility. (#698)
- Rebuild solver sources correctly from top-level builds. (#697)
- Improve linking of threads, BLAS/LAPACK, and VTK. (#712)
- Remove obsolete solver functions, includes, and compiler flags. (#711, #716)
- Correct PETSc CMake option syntax and library discovery. (#687, #714)
- Validate builds and tests on Ubuntu 26.04. (#742)

## Packaging and Release Infrastructure

- Embed the release version, source tag, source commit, and solver configuration in the executable. (#781)
- Add reproducible Ubuntu 26 dependency images. (#774)
- Build and validate separate Ubuntu 26 packages for the default and Trilinos configurations. (#783)
- Add macOS Tahoe ARM package paths, identifiers, filenames, and metadata for the default and Trilinos configurations. (#800)
- Build final svMultiPhysics container images from the validated Ubuntu packages. (#788)
- Install the pinned svZeroDSolver release in each container image. (#788)
- Record dependency-image digests, installed package checksums, source revisions, and OCI provenance metadata. (#788)
- Validate each image using `svmultiphysics --version`, shared-library checks, MPI checks, and the integration test suite. (#788)
- Publish the default and Trilinos images only after both configurations pass validation. (#788)
- Protect image publication with the GitHub `release` environment. (#788)
- Make image publication safe to rerun when identical images already exist. (#788)
- Record package checksums, runtime dependencies, container digests, and release provenance.

The executable reports its release identity with:

```console
svmultiphysics --version
```

## Packages and Container Images

The release provides separate Release builds for the default and Trilinos configurations.

The native packages follow the same dependency model as the June 2026 packages: they install
the `svmultiphysics` executable, while platform runtime prerequisites are installed separately.

### macOS Tahoe ARM installation

Install [Homebrew](https://brew.sh/) and the common runtime prerequisites:

```console
brew install open-mpi vtk
```

The Trilinos configuration also requires:

```console
brew install hypre lapack openblas
```

Install one svMultiPhysics configuration:

```console
sudo installer -pkg svMultiPhysics-MacOS-Tahoe-arm-default-1.0.0.3.pkg -target /
```

or:

```console
sudo installer -pkg svMultiPhysics-MacOS-Tahoe-arm-trilinos-1.0.0.3.pkg -target /
```

Verify the installation:

```console
/usr/local/sv/svMultiPhysics/1.0.0.3/bin/svmultiphysics --version
```

Simulations using 3D-0D coupling additionally require the
[svZeroDSolver `v3.2` package](https://simtk.org/projects/simvascular):

```console
sudo installer -pkg svZeroDSolver-MacOS-Tahoe-arm-v3.2.pkg -target /
```

The package installs the coupling library at:

```text
/usr/local/sv/svZeroDSolver/v3.2/lib/libsvzero_interface.dylib
```

Set the solver input file's `Shared_library` value to this installed path.

The locally built default and Trilinos svMultiPhysics packages and the final signed and
notarized svZeroDSolver package were validated on a clean macOS Tahoe ARM virtual machine.
Validation included a standalone svZeroDSolver case, a default 3D fluid case, a case using the
Trilinos backend, and coupled 3D-0D cases with both svMultiPhysics configurations using 1, 3,
and 4 MPI ranks. The official svMultiPhysics packages will be retested before publication.

### Ubuntu 26.04 installation

Install one svMultiPhysics package with `apt`. The package declares its LAPACK and Open MPI
runtime prerequisites, and `apt` installs them automatically. No separate prerequisite
installation command is required when using `apt`:

```console
sudo apt-get update
sudo apt-get install ./svMultiPhysics-Ubuntu-26-default-1.0.0.3.deb
```

or:

```console
sudo apt-get update
sudo apt-get install ./svMultiPhysics-Ubuntu-26-trilinos-1.0.0.3.deb
```

For simulations using 3D-0D coupling, install the pinned
[svZeroDSolver `v3.2` package](https://simtk.org/projects/simvascular)
alongside svMultiPhysics:

```console
sudo apt-get install \
  ./svZeroDSolver-Ubuntu-26-v3.2.deb \
  ./svMultiPhysics-Ubuntu-26-default-1.0.0.3.deb
```

Verify the installation:

```console
/usr/local/sv/svMultiPhysics/1.0.0.3/bin/svmultiphysics --version
```

### macOS Tahoe ARM

- `svMultiPhysics-MacOS-Tahoe-arm-default-1.0.0.3.pkg`
- `svMultiPhysics-MacOS-Tahoe-arm-trilinos-1.0.0.3.pkg`
- `svZeroDSolver-MacOS-Tahoe-arm-v3.2.pkg`

### Ubuntu 26.04

- `svMultiPhysics-Ubuntu-26-default-1.0.0.3.deb`
- `svMultiPhysics-Ubuntu-26-trilinos-1.0.0.3.deb`
- `svZeroDSolver-Ubuntu-26-v3.2.deb`

### Container images

- [`simvascular/svmultiphysics:1.0.0.3-default`](https://hub.docker.com/r/simvascular/svmultiphysics)
- [`simvascular/svmultiphysics:1.0.0.3-trilinos`](https://hub.docker.com/r/simvascular/svmultiphysics)

Checksums, container digests, and release provenance will be published with the artifacts.

## Continuous Integration

- Separate fast required pull-request tests from full main-branch validation. (#728)
- Reduce integration-test execution time. (#719)
- Avoid duplicate pull-request test runs. (#689)
- Skip full tests for documentation-only changes. (#723)
- Add review policies for the object-oriented refactoring effort. (#673, #748, #760)

## Notes for Developers

This release includes substantial internal reorganization of finite-element, physics, linear-solver, coupled-boundary-condition, and common math code. These changes establish clearer component boundaries while preserving the public solver workflow.

The final source commit, package checksums, and container digests will be added after the release artifacts have passed validation.
