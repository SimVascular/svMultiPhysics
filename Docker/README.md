# Docker images

Docker Engine is required to build and run these images. On macOS, Docker Desktop
provides the engine; on Linux, Docker Engine can be installed directly. The
release images currently target `linux/amd64`.

## Image roles

- `Docker/ubuntu/dockerfile` builds the Ubuntu 26.04 dependency image containing
  the libraries and tools required to build and run svMultiPhysics. It does not
  contain svMultiPhysics itself. The `Build Ubuntu 26 dependency image` workflow
  builds, validates, and publishes this image to GitHub Container Registry.
- `Docker/svmultiphysics/dockerfile` builds a final svMultiPhysics runtime image
  from a digest-pinned dependency image. It installs verified svMultiPhysics and
  svZeroDSolver Ubuntu packages rather than compiling source in the image.
- `Docker/solver/dockerfile` is the legacy Docker Hub recipe. It builds PETSc and
  Trilinos configurations from source using `simvascular/libraries:latest` and
  is maintained by the existing `Build and Push Docker Image` workflow.

## Final svMultiPhysics images

Each final release publishes two separate images:

```text
ghcr.io/simvascular/svmultiphysics:<version>-default
ghcr.io/simvascular/svmultiphysics:<version>-trilinos
```

The `default` image contains the default svMultiPhysics configuration. The
`trilinos` image contains the Trilinos configuration. Both images:

- use an Ubuntu 26.04 AMD64 dependency image pinned by digest;
- install the corresponding Release svMultiPhysics package;
- install the svZeroDSolver package pinned by `release.json`;
- expose `svmultiphysics` on `PATH`; and
- record source, package, dependency-image, and workflow provenance in image
  labels.

No floating `latest` tag is published by this workflow. Select an explicit
release version and configuration.

### Pull and inspect an image

Replace `<version>` with a four-part svMultiPhysics release version:

```sh
docker pull ghcr.io/simvascular/svmultiphysics:<version>-default
docker run --rm --platform linux/amd64 \
  ghcr.io/simvascular/svmultiphysics:<version>-default \
  svmultiphysics --version
```

Use the `-trilinos` tag instead to run the Trilinos configuration.

### Run a case

Mount the case directory at `/case`, make it the working directory, and invoke
the installed executable directly:

```sh
docker run --rm --platform linux/amd64 \
  -v /absolute/path/to/case:/case \
  -w /case \
  ghcr.io/simvascular/svmultiphysics:<version>-default \
  mpirun --allow-run-as-root --oversubscribe -n 4 \
    svmultiphysics solver.xml
```

Files written under `/case` remain on the host after the container exits.

## Release workflow

The `Build svMultiPhysics images` workflow accepts a dependency-image reference
and a final svMultiPhysics release tag. It:

1. resolves and verifies the source tag, commit, release metadata, dependency
   image digest, and package checksums;
2. builds the `default` and `trilinos` images in parallel;
3. saves each complete image as a short-lived workflow artifact;
4. loads and validates both artifacts independently, including version output,
   labels, package provenance, shared libraries, MPI, and integration tests; and
5. publishes both images to GHCR only after both validations pass.

Publication is protected by the `release` GitHub Environment. Only the
publication job receives `packages: write`. If both requested tags already
contain images with matching provenance, a rerun succeeds without rebuilding,
retesting, or republishing them. An existing tag with different content causes
the workflow to fail.

The Dockerfile is intentionally parameterized for this workflow. Local builds
must provide the digest-pinned dependency image, the downloaded package files
and checksums, source provenance, and image metadata as build arguments. See
`.github/workflows/svmultiphysics-image.yml` for the canonical invocation.

## Build the Ubuntu 26 dependency image locally

From the repository root, run:

```sh
docker build --platform linux/amd64 \
  -f Docker/ubuntu/dockerfile \
  -t svmp-libraries:ubuntu26 Docker/ubuntu
```

`DEPENDENCY_BUILD_TYPE` defaults to `Release` and accepts only `Release`,
`RelWithDebInfo`, or `Debug`. CMake-built dependencies use the selected type.
Boost and PETSc use their debug settings only for `Debug`. Both `Release` and
`RelWithDebInfo` use Boost's release variant and the same optimized PETSc flags
(`-g -O2`). The release process promotes only `Release` dependency images.

The Dockerfile downloads an x86-64 CMake binary, so this build targets
`linux/amd64`. On Apple Silicon, Docker runs it under emulation. The image
contains `/usr/local/share/svMultiPhysics/dependency-versions.txt` with resolved
source commits and installed Ubuntu, Conda, and Python package versions.

The `svmp-libraries:ubuntu26` tag is local unless explicitly pushed. Building it
does not change any published image.

## Legacy Docker Hub image

The legacy recipe remains available while its separate workflow is maintained:

```sh
docker build --platform linux/amd64 \
  -f Docker/solver/dockerfile \
  -t svmp-solver:local .
```

Its PETSc and Trilinos executables are located at:

```text
/build-petsc/svMultiPhysics-build/bin/svmultiphysics
/build-trilinos/svMultiPhysics-build/bin/svmultiphysics
```

The legacy image uses an Ubuntu 24.04-based dependency image and is independent
of the Ubuntu 26 package-based release images described above.
