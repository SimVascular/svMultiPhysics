# Docker images

Docker Engine is required to build and run these images. On macOS, Docker Desktop
provides the engine; on Linux, Docker Engine can be installed directly. Building
the dependency image requires substantial CPU, memory, and disk space.

## Image roles

- `Docker/ubuntu/dockerfile` builds an Ubuntu 26.04 image with the libraries and
  tools needed to build and test svMultiPhysics. It does not contain a compiled
  svMultiPhysics solver. The manually triggered dependency-image workflow stages
  the image in GitHub Container Registry, then validates it before promotion.
- `Docker/solver/dockerfile` builds svMultiPhysics from the SimVascular GitHub
  repository in PETSc and Trilinos configurations. Both its builder and final
  stages use the published `simvascular/libraries:latest` dependency image, not
  the local Ubuntu 26 image. The repository's Docker publishing workflow uses
  this Dockerfile and publishes `simvascular/solver:latest` on pushes to `main`.

Despite its tag name, `simvascular/solver:latest` is an svMultiPhysics image;
it is not the separate Fortran CFD solver. The separately published
`simvascular/svmultiphysics:latest` also contains svMultiPhysics PETSc and
Trilinos builds, but the current repository workflow does not publish that tag.

## Build the Ubuntu 26 dependency image

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
(`-g -O2`); this is intentional. The release process promotes only `Release`
dependency images to `rc` or final. Other build-type spellings fail before
dependencies are installed.

The current Dockerfile downloads an x86-64 CMake binary, so this build targets
`linux/amd64`. On Apple Silicon, Docker runs it under emulation. The Dockerfile
uses `bsdtar` in archive-extracting stages under arm64-to-amd64 emulation to
avoid a GNU tar failure.

The `svmp-libraries:ubuntu26` tag is local to your machine unless you explicitly
push it to a registry. Building it does not change either published Docker Hub
image or the base image used by `Docker/solver/dockerfile`.

The Dockerfile defaults to the Ubuntu 26.04 base tag; the dependency-image
workflow resolves its digest before building and records the exact base image
reference in the version inventory. Local builds can pass `UBUNTU_BASE_IMAGE`
with a digest to get the same provenance. Source-built libraries use release
tags where available (or Git commits otherwise). The Conda environment
requests Python 3.12; Ubuntu, Conda, and pip package versions are not pinned.
The image contains `/usr/local/share/svMultiPhysics/dependency-versions.txt`
with the resolved source commits and installed Ubuntu, Conda, and Python package
versions.
The dependency-image workflow also uploads this inventory as an artifact.
Changing any pin creates a new image candidate and requires the default,
Trilinos, and PETSc validation jobs to pass again. Do not assume that a pin
selected from an older image reproduces its exact binaries when the older build
did not record the source commit.

## Build the existing solver image

The solver Dockerfile's current steps clone svMultiPhysics and build its PETSc
and Trilinos configurations using `simvascular/libraries:latest`. From the
repository root, build it without changing those steps:

```sh
docker build --platform linux/amd64 \
  -f Docker/solver/dockerfile \
  -t svmp-solver:local .
```

The PETSc and Trilinos executables are in:

```text
/build-petsc/svMultiPhysics-build/bin/svmultiphysics
/build-trilinos/svMultiPhysics-build/bin/svmultiphysics
```

This solver build uses the published Ubuntu 24.04-based dependency image. It
does not validate the new Ubuntu 26 dependency image. To test Ubuntu 26, run
builds and tests inside `svmp-libraries:ubuntu26` instead; that does not require
changing the existing solver Dockerfile or its publishing workflow.

## Run an image

To open a shell and mount a host directory at `/workspace`:

```sh
docker run --rm -it --platform linux/amd64 \
  -v /absolute/path/on/host:/workspace svmp-libraries:ubuntu26
```

Files written under `/workspace` remain on the host after the container exits.
