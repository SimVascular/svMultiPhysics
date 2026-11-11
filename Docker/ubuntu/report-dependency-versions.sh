#!/usr/bin/env bash
set -euo pipefail

printf 'Ubuntu: '
. /etc/os-release
printf '%s\n' "$PRETTY_NAME"
printf 'Build type: %s\n' "$DEPENDENCY_BUILD_TYPE"
/cmake/bin/cmake --version | head -n 1
/openmpi/bin/ompi_info --version | head -n 1
grep '^set(PACKAGE_VERSION ' /vtk/lib/cmake/vtk-9.3/vtk-config-version.cmake | head -n 1
grep '^#define BOOST_LIB_VERSION ' /boost/include/boost/version.hpp
ls /lapack/cmake/lapack-*/lapack-config-version.cmake
grep -E '^#define H5_VERS_(MAJOR|MINOR|RELEASE) ' /hdf5/include/H5public.h
grep '^#define HYPRE_RELEASE_VERSION ' /hypre/include/HYPRE_config.h
grep '^#define TRILINOS_VERSION_STRING ' /trilinos/include/Trilinos_version.h
grep -E '^#define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR) ' /petsc/include/petscversion.h
printf 'GoogleTest: v1.14.0\n'
grep '^# VER:' /miniconda.sh
printf 'Miniconda installer: '
sha256sum /miniconda.sh
conda run -n svmultiphysics python --version
printf 'Source commits:\n'
for source in /source-revisions/*; do
    printf '  %s: %s\n' "${source##*/}" "$(cat "$source")"
done

if [[ "${1:-}" == --summary ]]; then
    exit 0
fi

printf '\nUbuntu packages:\n'
dpkg-query -W -f='${binary:Package}\t${Version}\n'
printf '\nConda packages (resolved versions):\n'
conda list -n svmultiphysics --explicit
printf '\nPython packages:\n'
conda run -n svmultiphysics python -m pip freeze
