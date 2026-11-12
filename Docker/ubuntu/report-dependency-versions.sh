#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

macro_value() {
    awk -v name="$2" '$1 == "#define" && $2 == name { gsub(/"/, "", $3); print $3; found = 1; exit } END { if (!found) exit 1 }' "$1"
}

vtk_version=$(awk '/^set\(PACKAGE_VERSION / { gsub(/["\)]/, "", $2); print $2; found = 1; exit } END { if (!found) exit 1 }' \
    /vtk/lib/cmake/vtk-9.3/vtk-config-version.cmake)
boost_number=$(macro_value /boost/include/boost/version.hpp BOOST_VERSION)
lapack_configs=(/lapack/cmake/lapack-*/lapack-config-version.cmake)
test "${#lapack_configs[@]}" -eq 1
lapack_version=$(basename "$(dirname "${lapack_configs[0]}")")

printf 'Ubuntu: '
. /etc/os-release
printf '%s\n' "$PRETTY_NAME"
printf 'Ubuntu base image: %s\n' "$UBUNTU_BASE_IMAGE"
printf 'Build type: %s\n' "$DEPENDENCY_BUILD_TYPE"
printf 'CMake: %s\n' "$(/cmake/bin/cmake --version | awk 'NR == 1 {print $3}')"
printf 'Open MPI: %s\n' "$(/openmpi/bin/ompi_info --version | awk 'NR == 1 {print $NF}')"
printf 'VTK: %s\n' "$vtk_version"
printf 'Boost: %s.%s.%s\n' "$((boost_number / 100000))" "$(((boost_number / 100) % 1000))" "$((boost_number % 100))"
printf 'LAPACK: %s (installed CMake package)\n' "${lapack_version#lapack-}"
printf 'HDF5: %s.%s.%s\n' \
    "$(macro_value /hdf5/include/H5public.h H5_VERS_MAJOR)" \
    "$(macro_value /hdf5/include/H5public.h H5_VERS_MINOR)" \
    "$(macro_value /hdf5/include/H5public.h H5_VERS_RELEASE)"
printf 'HYPRE: %s\n' "$(macro_value /hypre/include/HYPRE_config.h HYPRE_RELEASE_VERSION)"
printf 'Trilinos: %s\n' "$(macro_value /trilinos/include/Trilinos_version.h TRILINOS_VERSION_STRING)"
printf 'PETSc: %s.%s.%s\n' \
    "$(macro_value /petsc/include/petscversion.h PETSC_VERSION_MAJOR)" \
    "$(macro_value /petsc/include/petscversion.h PETSC_VERSION_MINOR)" \
    "$(macro_value /petsc/include/petscversion.h PETSC_VERSION_SUBMINOR)"
printf 'GoogleTest: 1.14.0\n'
printf 'Miniconda installer SHA256: %s\n' "$(sha256sum /miniconda.sh | awk '{print $1}')"
printf 'Python: %s\n' "$(conda run -n svmultiphysics python -c 'import platform; print(platform.python_version())')"
printf 'Source commits:\n'
source_revisions=(/source-revisions/*)
if (( ${#source_revisions[@]} == 0 )); then
    echo 'Source revisions are missing from the dependency image' >&2
    exit 1
fi
for source in "${source_revisions[@]}"; do
    revision=$(cat "$source")
    if [[ ! "$revision" =~ ^[0-9a-f]{40}$ ]]; then
        echo "Invalid source revision in $source" >&2
        exit 1
    fi
    printf '  %s: %s\n' "${source##*/}" "$revision"
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
