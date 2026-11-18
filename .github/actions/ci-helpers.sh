ci_release_version() {
  local manifest="${1:-${GITHUB_WORKSPACE:-.}/release.json}"

  python3 - "$manifest" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as manifest:
    release = json.load(manifest)
parts = [release[key] for key in ("major", "minor", "patch", "build")]
if any(type(part) is not int or part < 0 for part in parts):
    raise SystemExit("release.json needs four non-negative integer version fields")
print(".".join(map(str, parts)))
PY
}

ci_resolve_remote_tag() {
  local remote="$1"
  local tag="$2"
  local commit

  commit=$(git ls-remote "$remote" "refs/tags/$tag^{}" | awk 'NR == 1 { print $1 }') || return 1
  if [[ -z "$commit" ]]; then
    commit=$(git ls-remote --refs "$remote" "refs/tags/$tag" | awk 'NR == 1 { print $1 }') || return 1
  fi
  if [[ ! "$commit" =~ ^[0-9a-f]{40}$ ]]; then
    echo "Tag $tag does not resolve to a full commit SHA at $remote" >&2
    return 1
  fi

  printf '%s\n' "$commit"
}

ci_commit_is_on_origin_main() {
  local commit="$1"

  git fetch --no-tags origin refs/heads/main:refs/remotes/origin/main || return 1
  git merge-base --is-ancestor "$commit" refs/remotes/origin/main
}

ci_resolve_ghcr_image() {
  local requested_image="$1"
  local digest inspect_output repository

  if [[ ! "$requested_image" =~ ^ghcr\.io/[a-z0-9._/-]+(:[A-Za-z0-9_.-]+|@sha256:[0-9a-f]{64})$ ]]; then
    echo 'Image must be a fully qualified GHCR tag or SHA-256 digest' >&2
    return 1
  fi
  inspect_output=$(docker buildx imagetools inspect "$requested_image") || return 1
  digest=$(awk '$1 == "Digest:" && !found++ { print $2 }' <<< "$inspect_output")
  if [[ ! "$digest" =~ ^sha256:[0-9a-f]{64}$ ]]; then
    echo "Could not resolve $requested_image to an image digest" >&2
    return 1
  fi
  if [[ "$requested_image" == *@* && "${requested_image##*@}" != "$digest" ]]; then
    echo "Requested digest ${requested_image##*@} resolved as $digest" >&2
    return 1
  fi
  repository=${requested_image%%@*}
  if [[ "$repository" == "$requested_image" ]]; then
    repository=${requested_image%:*}
  fi

  printf '%s %s %s@%s\n' "$repository" "$digest" "$repository" "$digest"
}

ci_svzerodsolver_identity() {
  local manifest="${1:-${GITHUB_WORKSPACE:-.}/release.json}"

  python3 - "$manifest" <<'PY'
import json
import re
import sys

with open(sys.argv[1], encoding="utf-8") as manifest:
    dependency = json.load(manifest)["dependencies"]["svZeroDSolver"]
repository, tag, commit = (dependency[key] for key in ("repository", "tag", "commit"))
if (repository != "SimVascular/svZeroDSolver"
        or not isinstance(tag, str)
        or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", tag)
        or not isinstance(commit, str)
        or not re.fullmatch(r"[0-9a-f]{40}", commit)):
    raise SystemExit("release.json must pin the official svZeroDSolver repository, a tag, and a full commit SHA")
print(tag, commit)
PY
}

ci_svzerodsolver_package_identity() {
  local manifest="${1:-${GITHUB_WORKSPACE:-.}/release.json}"
  local checksum commit tag

  read -r tag commit < <(ci_svzerodsolver_identity "$manifest") || return 1
  checksum=$(python3 - "$manifest" <<'PY'
import json
import re
import sys

with open(sys.argv[1], encoding="utf-8") as manifest:
    checksum = json.load(manifest)["dependencies"]["svZeroDSolver"]["ubuntu26_package_sha256"]
if not isinstance(checksum, str) or not re.fullmatch(r"[0-9a-f]{64}", checksum):
    raise SystemExit("release.json must pin the svZeroDSolver Ubuntu 26 package checksum")
print(checksum)
PY
  ) || return 1

  printf 'SimVascular/svZeroDSolver %s %s %s\n' "$tag" "$commit" "$checksum"
}

ci_resolve_svzerodsolver() {
  local manifest="${1:-${GITHUB_WORKSPACE:-.}/release.json}"
  local tag commit remote_commit pinned

  pinned=$(ci_svzerodsolver_identity "$manifest") || return 1
  read -r tag commit <<< "$pinned"

  local remote="https://github.com/SimVascular/svZeroDSolver.git"
  remote_commit=$(ci_resolve_remote_tag "$remote" "$tag") || return 1
  if [[ "$remote_commit" != "$commit" ]]; then
    echo "svZeroDSolver tag $tag resolves to ${remote_commit:-nothing}, not pinned commit $commit" >&2
    return 1
  fi

  printf '%s\n' "$commit"
}

ci_build_svzerodsolver() {
  local jobs="$1"
  local ref="${2:-}"
  local cmake_args=()
  if [ "${SVMP_OPTIMIZED_CI:-false}" = true ]; then
    cmake_args+=(-DCMAKE_BUILD_TYPE=Release)
  fi

  git clone https://github.com/SimVascular/svZeroDSolver.git
  if [ -n "$ref" ]; then
    git -C svZeroDSolver checkout --detach "$ref"
  fi
  mkdir svZeroDSolver/build
  (
    cd svZeroDSolver/build
    cmake .. "${cmake_args[@]}" \
      -DCMAKE_PROJECT_svZeroDSolver_INCLUDE="$GITHUB_WORKSPACE/Code/ThirdParty/eigen/share/eigen3/cmake/Eigen3Config.cmake" \
      -DFETCHCONTENT_SOURCE_DIR_EIGEN="$GITHUB_WORKSPACE/Code/ThirdParty/eigen"
    make -j"$jobs"
  )
}

ci_build_solver() {
  local build_dir="$1"
  local jobs="$2"
  shift 2

  local cmake_args=(-DENABLE_COVERAGE=ON -DENABLE_ARRAY_INDEX_CHECKING=ON)
  if [ "${SVMP_OPTIMIZED_CI:-false}" = true ]; then
    cmake_args=(-DCMAKE_BUILD_TYPE="${SVMP_CMAKE_BUILD_TYPE:-Release}" -DENABLE_COVERAGE=OFF -DENABLE_ARRAY_INDEX_CHECKING=OFF)
  fi

  mkdir "$build_dir"
  (
    cd "$build_dir"
    cmake "${cmake_args[@]}" -DENABLE_UNIT_TEST=ON "$@" ..
    make -j"$jobs"
  )
}

ci_run_rank_tests() {
  local rank="$1"
  local workers="$2"
  shift 2

  echo "Running rank-$rank tests with $workers worker(s)"
  SVMP_TEST_PROCS="$rank" conda run --no-capture-output -n svmultiphysics pytest -v -rfEP --durations=0 "$@" -n "$workers" --dist=loadfile
}

ci_run_integration_tests() {
  local single_rank_workers="$1"
  local multi_rank_workers="$2"
  local rank
  local ranks=()

  git lfs pull
  cd tests
  if [ "${SVMP_PARALLEL_INTEGRATION_TESTS:-false}" = true ]; then
    conda run -n svmultiphysics pip install pytest-xdist
    IFS=, read -r -a ranks <<< "${SVMP_TEST_PROCS-1,3,4}"
    for rank in "${ranks[@]}"; do
      if [ "$rank" = 1 ]; then
        ci_run_rank_tests "$rank" "$single_rank_workers"
      else
        ci_run_rank_tests "$rank" "$multi_rank_workers" --ignore=test_shell.py
      fi
    done
  else
    conda run --no-capture-output -n svmultiphysics pytest -v -rfEP --durations=0
  fi
}

ci_installed_svmultiphysics_paths() {
  local executable="${SVMP_SOLVER_EXECUTABLE:-}"
  local interface

  if [ -z "$executable" ]; then
    executable=$(find /usr/local/sv/svMultiPhysics -type f -name svmultiphysics -executable -print -quit) || return 1
  fi
  interface=$(find /usr/local/sv/svZeroDSolver -type f -name libsvzero_interface.so -print -quit) || return 1
  if [ ! -x "$executable" ]; then
    echo 'Could not locate the installed svMultiPhysics executable' >&2
    return 1
  fi
  if [ ! -f "$interface" ]; then
    echo 'Could not locate the installed svZeroDSolver interface library' >&2
    return 1
  fi

  printf '%s %s\n' "$executable" "$interface"
}

ci_check_installed_svmultiphysics_version() {
  local version="$1"
  local tag="$2"
  local commit="$3"
  local configuration="$4"
  local executable interface actual expected

  read -r executable interface < <(ci_installed_svmultiphysics_paths) || return 1
  actual=$("$executable" --version) || return 1
  expected="[svMultiPhysics] version $version tag $tag commit $commit configuration $configuration"
  if [ "$actual" != "$expected" ]; then
    echo "Expected: $expected" >&2
    echo "Actual:   $actual" >&2
    return 1
  fi
  printf '%s\n' "$actual"
}

ci_check_installed_svmultiphysics_libraries() {
  local executable interface

  read -r executable interface < <(ci_installed_svmultiphysics_paths) || return 1
  if ldd "$executable" | grep 'not found'; then
    echo 'svMultiPhysics has missing shared-library dependencies' >&2
    return 1
  fi
  if ldd "$interface" | grep 'not found'; then
    echo 'svZeroDSolver has missing shared-library dependencies' >&2
    return 1
  fi
}

ci_check_installed_svmultiphysics_mpi() {
  local rank
  local ranks=()

  IFS=, read -r -a ranks <<< "${SVMP_TEST_PROCS:-1,3,4}"
  for rank in "${ranks[@]}"; do
    mpirun --oversubscribe -np "$rank" /bin/true || return 1
  done
}

ci_run_installed_svmultiphysics_tests() {
  local result_file="$1"
  local conda_executable="${2:-conda}"
  local executable interface

  read -r executable interface < <(ci_installed_svmultiphysics_paths) || return 1
  mkdir -p "$(dirname "$result_file")" "${GITHUB_WORKSPACE}/svZeroDSolver/build/src/interface"
  ln -sfn "$interface" "${GITHUB_WORKSPACE}/svZeroDSolver/build/src/interface/libsvzero_interface.so"
  export SVMP_SOLVER_EXECUTABLE="$executable"
  (
    cd "${GITHUB_WORKSPACE}/tests"
    "$conda_executable" run --no-capture-output -n svmultiphysics \
      pytest -v -rfEP --durations=0 --ignore=test_release_identity.py \
      --junitxml="$result_file"
  )
}
