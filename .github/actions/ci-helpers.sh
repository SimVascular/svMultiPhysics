ci_resolve_svzerodsolver() {
  local manifest="${1:-${GITHUB_WORKSPACE:-.}/release.json}"
  local repository tag commit remote_commit

  repository=$(jq -er '.dependencies.svZeroDSolver.repository | select(type == "string")' "$manifest") || return 1
  tag=$(jq -er '.dependencies.svZeroDSolver.tag | select(type == "string")' "$manifest") || return 1
  commit=$(jq -er '.dependencies.svZeroDSolver.commit | select(type == "string")' "$manifest") || return 1

  if [[ "$repository" != SimVascular/svZeroDSolver || ! "$tag" =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ || ! "$commit" =~ ^[0-9a-f]{40}$ ]]; then
    echo "release.json must pin the official svZeroDSolver repository, a tag, and a full commit SHA" >&2
    return 1
  fi

  local remote="https://github.com/$repository.git"
  remote_commit=$(git ls-remote "$remote" "refs/tags/$tag^{}" | awk 'NR == 1 { print $1 }') || return 1
  if [[ -z "$remote_commit" ]]; then
    remote_commit=$(git ls-remote --refs "$remote" "refs/tags/$tag" | awk 'NR == 1 { print $1 }') || return 1
  fi
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
    cmake_args=(-DCMAKE_BUILD_TYPE=Release -DENABLE_COVERAGE=OFF -DENABLE_ARRAY_INDEX_CHECKING=OFF)
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
