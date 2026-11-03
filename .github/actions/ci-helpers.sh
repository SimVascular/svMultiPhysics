ci_build_svzerodsolver() {
  local jobs="$1"
  local cmake_args=()
  if [ "${SVMP_OPTIMIZED_CI:-false}" = true ]; then
    cmake_args+=(-DCMAKE_BUILD_TYPE=Release)
  fi

  git clone https://github.com/SimVascular/svZeroDSolver.git
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
