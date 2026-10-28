# GitHub Workflows for svMultiPhysics

This directory contains GitHub Actions workflows for tests, checks, Docker builds, documentation, and code coverage. Pull requests run the required fast-test workflow. Pushes to `main` run the full-test workflow instead of repeating both test configurations.

## Workflows

### 1. Fast Tests Workflow (`workflows/fast-tests.yml`)

**Purpose**: Provides the required pull-request validation with optimized Release builds on Ubuntu and macOS.

**Triggers**:
- Pull requests

**Key Features**:
- **Smart Change Detection**: Only runs integration tests when changes are made to important directories (`Code/`, `tests/`, `.github/`)
- **Multi-Platform Testing**: Runs tests on both Ubuntu (using Docker container) and macOS
- **Optimized Test Configuration**: Runs unit tests and rank-1 and rank-3 integration tests using Release builds
- **No Coverage Reporting**: Coverage and array-index checking are reserved for the full workflow
- **Required Summary**: Reports the single `Required Tests Summary` check used by branch protection, including when platform tests are skipped

**Jobs**:
1. `check-for-changes`: Determines if fast tests should run based on file changes in important directories. Uses `dorny/paths-filter@v3` to filter directories.
2. `test-ubuntu`: Runs tests in Ubuntu environment through `actions/test-ubuntu/action.yml`. Triggered only if `check-for-changes` returns `run_tests = true`
3. `test-macos`: Runs tests on macOS through `actions/test-macos/action.yml`. Triggered only if `check-for-changes` returns `run_tests = true`
4. `required-tests-summary`: Coordinates test results for branch protection. It succeeds when both platform jobs pass or when they are intentionally skipped, and fails when a required job fails or the execution decision cannot be determined.

### 2. Full Tests Workflow (`workflows/full-tests.yml`)

**Purpose**: Provides comprehensive validation with ranks 1, 3, and 4, coverage, and array-index checking on Ubuntu and macOS.

**Triggers**:
- Pushes to `main`
- Manual dispatch selecting `main` (jobs are skipped for other refs)
- Applying the `full-tests` label to a pull request; subsequent commits rerun the full tests while the label remains

For a risky pull request, a maintainer can apply the repository's `full-tests` label to test its merge result. The label must exist before this trigger can be used. Remove the label to stop full tests on later commits. PR-triggered runs use the pull request's security and cache scope; do not use a manual dispatch to check out PR code. The Ubuntu job saves its coverage report as a GitHub Actions artifact. Only pushes to `main` upload coverage to Codecov.

On pushes to `main`, only this workflow runs the solver test suite. The Fast Tests workflow is limited to pull requests to avoid repeating its unit tests and rank-1 and rank-3 integration tests after a merge.

### 3. Docker Solver Build Workflow (`workflows/docker_solver_build.yml`)

**Purpose**: Builds and publishes Docker images for the solver component.

**Triggers**:
- Pushes to the `main` branch

**Process**:
1. Sets up QEMU and Docker Buildx for multi-platform builds
2. Authenticates with DockerHub using stored credentials
3. Builds the Docker image using `Docker/solver/dockerfile`
4. Pushes the image to `simvascular/solver:latest`

### 4. Documentation Workflow (`workflows/documentation.yml`)

**Purpose**: Builds and deploys HTML documentation using Doxygen.

**Triggers**:
- All pushes and pull requests

**Process**:
1. Installs Doxygen
2. Generates documentation from `Documentation/Doxyfile`
3. Uploads documentation as an artifact
4. Deploys documentation to GitHub Pages (only for `main` branch)

### 5. Copyright Workflow (`workflows/copyright.yml`)

**Purpose**: Ensures all source files contain proper copyright headers.

**Triggers**:
- All pushes and pull requests

**Process**:
- Uses `kt3k/license_checker` to verify copyright headers in `.h` and `.cpp` files
- Fails the workflow if any files are missing proper copyright notices

## Testing Actions

### 1. Test Ubuntu Action (`actions/test-ubuntu/action.yml`)

**Purpose**: Shared Ubuntu build and test sequence used by both test workflows.

**Environment**: Uses `simvascular/libraries:ubuntu22` Docker container

**Steps**:
1. **Build svZeroDSolver**: Clones and builds the ZeroD solver dependency
2. **Build svMultiPhysics (Trilinos)**: Builds with Trilinos support and unit tests enabled
3. **Build svMultiPhysics (PETSc)**: Builds with PETSc support as an alternative
4. **Run Integration Tests**: Executes pytest tests in the `tests/` directory
5. **Run Unit Tests**: Executes CTest unit tests
6. **Generate Coverage**: Creates coverage reports only for Full Tests
7. **Upload Coverage**: Sends coverage data to Codecov only for Full Tests on `main`

### 2. Test macOS Action (`actions/test-macos/action.yml`)

**Purpose**: Shared macOS build and test sequence used by both test workflows.

**Steps**:
1. **Install Dependencies**: Uses Homebrew to install GCC, VTK, OpenBLAS, LAPACK, Mesa, OpenMPI, Qt, and lcov
2. **Install Miniconda**: Sets up Python environment
3. **Build svZeroDSolver**: Clones and builds the ZeroD solver
4. **Build svMultiPhysics**: Builds with unit tests enabled
5. **Install Test Dependencies**: Creates conda environment with Python testing libraries
6. **Run Integration Tests**: Executes pytest tests
7. **Run Unit Tests**: Executes CTest unit tests

## Configuration Files

### Codecov Configuration (`codecov.yml`)

Configures code coverage reporting with Codecov:
- Sets coverage status to "informational" for both project and patch coverage
- This means coverage reports won't block PRs but provide visibility into test coverage

## Usage

### For Contributors

1. **Making Changes**: When you push changes or create a pull request, the workflows will automatically run
2. **Test Requirements**: Ensure your changes don't break existing tests
3. **Copyright**: Make sure all new `.h` and `.cpp` files include proper copyright headers
4. **Documentation**: Update documentation if you add new features

### For Maintainers

1. **Branch Protection**: [Branch protection rules](https://github.com/SimVascular/svMultiPhysics/settings/branch_protection_rules/36012339) require successes from the following checks in order for a pull request to be merged.
   1. "check-license-lines" in `copyright.yml`
   2. "documentation" in `documentation.yml`
   3. "Required Tests Summary" in `fast-tests.yml`
   4. "review-policy" from the svMultiPhysics Review Policy app
2. **Docker Releases**: Docker images are automatically built and pushed when changes are merged to `main`
3. **Documentation**: Documentation is automatically deployed to GitHub Pages from the `main` branch

## Dependencies

### Required Secrets

The workflows require the following secrets to be configured in the repository settings:

- `CODECOV_TOKEN`: Token for uploading coverage reports to Codecov
- `DOCKER_PASSWORD`: Password for DockerHub authentication
- `GITHUB_TOKEN`: Automatically provided by GitHub for documentation deployment
