# GitHub Workflows for svMultiPhysics

This directory contains GitHub Actions workflows that run various actions (tests, checks, docker builds, documentation, code coverage assessment) on the repository. These workflows are run whenever new code is added to the repository, either upon push to a fork of svMultiPhysics, or upon a pull request to SimVascular/svMultiPhysics. Upon pull request, these actions must succeed for pull request to be merged.

## Workflows

### 1. Tests Workflow (`workflows/tests.yml`)

**Purpose**: Runs comprehensive integration tests on both Ubuntu and macOS platforms when relevant code changes are detected. The integration tests can be found in the `tests/` directory and use `pytest`.

**Triggers**: 
- Pull requests
- Pushes to any branch

**Key Features**:
- **Smart Change Detection**: Only runs full test suite when changes are made to important directories (`Code/`, `tests/`, `.github/`)
- **Multi-Platform Testing**: Runs tests on both Ubuntu (using Docker container) and macOS
- **Coverage Reporting**: Generates and uploads code coverage reports to Codecov

**Jobs**:
1. `check-for-changes`: Determines if full test suite should run based on file changes in important directories. Uses `dorny/paths-filter@v3` to filter directories. 
2. `test-ubuntu`: Runs tests in Ubuntu environment through `actions/test-ubuntu/action.yml`. Triggered only if `check-for-changes` returns `run_full_tests = true`
3. `test-macos`: Runs tests on macOS through `actions/test-macos/action.yml`. Triggered only if `check-for-changes` returns `run_full_tests = true`
4. `all-tests`: Coordinates test results for branch protection rules. If `test-ubuntu` and `test-macos` run, returns success status only if both `test-ubuntu` and `test-macos` return success. If `test-ubuntu` and `test-macos` are skipped, always returns success. This test must pass for a pull request to be merged, as set by branch protection rules.

### 2. Docker Solver Build Workflow (`workflows/docker_solver_build.yml`)

**Purpose**: Builds and publishes Docker images for the solver component.

**Triggers**: 
- Pushes to the `main` branch

**Process**:
1. Sets up QEMU and Docker Buildx for multi-platform builds
2. Authenticates with DockerHub using stored credentials
3. Builds the Docker image using `Docker/solver/dockerfile`
4. Pushes the image to `simvascular/solver:latest`

### 3. Documentation Workflow (`workflows/documentation.yml`)

**Purpose**: Builds and deploys HTML documentation using Doxygen.

**Triggers**: 
- All pushes and pull requests

**Process**:
1. Installs Doxygen
2. Generates documentation from `Documentation/Doxyfile`
3. Uploads documentation as an artifact
4. Deploys documentation to GitHub Pages (only for `main` branch)

### 4. Copyright Workflow (`workflows/copyright.yml`)

**Purpose**: Ensures all source files contain proper copyright headers.

**Triggers**: 
- All pushes and pull requests

**Process**:
- Uses `kt3k/license_checker` to verify copyright headers in `.h` and `.cpp` files
- Fails the workflow if any files are missing proper copyright notices

## Testing Actions

### 1. Test Ubuntu Action (`actions/test-ubuntu/action.yml`)

**Purpose**: Comprehensive integration testing on Ubuntu with code coverage.

**Environment**: Uses `simvascular/libraries:ubuntu22` Docker container

**Steps**:
1. **Build svZeroDSolver**: Clones and builds the ZeroD solver dependency
2. **Build svMultiPhysics (Trilinos)**: Builds with Trilinos support, coverage, and unit tests enabled
3. **Build svMultiPhysics (PETSc)**: Builds with PETSc support as an alternative
4. **Run Integration Tests**: Executes pytest tests in the `tests/` directory
5. **Run Unit Tests**: Executes CTest unit tests
6. **Generate Coverage**: Creates coverage reports
7. **Upload Coverage**: Sends coverage data to Codecov

### 2. Test macOS Action (`actions/test-macos/action.yml`)

**Purpose**: Comprehensive testing on macOS with native dependencies.

**Steps**:
1. **Install Dependencies**: Uses Homebrew to install GCC, VTK, OpenBLAS, LAPACK, Mesa, OpenMPI, Qt, and lcov
2. **Install Miniconda**: Sets up Python environment
3. **Build svZeroDSolver**: Clones and builds the ZeroD solver
4. **Build svMultiPhysics**: Builds with coverage and unit tests enabled
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
   3. "All Tests" in `tests.yml`
2. **Docker Releases**: Docker images are automatically built and pushed when changes are merged to `main`
3. **Documentation**: Documentation is automatically deployed to GitHub Pages from the `main` branch

## Dependencies

### Required Secrets

The workflows require the following secrets to be configured in the repository settings:

- `CODECOV_TOKEN`: Token for uploading coverage reports to Codecov
- `DOCKER_PASSWORD`: Password for DockerHub authentication
- `GITHUB_TOKEN`: Automatically provided by GitHub for documentation deployment
