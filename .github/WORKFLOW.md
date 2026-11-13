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

### 3. Ubuntu 26 Dependency Image Workflow (`workflows/dependency-image.yml`)

**Purpose**: Build, validate, and stage the Ubuntu 26 AMD64 dependency image in GitHub Container Registry (GHCR). This image is a build environment for svMultiPhysics, not the published solver image. Its four-part version comes from the numeric `major`, `minor`, `patch`, and `build` fields in `release.json`.

**Trigger**: Manual `workflow_dispatch` from the Actions tab on `main`. GitHub requires the workflow file to exist on the repository's default branch before manual dispatch is available. The selected source commit must belong to the current repository's `main` branch.

**Access control**: The upstream repository's GitHub Actions policy restricts this workflow to members of the Refactoring Team. The workflow file itself does not restrict who can click **Run workflow**; keep the policy scoped to this workflow and the `workflow_dispatch` event. Do not treat the documentation or a job-level `if` as an access control. See [Control workflow execution](https://docs.github.com/en/actions/how-tos/administer/control-workflow-execution) and [Manually run a workflow](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow).

**Manual inputs**:

| Input | Default | Meaning |
| --- | --- | --- |
| `release_stage` | `dev` | `dev` builds from scratch and runs validation. `rc` promotes the matching `dev` image; `final` promotes the matching `rc` image. Promotion does not rebuild. |
| `source_commit` | Current workflow ref's commit for `dev` | Full 40-character svMultiPhysics commit SHA to build or promote. Required for `rc` and `final`; those stages also require the preceding stage's Git tag to point to this commit. |
| `build_type` | `Release` | Dependency CMake build type: `Release`, `RelWithDebInfo`, or `Debug`; other spellings fail. Boost uses its release variant for both `Release` and `RelWithDebInfo`; PETSc intentionally uses the same optimized `-g -O2` flags for both. Only `Release` can be promoted to `rc` or `final`. |

The workflow resolves the Ubuntu 26.04 base image digest before building and records the exact base reference in the dependency inventory.

The workflow creates a versioned Git tag and GHCR image tag only after a `dev` build passes its image and solver tests. With `major`, `minor`, `patch`, and `build` set to `1`, `0`, `0`, and `1`, Release tags are `v1.0.0.1-dev`, `v1.0.0.1-rc`, and `v1.0.0.1`; non-Release development tags add the lowercase build type. The source image is initially staged under a commit-based tag so failed tests cannot publish a versioned `dev` tag. If the `dev` Git tag and image already exist for the same commit, version, and build type, a rerun succeeds without rebuilding or retesting. If the tag or image belongs to another commit, increment and commit the `build` field in `release.json` before starting a new build. The `rc` and `final` stages resolve the preceding image tag to a digest once, verify its source commit, version, and build type, then promote by digest rather than tag. This prevents a tag change during promotion, but does not prove the tag was unchanged since the preceding stage's tests. Each stage checks that its versioned tag has the same digest as its source and records the digest in the run summary. Repeating an existing promotion verifies the image and leaves it unchanged.

### 4. Ubuntu 26 Dependency Image Tests (`workflows/dependency-image-tests.yml`)

**Purpose**: Validate the staged dependency image before it receives a versioned `dev` tag. This reusable workflow is called by the image workflow and cannot be started directly from the Actions tab.

It receives the image repository, immutable image digest, svMultiPhysics source commit, and build type from the caller. It reads the svZeroDSolver tag and commit from that source commit's `release.json` and verifies that the tag still points to the pinned commit. On separate Ubuntu runners, it builds the default, Trilinos, and PETSc solver configurations against the exact image, builds the pinned svZeroDSolver commit, and runs unit and integration tests. Successful job summaries record the image digest, source commits, configuration, and build type; outcomes are shown in GitHub Actions. A failure prevents the caller from creating the versioned `dev` Git and image tags.

### 5. Ubuntu 26 Package Build (`workflows/ubuntu26-packages.yml`)

**Purpose**: Build default and Trilinos `.deb` packages using a selected Ubuntu 26 dependency image. The workflow resolves the image to a digest before building, and the package version comes from `release.json`.

**Trigger**: Manual `workflow_dispatch`. The workflow must be on the repository's default branch before it can be started manually.

**Manual inputs**:

| Input | Default | Meaning |
| --- | --- | --- |
| `source_ref` | Selected workflow ref's commit | Full 40-character commit SHA or complete `-dev`, `-rc`, or final Git tag, including the `v` prefix (for example, `v1.0.0.2-dev`, `v1.0.0.2-rc`, or `v1.0.0.2`). The selected commit supplies the solver, version, and test fixtures and must belong to the current repository's `main`. |
| `dependency_image` | `ghcr.io/simvascular/svmultiphysics-dependencies:1.0.0.1-dev` | Fully qualified GHCR image tag or digest. It may identify a dependency image from an earlier release. The workflow resolves it to one digest for both package jobs. |
| `build_type` | `Release` | Solver CMake build type: `Release`, `RelWithDebInfo`, or `Debug`. Only a successful `Release` run can create the final Git tag. |

The package filename and embedded version use the four-part version from `release.json` without a `-dev` suffix; selecting a different build type does not change that filename.

The selected solver commit's `release.json` supplies the final four-part version. The final `v<version>` tag need not exist when the workflow starts. If it already exists, it must point to the selected commit; the workflow never moves it. The package job checks out the selected commit, passes `SV_RELEASE_PACKAGE=ON` and `SV_TAG_AFTER_VALIDATION=ON` to CMake, and checks that the build-tree solver reports the intended final tag and commit. CPack generates the final package filename and version metadata; the workflow verifies both. The default and Trilinos configurations each run unit and integration tests before packaging. A separate clean Ubuntu 26.04 container then verifies each installed solver's `--version` output and runs its integration tests against that executable.

Only after the default and Trilinos validation jobs pass does a manual Release run on `main` create the final Git tag on the tested commit and upload the two `.deb` packages, their SHA-256 checksums, and a full third-party dependency report to a GitHub Release under that tag. The release notes record the workflow commit, built solver commit, resolved dependency image digest, pinned svZeroDSolver revision, and key dependency versions. The release starts as a draft and becomes public only after all assets are present. A rerun verifies existing assets and provenance and refuses to replace different content. The dependency image is selected independently: solver releases may reuse a previously built dependency image.

**Access control**: Only members of the Refactoring Team may run this workflow in the upstream repository. A repository administrator must enforce this with an Actions policy scoped to this workflow and the `workflow_dispatch` event; the workflow file alone does not restrict who can start a manual run.

Each package job builds one solver configuration, runs its unit and integration tests with the pinned svZeroDSolver source, checks the Debian package version, and uploads the package and SHA-256 checksum as an Actions artifact. The default job also captures the dependency image's full version inventory, including Ubuntu, Conda, and Python packages. A dependent validation job for each configuration runs in a clean, digest-pinned Ubuntu 26.04 container. It verifies the pinned svZeroDSolver release package against the checksum in `release.json`, installs both packages, checks the installed solver identity, and runs the integration suite against the installed solver using Python 3.12. Test results are uploaded as separate artifacts. After validation, the Ubuntu packages, checksums, and dependency inventory are also published to GitHub Releases for external testing. SimTK publication remains separate.

### 6. Docker Solver Build Workflow (`workflows/docker_solver_build.yml`)

**Purpose**: Builds and publishes Docker images for the solver component.

**Triggers**:
- Pushes to the `main` branch

**Process**:
1. Sets up QEMU and Docker Buildx for multi-platform builds
2. Authenticates with DockerHub using stored credentials
3. Builds the Docker image using `Docker/solver/dockerfile`
4. Pushes the image to `simvascular/solver:latest`

### 7. Documentation Workflow (`workflows/documentation.yml`)

**Purpose**: Builds and deploys HTML documentation using Doxygen.

**Triggers**:
- All pushes and pull requests

**Process**:
1. Installs Doxygen
2. Generates documentation from `Documentation/Doxyfile`
3. Uploads documentation as an artifact
4. Deploys documentation to GitHub Pages (only for `main` branch)

### 8. Copyright Workflow (`workflows/copyright.yml`)

**Purpose**: Ensures all source files contain proper copyright headers.

**Triggers**:
- All pushes and pull requests

**Process**:
- Uses `kt3k/license_checker` to verify copyright headers in `.h` and `.cpp` files
- Fails the workflow if any files are missing proper copyright notices

### 9. Review Policy Workflow (`workflows/review-policy.yml`)

The required `review-policy` check normally requires two approvals from eligible reviewers with write access. A Refactoring Team member can apply one of the following labels to allow one approval from an eligible member of that team other than the PR author:

| Label | Ordinary source changes | Workflow, action, or review-policy changes |
| --- | --- | --- |
| `OOP Refactor` | One approval | Two approvals |
| `software engineering` | One approval | One approval |

The `software engineering` exception includes changes to the review-policy script and its own workflows. When both labels are present, each exception retains its own requirements, and either authorized applicable exception can be used. Labels applied by someone outside the Refactoring Team are removed by the policy workflow. The PR author's approval never counts, and an eligible review requesting changes blocks the check. A potentially truncated changed-file list still requires two approvals.

The workflow reevaluates open PRs when labels, PRs, or reviews change, when policy changes reach `main`, and on its scheduled refresh. It runs the policy code from `main`, so the PR introducing a policy change is evaluated under the policy already on `main`.

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
