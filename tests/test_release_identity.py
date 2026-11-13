import json
import re
import shutil
import subprocess
from pathlib import Path

from .conftest import cpp_exec


def test_release_manifest_and_tag_must_match(tmp_path):
    source = tmp_path / "source"
    module_dir = source / "Code" / "CMake"
    module_dir.mkdir(parents=True)
    repo = Path(__file__).resolve().parents[1]
    shutil.copyfile(repo / "Code/CMake/SimVascularReleaseIdentity.cmake",
                    module_dir / "SimVascularReleaseIdentity.cmake")
    (source / "release.json").write_text(json.dumps({
        "major": 1, "minor": 0, "patch": 0, "build": 2,
    }))
    check = source / "check.cmake"
    check.write_text(
        f'include("{module_dir / "SimVascularReleaseIdentity.cmake"}")\n'
        'if(NOT SV_PACKAGE_VERSION STREQUAL "1.0.0.2")\n'
        '  message(FATAL_ERROR "Incorrect release version")\n'
        'endif()\n'
    )
    subprocess.run(["git", "init", "-q", str(source)], check=True)
    subprocess.run(["git", "-C", str(source), "add", "."], check=True)
    subprocess.run(["git", "-C", str(source), "-c", "user.name=Test",
                    "-c", "user.email=test@example.com", "commit", "-qm", "fixture"], check=True)
    command = ["cmake", "-DSV_ENABLE_DISTRIBUTION=ON",
               "-DSV_RELEASE_PACKAGE=ON", "-P", str(check)]

    developer_build = subprocess.run(
        ["cmake", "-P", str(check)], capture_output=True, text=True,
    )
    assert developer_build.returncode == 0, developer_build.stderr

    local_package = subprocess.run(
        ["cmake", "-DSV_ENABLE_DISTRIBUTION=ON", "-P", str(check)],
        capture_output=True, text=True,
    )
    assert local_package.returncode == 0, local_package.stderr

    missing_tag = subprocess.run(command, capture_output=True, text=True)
    assert missing_tag.returncode != 0
    assert "must exist and point to source commit" in missing_tag.stderr

    subprocess.run(["git", "-C", str(source), "tag", "v1.0.0.2"], check=True)
    assert subprocess.run(command, capture_output=True, text=True).returncode == 0
    (source / "release.json").write_text((source / "release.json").read_text() + "\n")
    assert subprocess.run(
        ["cmake", "-DSV_ENABLE_DISTRIBUTION=ON", "-P", str(check)],
        capture_output=True, text=True,
    ).returncode == 0
    dirty_source = subprocess.run(command, capture_output=True, text=True)
    assert dirty_source.returncode != 0
    assert "clean source tree" in dirty_source.stderr
    subprocess.run(["git", "-C", str(source), "restore", "release.json"], check=True)

    subprocess.run(["git", "-C", str(source), "-c", "user.name=Test",
                    "-c", "user.email=test@example.com", "commit", "--allow-empty",
                    "-qm", "later commit"], check=True)
    moved_source = subprocess.run(command, capture_output=True, text=True)
    assert moved_source.returncode != 0
    assert "must exist and point to source commit" in moved_source.stderr


def test_ubuntu_package_name_and_metadata(tmp_path):
    repo = Path(__file__).resolve().parents[1]
    cmake_lists = (
        'cmake_minimum_required(VERSION 3.22)\n'
        'project(release_identity C CXX)\n'
        'set(APPLE FALSE)\n'
        'set(LINUX TRUE)\n'
        'set(LSB_DISTRIB Ubuntu)\n'
        'set(LSB_VERSION 26.04)\n'
        'set(SV_INSTALL_RUNTIME_DIR bin)\n'
        'set(SV_MAJOR_VERSION 1)\n'
        'set(SV_MINOR_VERSION 0)\n'
        'set(SV_PATCH_VERSION 0)\n'
        'set(SV_PACKAGE_VERSION 1.0.0.2)\n'
        'set(SV_FULL_VERSION 1.0.0.2)\n'
        f'add_subdirectory("{repo / "Distribution"}" distribution)\n'
    )
    (tmp_path / "CMakeLists.txt").write_text(cmake_lists)
    for release_package in (False, True):
        for configuration in ("default", "trilinos"):
            build = tmp_path / (configuration + ("-release" if release_package else "-local"))
            result = subprocess.run(
                ["cmake", "-S", str(tmp_path), "-B", str(build),
                 f"-DSV_RELEASE_PACKAGE={'ON' if release_package else 'OFF'}",
                 f"-DSV_USE_TRILINOS={'ON' if configuration == 'trilinos' else 'OFF'}"],
                capture_output=True, text=True,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            config = (build / "CPackConfig.cmake").read_text()
            suffix = "" if release_package else "-local"
            version = "1.0.0.2" if release_package else "1.0.0.2~local"
            assert (f'set(CPACK_PACKAGE_FILE_NAME '
                    f'"svMultiPhysics-Ubuntu-26-{configuration}-1.0.0.2{suffix}")') in config
            assert f'set(CPACK_DEBIAN_PACKAGE_VERSION "{version}")' in config


def test_solver_reports_tag_and_commit():
    result = subprocess.run(
        [cpp_exec, "--version"], capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    assert re.fullmatch(
        r"\[svMultiPhysics\] version \S+ commit ([0-9a-f]{40}|unknown)\n",
        result.stdout,
    )


def test_normal_solver_run_does_not_report_release_identity(tmp_path):
    result = subprocess.run(
        [cpp_exec, str(tmp_path / "missing.xml")],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode != 0
    assert "[svMultiPhysics] version " not in result.stdout
