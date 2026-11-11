"""Report third-party versions from the image used to link svMultiPhysics."""

import glob
import re
import subprocess
import sys
from pathlib import Path


def cmake_version(patterns, variables):
    for pattern in patterns:
        for filename in sorted(glob.glob(pattern)):
            contents = Path(filename).read_text(errors="replace")
            for variable in variables:
                match = re.search(
                    rf"(?im)^\s*set\s*\(\s*{variable}\s+\"?([0-9][^\s\"\)]+)",
                    contents,
                )
                if match:
                    return match.group(1)
    return "not recorded"


def shared_library_version(name):
    for filename in sorted(Path("/lapack").glob(f"lib{name}.so*")):
        match = re.search(r"\.so\.(\d+(?:\.\d+)+)", filename.resolve().name)
        if match:
            return match.group(1)
    return "not recorded"


def command_version(command, pattern):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError:
        return "not recorded"
    if result.returncode == 0:
        match = re.search(pattern, result.stdout, re.MULTILINE)
        if match:
            return match.group(1)
    return "not recorded"


def main():
    if len(sys.argv) != 2 or sys.argv[1] not in {"default", "trilinos"}:
        raise SystemExit("usage: report-linked-dependencies.py default|trilinos")
    configuration = sys.argv[1]
    versions = [
        (
            "Open MPI",
            command_version(
                ["/openmpi/bin/ompi_info", "--version"],
                r"Open MPI v?(\d+(?:\.\d+)+)",
            ),
        ),
        ("BLAS", shared_library_version("blas")),
        ("LAPACK", shared_library_version("lapack")),
        (
            "VTK (CommonCore, CommonDataModel, CommonMisc, FiltersGeometry, IOXML)",
            cmake_version(
                [
                    "/vtk/lib/cmake/vtk-*/VTKConfigVersion.cmake",
                    "/vtk/lib64/cmake/vtk-*/VTKConfigVersion.cmake",
                ],
                ["PACKAGE_VERSION"],
            ),
        ),
        (
            "zlib",
            command_version(
                ["dpkg-query", "-W", "-f=${Version}", "zlib1g"],
                r"^(\S+)",
            ),
        ),
    ]
    if configuration == "trilinos":
        versions.append(
            (
                "Trilinos",
                cmake_version(
                    [
                        "/trilinos/lib/cmake/Trilinos/TrilinosConfigVersion.cmake",
                        "/trilinos/lib/cmake/Trilinos/TrilinosConfig.cmake",
                        "/trilinos/lib64/cmake/Trilinos/TrilinosConfigVersion.cmake",
                        "/trilinos/lib64/cmake/Trilinos/TrilinosConfig.cmake",
                    ],
                    ["PACKAGE_VERSION", "Trilinos_VERSION"],
                ),
            )
        )

    print("\n#### Third-party link inputs\n")
    print(
        "The solver target links these libraries for this configuration. "
        "Versions are read from the pinned dependency image, not the Ubuntu "
        "packages installed during validation.\n"
    )
    print("| Library | Version in pinned build image |")
    print("| --- | --- |")
    for library, version in versions:
        print(f"| {library} | {version} |")
    print(
        "\nThe solver also links bundled METIS, ParMETIS, TetGen, and TinyXML; "
        "their source is identified by the source commit above. PETSc is not "
        "enabled for these packages."
    )


if __name__ == "__main__":
    main()
