#!/usr/bin/env python3
"""
Script to replace old license headers with new SPDX-format headers.
Updates .cpp, .h, .py, and CMakeLists.txt files.
"""

import os
import sys
import re
from pathlib import Path

# SPDX headers for different file types
SPDX_HEADER_C_STYLE = """// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause
"""

SPDX_HEADER_HASH_STYLE = """# SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
# SPDX-License-Identifier: BSD-3-Clause
"""

SPDX_HEADER_FORTRAN_STYLE = """! SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
! SPDX-License-Identifier: BSD-3-Clause
"""

# Pattern to match the old license header (C-style /* ... */)
OLD_LICENSE_PATTERN_C = re.compile(
    r'^/\*\s*Copyright.*?SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE\.\s*\*/',
    re.DOTALL | re.MULTILINE
)

# Pattern to match the old license header (hash-style)
OLD_LICENSE_PATTERN_HASH = re.compile(
    r'^#\s*Copyright.*?SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE\.',
    re.DOTALL | re.MULTILINE
)

# Pattern to match the old license header (Fortran-style)
OLD_LICENSE_PATTERN_FORTRAN = re.compile(
    r'^!\s*Copyright.*?SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE\.\s*!',
    re.DOTALL | re.MULTILINE
)


def replace_license_header(file_path):
    """Replace old license header with SPDX format."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return False

    # Determine file type and appropriate header
    suffix = Path(file_path).suffix.lower()
    name = Path(file_path).name

    if suffix in ['.cpp', '.h', '.c', '.hpp', '.cxx', '.cc']:
        new_header = SPDX_HEADER_C_STYLE
        pattern = OLD_LICENSE_PATTERN_C
    elif suffix in ['.py'] or name == 'CMakeLists.txt':
        new_header = SPDX_HEADER_HASH_STYLE
        pattern = OLD_LICENSE_PATTERN_HASH
    elif suffix in ['.f', '.f90', '.for']:
        new_header = SPDX_HEADER_FORTRAN_STYLE
        pattern = OLD_LICENSE_PATTERN_FORTRAN
    else:
        print(f"Skipping unknown file type: {file_path}")
        return False

    # Check if file already has SPDX header
    if 'SPDX-License-Identifier' in content:
        print(f"Already has SPDX header: {file_path}")
        return False

    # Replace old header with new SPDX header
    new_content, num_replacements = pattern.subn(new_header.rstrip(), content)

    if num_replacements == 0:
        print(f"No old license header found in: {file_path}")
        return False

    if num_replacements > 1:
        print(f"WARNING: Multiple license headers found in: {file_path}")

    # Write the updated content
    try:
        with open(file_path, 'w', encoding='utf-8') as f:
            f.write(new_content)
        print(f"Updated: {file_path}")
        return True
    except Exception as e:
        print(f"Error writing {file_path}: {e}")
        return False


def find_source_files(root_dir):
    """Find all source files that need license header updates."""
    files_to_update = []

    # Walk through Code/Source directory
    source_dir = Path(root_dir) / 'Code' / 'Source'
    if not source_dir.exists():
        print(f"Warning: {source_dir} does not exist")
        return files_to_update

    for file_path in source_dir.rglob('*'):
        if file_path.is_file():
            suffix = file_path.suffix.lower()
            name = file_path.name

            if suffix in ['.cpp', '.h', '.c', '.hpp', '.cxx', '.cc', '.py', '.f', '.f90', '.for'] or name == 'CMakeLists.txt':
                files_to_update.append(str(file_path))

    return sorted(files_to_update)


def main():
    """Main function to update all license headers."""
    # Get repository root (script should be run from there)
    repo_root = Path(__file__).parent.absolute()

    print(f"Searching for source files in: {repo_root}")
    files_to_update = find_source_files(repo_root)

    print(f"\nFound {len(files_to_update)} files to check")
    print("=" * 80)

    updated_count = 0
    skipped_count = 0
    error_count = 0

    for file_path in files_to_update:
        result = replace_license_header(file_path)
        if result:
            updated_count += 1
        elif result is False:
            skipped_count += 1
        else:
            error_count += 1

    print("=" * 80)
    print(f"\nSummary:")
    print(f"  Updated: {updated_count}")
    print(f"  Skipped: {skipped_count}")
    print(f"  Errors:  {error_count}")
    print(f"  Total:   {len(files_to_update)}")

    return 0 if error_count == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
