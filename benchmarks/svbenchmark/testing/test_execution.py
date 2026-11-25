import pytest
from pathlib import Path
from svbenchmark.execution import run_command, Command
from subprocess import CalledProcessError


def test_run_cmd():
    pyfile_path = Path(__file__).parent / "fixtures" / "success_run_cmd.py"
    cmd = Command(argv=("python", pyfile_path.as_posix()), cwd=Path("./"))
    run_result = run_command(cmd, cmd.cwd)
    assert run_result.returncode == 0


def test_run_cmd_error():
    pyfile_path = Path(__file__).parent / "fixtures" / "error_run_cmd.py"
    cmd = Command(argv=("python", pyfile_path.as_posix()), cwd=Path("./"))
    with pytest.raises(CalledProcessError):
        run_result = run_command(cmd, cmd.cwd)
