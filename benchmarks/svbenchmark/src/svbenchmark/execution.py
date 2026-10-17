


from dataclasses import dataclass
from pathlib import Path
import subprocess
import time
from typing import Self

from ._expect import expect_str, expect_str_array

def run_solver(solver : SolverRun, context : SolveExecutionContext, cwd : Path) -> RunResult:
    """
    Runs the solver using the execution context and a working directory
    """
    command = solver.build_command(context)
    return run_command(command, cwd)

@dataclass(frozen=True, kw_only=True)
class SolveExecutionContext:
    """
    Provides the "context" or environment defintions that allow to be run. Specifically this is solver call paths
    and settings.

    This is used in combination with the dataset TOML defined definitions to realize an actual command call to run
    with the TOML defined intent
    """
    solver_executable: Path
    mpi_launcher: str = "mpirun"
    mpi_args: tuple[str, ...] = ()
    mpi_ranks: int = 1

@dataclass(frozen=True)
class SolverRun:
    name : str
    input : Path
    cwd : Path

    @classmethod
    def from_dict(cls, name : str, data : dict[str, object]) -> Self:
        input_str = expect_str(data['input'], 'input')
        cwd = Path(expect_str(data.get('cwd', ''), 'cwd'))
        
        return cls(name = name,
                   input = Path(input_str),
                   cwd = cwd)

    def build_command(
        self,
        context: SolveExecutionContext,
    ) -> Command:

        return Command(
            argv=(
                context.mpi_launcher,
                *context.mpi_args,
                "-np",
                str(context.mpi_ranks),
                str(context.solver_executable),
                str(self.input),
            ),
            cwd = self.cwd
        )
    def build_named_command(self, context : SolveExecutionContext):
        return NamedCommand(self.name, self.build_command(context))

@dataclass(frozen=True)
class Command:
    argv : tuple[str, ...]
    cwd : Path

    def to_str(self) -> str:
        return ' '.join(self.argv)

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> Self:
        """ 
        Creates a Command from a dictionary, raises errors if the argv does not exist or is incorrectly
        formatted
        """
        argv = expect_str_array(data["argv"], 'argv')
        cwd = Path(expect_str(data.get('cwd', ''), 'cwd'))
        return cls(argv=tuple(argv), cwd=cwd)


@dataclass(frozen=True)
class NamedCommand:
    name : str
    command : Command


    @classmethod
    def from_dict(cls, name : str, data: dict[str, object]) -> Self:
        """ 
        Creates a Command from a dictionary, raises errors if the argv does not exist or is incorrectly
        formatted
        """
        command = Command.from_dict(data)

        return cls(name, command)

@dataclass(frozen=True)
class RunResult:
    """
    Holds the information pertaining to a single command run result
    """
    run_cmd : str
    run_time : float
    returncode : int

@dataclass(frozen=True)
class CaptureResult:
    run_result : RunResult
    stdout : str
    stderr : str

def run_command(command : Command, cwd : Path) -> RunResult:
    """ 
    Runs a step based on its context which provides a run.
    """

    start_time = time.perf_counter()
    result = subprocess.run(command.argv, cwd=cwd.joinpath(command.cwd))
    end_time = time.perf_counter()

    return(RunResult(run_cmd=command.to_str(), 
                     run_time=end_time-start_time, 
                     returncode=result.returncode,)
                     )

def capture_command(command : Command, cwd : Path) -> CaptureResult:
    """ 
    Captures a command run, the command call does not print to terminal.
    """

    start_time = time.perf_counter()
    result = subprocess.run(command.argv, cwd=cwd.joinpath(command.cwd), capture_output=True, text=True)
    end_time = time.perf_counter()

    run_result = RunResult(run_cmd=command.to_str(), 
                     run_time=end_time-start_time, 
                     returncode=result.returncode,)

    return(CaptureResult(run_result=run_result, stdout=result.stdout, stderr=result.stderr)
                     )