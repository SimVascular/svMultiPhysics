import json
from pathlib import Path
from unittest.mock import mock_open, patch

from svbenchmark.evaluation import (
    AbsoluteTolerance,
    EvaluationCriterionResult,
    EvaluationStatus,
)
from svbenchmark.execution import (
    ExecutionEnvironment,
    RunResult,
    SolveExecutionContext,
    SolverRunResult,
)
from datetime import datetime

from svbenchmark.manifest import BenchmarkRunResult
from svbenchmark.report import benchmark_to_json


def generate_benchmark_result() -> BenchmarkRunResult:
    """
    Creates the benchmark run result outside the test_report_json. This is mainly so that it is easy to manually
    generate the reference fixture.

    """
    solver_context = SolveExecutionContext(solver_executable=Path("fake_path"))
    execution_environment = ExecutionEnvironment(
        system="FakeSystem",
        release="FakeRelease",
        version="FakeVersion",
        machine="FakeMachine",
        processor="FakeProcessor",
        cpus=1,
    )

    FIXED_TIME = datetime(2026, 10, 27, 12, 0, 0)

    run_result = RunResult(
        run_cmd="",
        run_time_seconds=1,
        returncode=0,
        execution_environment=execution_environment,
        utc=FIXED_TIME,
    )
    eval_result = EvaluationCriterionResult(
        reference=0.0,
        tolerance=AbsoluteTolerance(0.001),
        value=0.0001,
        status=EvaluationStatus.PASSED,
        run_result=run_result,
    )
    solver_result = SolverRunResult.from_run_result(run_result, solver_context)

    benchmark_result = BenchmarkRunResult(
        name="Test",
        base_directory=Path("FakeDirectory"),
        preprocess={},
        solver={"fake_run": solver_result},
        evaluation={"fake_eval": eval_result},
        postprocess={},
        cleanup={},
    )
    return benchmark_result


def test_report_json():
    """
    Create a minimal valid benchmark result, make sure that the json has all the necessary componentns
    """

    benchmark_result = generate_benchmark_result()

    expected = json.loads(
        (Path(__file__).parent / "fixtures" / "benchmark_result.json").read_text()
    )

    actual = json.loads(benchmark_to_json(benchmark_result))

    assert actual == expected
