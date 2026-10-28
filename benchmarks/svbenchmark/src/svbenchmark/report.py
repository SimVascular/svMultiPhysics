"""
Reporting for a benchmark
"""

from typing import TextIO
from pathlib import Path
from dataclasses import asdict
from datetime import datetime
import json

from .execution import RunResult
from .evaluation import (
    AbsoluteTolerance,
    EvaluationCriterionResult,
    RelativeTolerance,
)
from .manifest import BenchmarkRunResult


def benchmark_to_json(result: BenchmarkRunResult) -> str:
    return json.dumps(asdict_handles_as_str(result))


def write_benchmark_json(benchmark_result: BenchmarkRunResult, json_fp: Path):
    with open(json_fp, "w", encoding="utf-8") as f:
            f.write(benchmark_to_json(benchmark_result))


def render_benchmark_summary(benchmark_result: BenchmarkRunResult, out: TextIO):
    """
    Renders a compact version of the benchmark
    """
    raise NotImplementedError


def render_evaluation_criterion_summary(
    name: str, eval_result: EvaluationCriterionResult, out: TextIO
) -> None:
    """Render one line representation of an evaluation result"""

    match eval_result.tolerance:
        case AbsoluteTolerance(value=value):
            tolerance = f"atol={value:.6e}"
        case RelativeTolerance(value=value):
            tolerance = f"rtol={value:.6e}"
        case None:
            tolerance = f"N/A"  # case: tracking value (not implemented yet)
        case unexpected:
            tolerance = f"ERROR"  # define tolerance here so the type checker know tolerance is str
            assert_never(unexpected)

    status = str(eval_result.status)

    # print to the TextIO
    print(
        f"{status:<4} "
        f"{name[:24]:<24} "
        f"{eval_result.value:>12.5e}  "
        f"{eval_result.reference:>12.5e}  "
        f"{tolerance}",
        file=out,
    )


def render_run_summary(run_result: RunResult, out: TextIO) -> None:
    """Render oneline representation of a run result"""
    raise NotImplementedError


def asdict_handles_as_str(obj):
    """
    Small asdict helper so that any non-defined serialization (Path, Datetime) are converted to string
    """
    return asdict(
        obj,
        dict_factory=lambda x: {
            k: (str(v) if isinstance(v, (Path, datetime)) else v) for k, v in x
        },
    )
