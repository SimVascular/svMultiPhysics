"""
Reporting for a benchmark
"""

from pathlib import Path
from dataclasses import asdict
from datetime import datetime
import json

from .manifest import BenchmarkRunResult


def benchmark_to_json(result: BenchmarkRunResult) -> str:
    return json.dumps(asdict_handles_as_str(result))


def write_benchmark_json(benchmark_result: BenchmarkRunResult, json_fp: Path):
    with open(json_fp, "w", encoding="utf-8") as f:
        f.write(benchmark_to_json(benchmark_result))


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
