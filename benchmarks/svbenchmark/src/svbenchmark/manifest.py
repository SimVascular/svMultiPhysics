

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import tomllib
from typing import Callable, Mapping, Self

from .evaluation import EvaluationCriterion, EvaluationCriterionResult
from .execution import NamedCommand, RunResult, SolverRun
from ._expect import expect_int, expect_str, expect_table

@dataclass(frozen=True, kw_only=True)
class Benchmark:
    """
    Holds the manifest of declarative Benchmark definition
    """
    name : str
    schema : int
    reference_source : str
    base_directory : Path
    preprocess : tuple[NamedCommand,...]  = tuple()
    solver : tuple[SolverRun,...]
    postprocess : tuple[NamedCommand,...] = tuple()
    evaluation : tuple[EvaluationCriterion,...]
    cleanup : tuple[NamedCommand,...]     = tuple()

    @classmethod
    def from_toml_dict(cls, base_directory : Path, data : dict[str, object]) -> Self:
        """
        Creates a Benchmark from a TOML dict. The base_dir is the directory where the paths/files in the benchmark
        are relative to. For a TOML file, this is just the TOML parent directory itself
        """

        # Define helper functions for repetitive extraction + error handling
        def extract_dictionary(key : str, data : dict[str, object]) -> dict[str, object]:
            """ Function that extracts a dictionary, raises a ValueError if incorrectly formatted"""
            dict_data = data[key]
            return expect_table(dict_data, key)

        def unpack_to_tuple[T](
            data: dict[str, object],
            from_dict_fn: Callable[[str, dict[str, object]], T],
        ) -> tuple[T, ...]:
            """Unpack named TOML tables into typed tuples of object, creating the object with a passed constructor."""
            l = list[T]()
            for k in data:
                l.append(from_dict_fn(k, expect_table(data[k], k)))

            return tuple(l)

        # parse the toml_dictionary
        base_directory = base_directory.resolve()
        name = expect_str(data['name'], "name")
        schema = expect_int(data['schema'], "schema")
        reference_source = expect_str(data['reference_source'], "reference_source")

        ###  go through each phase, unpacking as needed
        # preprocess
        if "preprocess" in data:
            preprocess_data = unpack_to_tuple(extract_dictionary('preprocess', data), NamedCommand.from_dict)
        else:
            preprocess_data = tuple()

        # solver (required)
        if "solver" not in data:
            raise ValueError("'Solver' table needs to be in TOML")
        solver_data = unpack_to_tuple(extract_dictionary('solver', data), SolverRun.from_dict)
        

        # postprocess
        if "postprocess" in data:
            postprocess_data = unpack_to_tuple(extract_dictionary('postprocess', data), NamedCommand.from_dict)
        else:
            postprocess_data = tuple()

        # evaluation
        if "evaluation" not in data:
            raise ValueError("'Evaluation' table needs to be in TOML")
        evaluation_data = unpack_to_tuple(extract_dictionary('evaluation', data), EvaluationCriterion.from_dict)

        # cleanup
        if "cleanup" in data:
            cleanup_data = unpack_to_tuple(extract_dictionary('cleanup', data), NamedCommand.from_dict)
        else:
            cleanup_data = tuple()

        return cls(schema = schema,
                   name = name, 
                   reference_source = reference_source,
                   base_directory = base_directory, 
                   preprocess = preprocess_data,
                   solver = solver_data, 
                   postprocess = postprocess_data, 
                   evaluation = evaluation_data, 
                   cleanup = cleanup_data)

    @classmethod
    def read_toml_file(cls, filepath : Path) -> Self:
        with filepath.open('rb') as f:
            data = tomllib.load(f)

        return cls.from_toml_dict(filepath.parent, data)

def load_benchmark(path : Path | str) -> Benchmark:
    """
    Basically the Benchmark.read_toml_file except it normalizes to Path (Paths are used inside svbenchmark)
    """
    path = Path(path) # cast to path, if not already
    return Benchmark.read_toml_file(path)

@dataclass(frozen=True)
class BenchmarkRunResult:
    """
    A complete benchmark result needs all phases populated (they can be empty dictionaries)
    """
    # benchmark information
    name : str
    base_directory : Path

    # results
    preprocess : Mapping[str, RunResult]
    solver : Mapping[str, RunResult]
    postprocess : Mapping[str, RunResult]
    evaluation : Mapping[str, EvaluationCriterionResult]
    cleanup : Mapping[str, RunResult]

    @classmethod
    def from_results(cls, benchmark : Benchmark, 
                     preprocess : Mapping[str, RunResult], 
                     solver : Mapping[str, RunResult], 
                     postprocess : Mapping[str, RunResult], 
                     evaluation : Mapping[str, EvaluationCriterionResult], 
                     cleanup:Mapping[str, RunResult]) -> Self:
        return cls(
            name = benchmark.name,
            base_directory = benchmark.base_directory,
            preprocess = preprocess,
            solver = solver,
            postprocess = postprocess,
            evaluation = evaluation,
            cleanup = cleanup
        )

def get_utc_earliest(benchmark_result : BenchmarkRunResult) -> datetime:
    """
    Gets the earliest datetime (UTC) present in the benchmark
    """
    # loop through the run call (up to solvers, since those are always required), return the first one encountered
    if len(benchmark_result.preprocess) > 0:
        preprocess_results = sorted(benchmark_result.preprocess.values(), key=lambda r: r.utc)
        return preprocess_results[0].utc
    else:
        solver_results = sorted(benchmark_result.solver.values(), key=lambda r: r.utc)
        return solver_results[0].utc
