

from dataclasses import dataclass
from pathlib import Path
import tomllib
from typing import Callable, Self

from .evaluation import EvaluationCriterion
from .execution import NamedCommand, SolverRun
from ._expect import expect_int, expect_str, expect_table

@dataclass(frozen=True, kw_only=True)
class Benchmark:
    """
    Holds the manifest of declarative Benchmark definition
    """
    name : str
    schema : int
    reference_source : str
    base_dir : Path
    preprocess : tuple[NamedCommand,...]  = tuple()
    solver : tuple[SolverRun,...]
    postprocess : tuple[NamedCommand,...] = tuple()
    evaluation : tuple[EvaluationCriterion,...]
    cleanup : tuple[NamedCommand,...]     = tuple()

    @classmethod
    def from_toml_dict(cls, base_dir : Path, data : dict[str, object]) -> Self:
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
        base_dir = base_dir.resolve()
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
                   base_dir = base_dir, 
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
   