
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Mapping, Self

from .manifest import Benchmark, BenchmarkRunResult
from .execution import SolveExecutionContext, NamedCommand, RunResult, run_command
from .evaluation import EvaluationCriterionResult, evaluate_criterion



def run_benchmark(benchmark : Benchmark, context : SolveExecutionContext) -> BenchmarkRunResult:

    """
    Runs all the defined steps in the benchmark
    """

    preprocess_results = run_preprocess(benchmark)
    solver_results = run_solver(benchmark, context)
    postprocess_results = run_postprocess(benchmark)
    evaluation_results = run_evaluation(benchmark)
    cleanup_results = run_cleanup(benchmark)

    return BenchmarkRunResult.from_results(
                              benchmark,
                              preprocess_results,
                              solver_results,
                              postprocess_results, 
                              evaluation_results, 
                              cleanup_results)

def _run_named_commands(named_commands : tuple[NamedCommand,...], cwd : Path) -> Mapping[str, RunResult]:
    """
    Helper function that runs through the named commands, and collates them into a dict[name : str, result : RunResult]
    """
    return_dict : dict[str, RunResult] = {}
    for _named_command in named_commands:
        return_dict[_named_command.name] = run_command(_named_command.command, cwd)
    return return_dict

def run_preprocess(benchmark : Benchmark) -> Mapping[str, RunResult]:
    """
    Runs the preprocessing definitions, returning a mapping between the named preprocess and its run result
    """
    return _run_named_commands(benchmark.preprocess, benchmark.base_directory)

def run_solver(benchmark : Benchmark, context : SolveExecutionContext) -> Mapping[str, RunResult]:
    """
    Runs the solver definitions, returning a mapping between the defined solve and its result
    
    """
    # collapse the solver tuple into a named_command
    named_commands : tuple[NamedCommand,...] = tuple(r.build_named_command(context) for r in benchmark.solver)
    return _run_named_commands(named_commands, benchmark.base_directory)

def run_postprocess(benchmark : Benchmark) -> Mapping[str, RunResult]:
    """
    Runs the postprocessing definitions, returning a mapping between the named preprocess and its run result
    """
    return _run_named_commands(benchmark.postprocess, benchmark.base_directory)

def run_evaluation(benchmark : Benchmark) -> Mapping[str, EvaluationCriterionResult]:
    """
    Runs the evaluation definitions, returning a mapping between the named preprocess and its run result
    """
    return_dict : dict[str, EvaluationCriterionResult] = {}
    for _criterion in benchmark.evaluation:
        return_dict[_criterion.name] = evaluate_criterion(_criterion, cwd=benchmark.base_directory)
    return return_dict

def run_cleanup(benchmark : Benchmark) -> Mapping[str, RunResult]:
    """
    Runs the cleanup definitions, returning a mapping between the named preprocess and its run result
    """
    return _run_named_commands(benchmark.cleanup, benchmark.base_directory)