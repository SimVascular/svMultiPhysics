
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
import subprocess
from typing import Self, assert_never
from .execution import Command, RunResult, capture_command
from ._expect import expect_float, expect_float_positive, expect_str

@dataclass(frozen=True, kw_only=True)
class EvaluationCriterion:
    """
    Datamodel of a single Evaluation Criteria
    """
    reference : float
    name : str
    command : Command
    tolerance : AbsoluteTolerance | RelativeTolerance
    cwd : Path

    @classmethod
    def from_dict(cls, name : str,  data : dict[str, object]) -> Self:
        reference = expect_float(data['reference'], 'reference')
        command = Command.from_dict(data)
        tolerance = cls.resolve_tolerance_from_dict(data)
        cwd = Path(expect_str(data.get('cwd', ''), 'cwd'))

        return cls(reference = reference, 
                   name = name, 
                   command = command, 
                   tolerance = tolerance,
                   cwd = cwd)

    @staticmethod
    def resolve_tolerance_from_dict(data : dict[str, object]) -> AbsoluteTolerance | RelativeTolerance:
        if 'atol' in data:
            atol = expect_float_positive(data['atol'], 'atol')
        else:
            atol = None

        if 'rtol' in data:
            rtol = expect_float_positive(data['rtol'], 'rtol')
        else:
            rtol = None

        # catch each case of atol, rtol being specified or absent. Support only atol xor rtol
        match (atol, rtol):
            case (value, None) if value is not None:
                tolerance = AbsoluteTolerance(value=value)

            case (None, value) if value is not None:
                tolerance = RelativeTolerance(value=value)

            case (None, None):
                raise ValueError("Either 'atol' or 'rtol' must be specified")

            case (value1, value2):
                raise ValueError("'atol' and 'rtol' cannot both be specified")

            case unexpected:
                assert_never(unexpected)

        return tolerance

@dataclass(frozen=True)
class AbsoluteTolerance:
    value : float

@dataclass(frozen=True)
class RelativeTolerance:
    value : float

@dataclass(frozen=True)
class EvaluationCriterionResult:
    """
    Data model of a comparison between a criterion against its measured value
    """

    reference : float
    tolerance : RelativeTolerance | AbsoluteTolerance | None
    value : float
    status : EvaluationStatus
    run_result : RunResult

class EvaluationStatus(StrEnum):
        PASSED = 'pass'
        FAILED = 'fail'

def evaluate_criterion(criterion : EvaluationCriterion, cwd : Path) -> EvaluationCriterionResult :
    """
    This function evaluates a criterion on a benchmark "end-to-end" using the Command provided on the Benchmark,
    collecting the result, and then evaluating it against the declared criterion
    """

    try:
        # extract value
        capture_result = capture_command(criterion.command, cwd = cwd)
    except subprocess.CalledProcessError as e:
        e.add_note(f"Error during script run for criterion evaluation.{criterion.name}")
        raise
    
    value = float(capture_result.stdout.strip()) # cast as a float

    # handle the two tolerance types
    match criterion.tolerance:

        case AbsoluteTolerance(value=atol):
            if abs(value - criterion.reference) < atol:
                status = EvaluationStatus.PASSED
            else:
                status = EvaluationStatus.FAILED

        case RelativeTolerance(value=rtol):
            if abs(value - criterion.reference) < abs(rtol*criterion.reference):
                status = EvaluationStatus.PASSED
            else:
                status = EvaluationStatus.FAILED

        # handle undefined cases        
        case None:
            raise NotImplementedError("No defined tolerance (likely for future metric tracking) is not yet supported")
        case unexpected:
            assert_never(unexpected)

    return EvaluationCriterionResult(reference = criterion.reference,
                                     tolerance = criterion.tolerance,
                                     value = value,
                                     status = status,
                                     run_result=capture_result.run_result)