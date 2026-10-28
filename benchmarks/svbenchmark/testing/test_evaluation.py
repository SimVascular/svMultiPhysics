import subprocess
from types import SimpleNamespace
from unittest.mock import Mock

import pytest

from svbenchmark.evaluation import EvaluationCriterion, EvaluationStatus, evaluate_criterion


class TestTolerance:

    @pytest.mark.parametrize('atol', [1e-4, 1.2e-8])
    def test_resolve_atol_from_dict(self, atol):
        """
        Normal creation of an absolute tolerance creation from dict
        """
        data_dict = {'atol' : atol}
        tolerance = EvaluationCriterion.resolve_tolerance_from_dict(data_dict)
        assert(atol == tolerance.value)

    @pytest.mark.parametrize('rtol', [0.99999, 0.00001])
    def test_resolve_rtol_from_dict(self, rtol):
        """
        Normal creation of a relative tolerance creation from dict
        """
        data_dict = {'rtol' : rtol}
        tolerance = EvaluationCriterion.resolve_tolerance_from_dict(data_dict)
        assert(rtol == tolerance.value)

    @pytest.mark.parametrize(['atol', 'rtol'], [(None, None), (1e-4, 1e-4), (2.e-4, 1.2e-3)])
    def test_resolve_tolerance_from_dict_invalid(self, atol, rtol):
        """
        All the invalid specifications (atol xor rtol must be specified)
        """
        data_dict = {}
        if atol is not None: data_dict['atol'] = atol
        if rtol is not None: data_dict['rtol'] = rtol

        with pytest.raises(ValueError):
            EvaluationCriterion.resolve_tolerance_from_dict(data_dict)

    @pytest.mark.parametrize(['atol', 'rtol'], [(-1e-4, None), (-1.2e-8,None), (-0.9999, -0.0001), (None, 0)])
    def test_resolve_tolerance_from_dict_neg_value(self, atol, rtol):
        """
        All tolerances must be positive
        """
        data_dict = {}
        if atol is not None: data_dict['atol'] = atol
        if rtol is not None: data_dict['rtol'] = rtol

        with pytest.raises(ValueError):
            EvaluationCriterion.resolve_tolerance_from_dict(data_dict)


class TestEvaluationCriterion:
    """
    Conditions:
    1. normal creation + evaluation (negative reference value)
        - rtol + atol
        - pass + fail
    2. normal creation + evaluation (negative reference value)
        - rtol + atol
        - pass + fail
    
    """

    @pytest.mark.parametrize(['atol', 'ref_value', 'value', 'expected'],
                            [(1e-4, 0, 1e-8, EvaluationStatus.PASSED),
                             (1e-4, 0, 1e-2, EvaluationStatus.FAILED),
                             (1e-2, -1, -1+1e-4, EvaluationStatus.PASSED),
                             (1e-2, 1, 1+1e-4, EvaluationStatus.PASSED),
                             (1e-8, 1, 1+1e-2, EvaluationStatus.FAILED),
                             (1e-8, -1, -1+1e-2, EvaluationStatus.FAILED)])
    def test_evaluate_criterion_atol(self, atol, ref_value, value, expected, monkeypatch):
        data_dict = {'atol' : atol,
                     'reference' : ref_value,
                     'argv' : [""]}

        # setup a mock capture for testing the return value only
        mock_capture = Mock(
            return_value=SimpleNamespace(stdout=f"{value}", run_result = None)
        )

        monkeypatch.setattr(
            'svbenchmark.evaluation.capture_command',
            mock_capture,
        )

        eval_criterion = EvaluationCriterion.from_dict("test", data_dict)
        eval_result = evaluate_criterion(eval_criterion, eval_criterion.cwd)

        assert(eval_result.status == expected)

    @pytest.mark.parametrize(
        ['rtol', 'ref_value', 'value', 'expected'],
        [
            (1e-4, 1, 1 + 1e-8, EvaluationStatus.PASSED),
            (1e-4, 1, 1 + 1e-2, EvaluationStatus.FAILED),
            (1e-2, -1, -1 + 1e-4, EvaluationStatus.PASSED),
            (1e-2, 1, 1 + 1e-4, EvaluationStatus.PASSED),
            (1e-8, 1, 1 + 1e-2, EvaluationStatus.FAILED),
            (1e-8, -1, -1 + 1e-2, EvaluationStatus.FAILED),
            (1e-2, 0, 0, EvaluationStatus.FAILED),
        ]
    )
    def test_evaluate_criterion_rtol(
        self, rtol, ref_value, value, expected, monkeypatch
    ):
        data_dict = {
            'rtol': rtol,
            'reference': ref_value,
            'argv': [""],
        }

        mock_capture = Mock(
            return_value=SimpleNamespace(stdout=f"{value}", run_result = None)
        )

        monkeypatch.setattr(
            'svbenchmark.evaluation.capture_command',
            mock_capture,
        )

        eval_criterion = EvaluationCriterion.from_dict("test", data_dict)
        eval_result = evaluate_criterion(eval_criterion, eval_criterion.cwd)

        assert eval_result.status == expected


    def test_evaluate_criterion_script_error(self, monkeypatch):
        mock_capture = Mock(
            side_effect=subprocess.CalledProcessError(
                returncode=1,
                cmd=["test-command"],
                stderr="Command failed",
            )
        )

        monkeypatch.setattr(
            'svbenchmark.evaluation.capture_command',
            mock_capture,
        )

        data_dict = {
            'atol': 1e-8,
            'reference': 0,
            'argv': [""],
        }

        eval_criterion = EvaluationCriterion.from_dict("test", data_dict)
        with pytest.raises(subprocess.CalledProcessError):
            evaluate_criterion(eval_criterion, eval_criterion.cwd)