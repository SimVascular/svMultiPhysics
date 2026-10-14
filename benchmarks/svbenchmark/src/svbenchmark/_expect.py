"""
Basic runtime checks to ensure correct types and help type hinting.

Follows semantics of Rust "expect(...)" to raise an error if the type is not correct.

"""

from typing import cast
from math import isfinite
from collections.abc import Callable


def expect_int(value: object, location: str) -> int:
    """
    Checks if value is int, returning it if so. Otherwise throws an error
    """
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(
            f"Expected '{location}' to be int, got {type(value).__name__}"
        )
    return value


def expect_float(value: object, location: str) -> float:
    """
    Checks if value is float, returning it if so. Otherwise throws an error
    """
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(
            f"Expected '{location}' to be numeric, got {type(value).__name__}"
        )
    return float(value)

def expect_float_positive(value: object, location: str) -> float:
    """
    Checks if value is positive float, returning it if so. Otherwise throws an error
    """
    result = expect_float(value, location)

    if not isfinite(result) or result <= 0.0:
        raise ValueError(
            f"Expected '{location}' to be a finite positive number, got {result!r}"
        )

    return result

def expect_str(value: object, location: str) -> str:
    """
    Checks if value is str, returning it if so. Otherwise throws an error
    """
    if not isinstance(value, str):
        raise ValueError(
            f"Expected '{location}' to be str, got {type(value).__name__}"
        )
    return value


def expect_table(value: object, location: str) -> dict[str, object]:
    """
    Checks if value is table (dict), returning it if so. Otherwise throws an error
    """
    if not isinstance(value, dict):
        raise ValueError(
            f"Expected '{location}' to be a table, got {type(value).__name__}"
        )

    if not all(isinstance(key, str) for key in value):
        raise ValueError(
            f"Expected '{location}' to have string keys"
        )

    return cast(dict[str, object], value)


def expect_str_array(value: object, location: str) -> tuple[str, ...]:
    """
    Checks if value is [str,...], returning it as a tuple if so. Otherwise throws an error
    """
    if not isinstance(value, list) or not all(
        isinstance(item, str) for item in value
    ):
        raise ValueError(
            f"Expected '{location}' to be an array of strings"
        )

    return tuple(value)