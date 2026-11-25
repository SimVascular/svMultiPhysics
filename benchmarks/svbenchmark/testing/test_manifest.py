from pathlib import Path
import tomllib
import pytest
from svbenchmark.manifest import Benchmark


# use the clean_manifest toml to create a clean data_dict
@pytest.fixture
def clean_manifest():
    path = Path(__file__).parent / "fixtures" / "clean_manifest.toml"

    with path.open("rb") as f:
        return tomllib.load(f)


def test_from_toml(clean_manifest):
    Benchmark.from_toml_dict(Path(""), clean_manifest)


@pytest.mark.parametrize(
    "missing_key", ("schema", "name", "reference_source", "solver", "evaluation")
)
def test_missing_key(clean_manifest, missing_key):
    data_dict = clean_manifest
    data_dict.pop(missing_key)

    with pytest.raises((KeyError, ValueError)):
        Benchmark.from_toml_dict(Path(""), data_dict)


@pytest.mark.parametrize("optional_key", ("preprocess", "postprocess", "cleanup"))
def test_optional_phases(clean_manifest, optional_key):
    """
    Test the optional phases, they should still be populated, even if omitted
    """
    data_dict = clean_manifest
    data_dict.pop(optional_key)

    # should still be populated
    b = Benchmark.from_toml_dict(Path(""), data_dict)
    assert b.__dict__[optional_key] == tuple()


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("schema", "1."),
        ("name", 2),
        ("reference_source", 2),
        ("solver", "text"),
        ("evaluation", "text"),
    ],
)
def test_wrong_type_schema(clean_manifest, key, value):
    data_dict = clean_manifest
    clean_manifest[key] = value
    with pytest.raises(ValueError):
        Benchmark.from_toml_dict(Path(""), data_dict)
