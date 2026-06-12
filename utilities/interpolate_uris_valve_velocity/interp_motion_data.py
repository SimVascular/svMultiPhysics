#!/usr/bin/env python3
"""Interpolate selected frame ranges from a valve-motion .dat file.

The input .dat file is expected to have this format:

    <n_steps> <n_nodes>
    x y z
    x y z
    ...

Coordinate rows are time-major: all nodes for frame 0, then all nodes for
frame 1, and so on.

By default, edit the values in ``DEFAULT_CONFIG`` and run:

    python interp_motion_data_clean.py

You can also override the input/output paths from the command line:

    python interp_motion_data.py --input NCC_motion.dat --output NCC_motion_interp.dat

Notes
-----
* Each segment request is ``(start_frame, end_frame, num_samples)``.
* ``num_samples`` includes both endpoints.
* If ``start_frame <= end_frame``, samples are taken along the original time
  sequence using the selected interpolation method.
* If ``start_frame > end_frame``, this linearly blends directly between the two keyframes.
  It does *not* wrap through the end of the sequence.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

COMMENT_PREFIXES = ("#", "//", "%", ";")
VALID_METHODS = {"linear", "quadratic"}


@dataclass(frozen=True)
class SegmentRequest:
    """Request ``num_samples`` frames between two source frames, inclusive."""

    start_frame: int
    end_frame: int
    num_samples: int
    use_two_frame_blend: bool = False


@dataclass(frozen=True)
class InterpolationConfig:
    """User-editable interpolation settings."""

    input_dat: Path
    output_dat: Path
    method: str = "linear"
    requests: tuple[SegmentRequest, ...] = ()
    join_touching_endpoints: bool = True
    precision: int = 6


# ---------------------------------------------------------------------------
# Default settings
# ---------------------------------------------------------------------------
# Edit this block for motion interpolation. Command-line arguments can override the
# input/output paths and interpolation method when needed.
DEFAULT_CONFIG = InterpolationConfig(
    input_dat=Path("NCC_motion.dat"),
    output_dat=Path("NCC_motion_interp.dat"),
    method="linear",
    requests=(
        SegmentRequest(start_frame=290, end_frame=329, num_samples=40),
    ),
    join_touching_endpoints=True,
    precision=6,
)

# Example close-cycle configuration. Copy into DEFAULT_CONFIG.requests when
# needed instead of keeping many commented-out blocks throughout the script.
EXAMPLE_CLOSE_REQUESTS = (
    SegmentRequest(329, 420, 10, use_two_frame_blend=True),
    SegmentRequest(420, 450, 10, use_two_frame_blend=True),
    SegmentRequest(450, 290, 22, use_two_frame_blend=False),
)


# ---------------------------------------------------------------------------
# File I/O
# ---------------------------------------------------------------------------
def read_dat(dat_path: str | Path) -> np.ndarray:
    """Read a motion .dat file as an array with shape ``(n_steps, n_nodes, 3)``."""

    dat_path = Path(dat_path)
    with dat_path.open("r", encoding="utf-8") as handle:
        header = _read_header(handle)
        n_steps, n_nodes = _parse_header(header)
        coords = _read_coordinates(handle)

    expected_rows = n_steps * n_nodes
    if coords.shape[0] != expected_rows:
        raise ValueError(
            f"Data length mismatch in {dat_path}: header says "
            f"{n_steps} * {n_nodes} = {expected_rows} coordinate rows, "
            f"but found {coords.shape[0]}."
        )

    return coords.reshape(n_steps, n_nodes, 3)


def _read_header(handle) -> str:
    """Return the first non-empty, non-comment line from an open file handle."""

    for raw_line in handle:
        line = raw_line.strip()
        if line and not line.startswith(COMMENT_PREFIXES):
            return line
    raise ValueError("Input file is empty or missing a '<n_steps> <n_nodes>' header.")


def _parse_header(header: str) -> tuple[int, int]:
    """Parse the ``<n_steps> <n_nodes>`` header line."""

    parts = header.replace(",", " ").split()
    if len(parts) < 2:
        raise ValueError("Header must contain '<n_steps> <n_nodes>'.")

    try:
        n_steps = int(float(parts[0]))
        n_nodes = int(float(parts[1]))
    except ValueError as exc:
        raise ValueError(f"Could not parse .dat header: {header!r}") from exc

    if n_steps <= 0 or n_nodes <= 0:
        raise ValueError("Header values n_steps and n_nodes must be positive.")

    return n_steps, n_nodes


def _read_coordinates(handle) -> np.ndarray:
    """Read all coordinate rows from an open .dat file handle."""

    rows: list[list[float]] = []
    for line_number, raw_line in enumerate(handle, start=2):
        line = raw_line.strip()
        if not line or line.startswith(COMMENT_PREFIXES):
            continue

        parts = line.replace(",", " ").split()
        if len(parts) != 3:
            raise ValueError(
                f"Expected three coordinate values on line {line_number}, got: {line!r}"
            )

        try:
            rows.append([float(value) for value in parts])
        except ValueError as exc:
            raise ValueError(
                f"Could not parse coordinate values on line {line_number}: {line!r}"
            ) from exc

    return np.asarray(rows, dtype=np.float64)


def write_dat(
    out_path: str | Path,
    frames: Sequence[np.ndarray],
    n_nodes: int,
    precision: int = 6,
) -> None:
    """Write interpolated frames to a motion .dat file."""

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fmt = f"%.{precision}f %.{precision}f %.{precision}f"
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write(f"{len(frames)} {n_nodes}\n")
        for frame in frames:
            np.savetxt(handle, frame, fmt=fmt)


# ---------------------------------------------------------------------------
# Interpolation
# ---------------------------------------------------------------------------
def interpolate_frames(coords: np.ndarray, config: InterpolationConfig) -> list[np.ndarray]:
    """Build all requested interpolated frames from ``coords``."""

    method = config.method.lower().strip()
    if method not in VALID_METHODS:
        raise ValueError(f"method must be one of {sorted(VALID_METHODS)}, got {config.method!r}.")

    if not config.requests:
        raise ValueError("At least one SegmentRequest is required.")

    frames_out: list[np.ndarray] = []
    previous_end: int | None = None

    for request in config.requests:
        _validate_request(request, n_steps=coords.shape[0])
        segment = sample_segment(coords, request, method)

        if (
            config.join_touching_endpoints
            and previous_end is not None
            and request.start_frame == previous_end
        ):
            segment = segment[1:]

        frames_out.extend(segment)
        previous_end = request.end_frame

    if not frames_out:
        raise RuntimeError("No frames produced. Check the segment requests.")

    return frames_out


def _validate_request(request: SegmentRequest, n_steps: int) -> None:
    """Validate a segment request against the source sequence length."""

    if not 0 <= request.start_frame < n_steps:
        raise ValueError(f"start_frame must be in [0, {n_steps - 1}], got {request.start_frame}.")
    if not 0 <= request.end_frame < n_steps:
        raise ValueError(f"end_frame must be in [0, {n_steps - 1}], got {request.end_frame}.")
    if request.num_samples < 2:
        raise ValueError("num_samples must be at least 2 because endpoints are included.")


def sample_segment(
    coords: np.ndarray,
    request: SegmentRequest,
    method: str,
) -> list[np.ndarray]:
    """Sample one configured segment."""

    if request.use_two_frame_blend or request.start_frame > request.end_frame:
        return sample_two_keyframe_blend(
            coords,
            request.start_frame,
            request.end_frame,
            request.num_samples,
        )

    return sample_along_sequence(
        coords,
        request.start_frame,
        request.end_frame,
        request.num_samples,
        method,
    )


def sample_along_sequence(
    coords: np.ndarray,
    start_frame: int,
    end_frame: int,
    num_samples: int,
    method: str,
) -> list[np.ndarray]:
    """Sample a non-wrapping segment along the original time sequence."""

    sample_times = np.linspace(float(start_frame), float(end_frame), num=num_samples)
    return [interpolate_at_time(coords, time, method) for time in sample_times]


def sample_two_keyframe_blend(
    coords: np.ndarray,
    start_frame: int,
    end_frame: int,
    num_samples: int,
) -> list[np.ndarray]:
    """Linearly blend between two keyframes, ignoring intermediate frames."""

    start_coords = coords[start_frame]
    end_coords = coords[end_frame]
    alphas = np.linspace(0.0, 1.0, num=num_samples)
    return [linear_blend(start_coords, end_coords, alpha) for alpha in alphas]


def interpolate_at_time(coords: np.ndarray, time: float, method: str) -> np.ndarray:
    """Interpolate the source sequence at a real-valued frame index."""

    rounded = round(time)
    if abs(time - rounded) < 1e-12:
        return coords[int(rounded)]

    if method == "linear":
        return interpolate_linear_sequence(coords, time)
    if method == "quadratic":
        return interpolate_quadratic_sequence(coords, time)

    raise ValueError(f"Unknown interpolation method: {method!r}")


def linear_blend(start_coords: np.ndarray, end_coords: np.ndarray, alpha: float) -> np.ndarray:
    """Linearly blend two coordinate arrays at fraction ``alpha`` in ``[0, 1]``."""

    return (1.0 - alpha) * start_coords + alpha * end_coords


def interpolate_linear_sequence(coords: np.ndarray, time: float) -> np.ndarray:
    """Piecewise-linear interpolation along adjacent source frames."""

    n_steps = coords.shape[0]
    if time <= 0:
        return coords[0]
    if time >= n_steps - 1:
        return coords[-1]

    left = int(np.floor(time))
    alpha = time - left
    return linear_blend(coords[left], coords[left + 1], alpha)


def interpolate_quadratic_sequence(coords: np.ndarray, time: float) -> np.ndarray:
    """Local quadratic interpolation through frames ``t-1``, ``t``, and ``t+1``.

    The function falls back to linear interpolation at the sequence boundaries,
    where a three-frame stencil is unavailable.
    """

    n_steps = coords.shape[0]
    if time <= 0:
        return coords[0]
    if time >= n_steps - 1:
        return coords[-1]

    center = int(np.floor(time))
    alpha = time - center

    if center - 1 < 0 or center + 1 >= n_steps:
        return interpolate_linear_sequence(coords, time)

    previous_frame = coords[center - 1]
    current_frame = coords[center]
    next_frame = coords[center + 1]

    # Lagrange basis evaluated at alpha for nodes -1, 0, +1.
    weight_previous = alpha * (alpha - 1.0) / 2.0
    weight_current = 1.0 - alpha * alpha
    weight_next = alpha * (alpha + 1.0) / 2.0

    return (
        weight_previous * previous_frame
        + weight_current * current_frame
        + weight_next * next_frame
    )


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    """Parse optional command-line overrides."""

    parser = argparse.ArgumentParser(
        description="Interpolate selected frame ranges from a valve-motion .dat file."
    )
    parser.add_argument("--input", type=Path, help="Input .dat file.")
    parser.add_argument("--output", type=Path, help="Output .dat file.")
    parser.add_argument(
        "--method",
        choices=sorted(VALID_METHODS),
        help="Interpolation method for along-sequence segments.",
    )
    parser.add_argument(
        "--precision",
        type=int,
        help="Number of decimal places written to the output file.",
    )
    return parser.parse_args()


def config_from_args(args: argparse.Namespace) -> InterpolationConfig:
    """Apply command-line overrides to ``DEFAULT_CONFIG``."""

    return InterpolationConfig(
        input_dat=args.input or DEFAULT_CONFIG.input_dat,
        output_dat=args.output or DEFAULT_CONFIG.output_dat,
        method=args.method or DEFAULT_CONFIG.method,
        requests=DEFAULT_CONFIG.requests,
        join_touching_endpoints=DEFAULT_CONFIG.join_touching_endpoints,
        precision=args.precision if args.precision is not None else DEFAULT_CONFIG.precision,
    )


def main() -> None:
    """Read the input file, interpolate requested segments, and write the output."""

    config = config_from_args(parse_args())
    coords = read_dat(config.input_dat)
    frames = interpolate_frames(coords, config)
    write_dat(config.output_dat, frames, n_nodes=coords.shape[1], precision=config.precision)
    print(f"Wrote {len(frames)} frames to {config.output_dat}")


if __name__ == "__main__":
    main()
