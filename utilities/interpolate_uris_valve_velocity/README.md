# Interpolate URIS Valve Velocity

This utility generates an interpolated valve motion dataset from an existing valve motion `.dat` file. The primary use case is to increase the temporal resolution of prescribed valve motion before computing the valve velocity used in the URIS forcing term.

## Background

The URIS formulation includes the prescribed valve velocity $u_{\Gamma_i}$,

$$
f_{\mathrm{URIS}} = \sum_{i=1}^{n} \left( u - u_{\Gamma_i} \right) \text{ .}
$$

The valve velocity is computed using finite differences of the prescribed valve motion. If the input motion is sampled with relatively large time steps, the resulting velocity can be noisy and may introduce interpolation errors. This utility creates additional intermediate frames to obtain a smoother estimate of the valve velocity.

## Input Format

The script expects a valve motion `.dat` file with the following format:

```text
<n_steps> <n_nodes>
x y z
x y z
...
```

where:

* `n_steps` is the number of time frames.
* `n_nodes` is the number of mesh nodes.
* Coordinates are stored in time-major order:

  * all nodes for frame 0,
  * followed by all nodes for frame 1,
  * etc.

For example:

```text
40 200
0.123 1.234 2.345
0.456 1.567 2.678
...
```

## Configuration

The script is configured through the `REQUESTS` list.

Each segment is specified as:

```python
Segment(start_step, end_step, num_output_frames, two_frame=False)
```

Parameters:

* `start_step`: first frame of the segment.
* `end_step`: last frame of the segment.
* `num_output_frames`: number of output frames to generate, including both endpoints.
* `two_frame`:

  * `False`: interpolate using the original motion trajectory.
  * `True`: linearly blend between the start and end frames only.

Example:

```python
REQUESTS = [
    Segment(290, 329, 40),
]
```

This generates 40 frames between timesteps 290 and 329, including both endpoints.

## Interpolation Methods

Two interpolation methods are available:

```python
METHOD = "linear"
```

or

```python
METHOD = "quadratic"
```

### Linear

Piecewise-linear interpolation between neighboring frames.

### Quadratic

Local quadratic interpolation using three neighboring frames. Near boundaries, the method automatically falls back to linear interpolation.

## Example

Input file:

```text
NCC_motion.dat
```

Configuration:

```python
leaflet_name = "NCC"
state = "open"

METHOD = "linear"

REQUESTS = [
    Segment(290, 329, 40),
]
```

Run:

```bash
python interp_motion_data.py
```

Output:

```text
NCC_motion_interp.dat
```

## Multi-Segment Example

```python
REQUESTS = [
    Segment(329, 420, 10, two_frame=True),
    Segment(420, 450, 10, two_frame=True),
    Segment(450, 290, 22),
]
```

The generated segments are concatenated into a single output motion file. If consecutive segments share an endpoint, duplicate frames are automatically removed.

## Notes

* All output segments include both endpoints.
* The output file preserves the original mesh topology.
* Only nodal coordinates are interpolated.
* The script is intended as a preprocessing utility for generating smoother valve velocity data used by URIS simulations.
