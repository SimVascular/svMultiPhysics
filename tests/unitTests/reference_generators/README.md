# Unit-test reference generators

These tools generate the trusted trajectory CSVs stored in
`tests/unitTests/reference_data`. They implement the cited model equations or
run a pinned authors' reference implementation; none reads svMultiPhysics
output.

Generators are grouped by tested interface:

```text
active_stress/   Nash-Panfilov and Regazzoni
ionic_model/     Aliev-Panfilov, FitzHugh-Nagumo, Bueno-Orovio, and TP06
```

Each model README records its source, protocol, caveats, and generation
command. Generators write only to an explicit `--output` path or stdout and do
not modify repository reference data.

## Requirements

The Python reference generators and verifier require Python 3.10+ and use only
the standard library.

The Regazzoni generator additionally requires Git, `patch`, a C++17 compiler,
Eigen, Boost `property_tree`, and a local clone of the authors'
`cardiac-activation` repository containing the pinned reference commit.

For example:

```bash
python3 ionic_model/bueno_orovio/generate_bueno_orovio.py \
  --profile epi --output /tmp/ionic_bueno_orovio_epi_trajectory.csv
```

Verify all generated references against an svMultiPhysics checkout with:

```bash
python3 verify_reference_data.py --repo /path/to/svMultiPhysics \
  --regazzoni-reference-repo /path/to/cardiac-activation
```

Use `--skip-regazzoni` if its external C++ dependencies are unavailable. The
verifier uses byte comparison for Python-generated references and strict
numerical comparison for Regazzoni; it never overwrites canonical files.
