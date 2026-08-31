# Unit-test reference generators

These tools generate the trusted trajectory CSVs stored in
`tests/unitTests/reference_data`. They implement the cited model equations;
none reads svMultiPhysics output. The separately produced Regazzoni reference
is documented under `active_stress/regazzoni` but is not regenerated in-repo.

Generators are grouped by tested interface:

```text
active_stress/   Nash-Panfilov generator and Regazzoni reference documentation
ionic_model/     Aliev-Panfilov, FitzHugh-Nagumo, Bueno-Orovio, and TP06
```

Generator READMEs record their source, protocol, caveats, and generation
command. The Regazzoni README instead records how its external reference was
produced. Generators write only to an explicit `--output` path or stdout and
do not modify repository reference data.

## Requirements

The Python reference generators and verifier require Python 3.10+ and use only
the standard library.

For example:

```bash
python3 ionic_model/bueno_orovio/generate_bueno_orovio.py \
  --profile epi --output /tmp/ionic_bueno_orovio_epi_trajectory.csv
```

Verify all generated references against an svMultiPhysics checkout with:

```bash
python3 verify_reference_data.py --repo /path/to/svMultiPhysics
```

The verifier checks the nine in-repository generators byte-for-byte and reports
that Regazzoni is not checked. It never overwrites canonical files.
