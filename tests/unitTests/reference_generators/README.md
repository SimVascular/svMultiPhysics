# Unit-test reference generators

These tools generate the trusted trajectory CSVs stored in
`tests/unitTests/reference_data`. They use the cited external model sources;
none reads svMultiPhysics output. The separately produced Bueno--Orovio and
Regazzoni references are documented in their model directories but are not
regenerated in-repo.

Generators are grouped by tested interface:

```text
active_stress/   Nash-Panfilov generator and Regazzoni reference documentation
ionic_model/     Aliev-Panfilov, FitzHugh-Nagumo, TP06, and BO documentation
```

Generator READMEs record their source, protocol, adaptations, and generation
command. The Bueno--Orovio and Regazzoni READMEs instead record how their
external references were produced. Generators write only to an explicit
`--output` path or stdout and do not modify repository reference data.

## Requirements

The Python reference generators and verifier require Python 3.10+ and use only
the standard library.

For example:

```bash
python3 ionic_model/aliev_panfilov/generate_aliev_panfilov.py \
  --output /tmp/ionic_aliev_panfilov_stimulated_trajectory.csv
```

Verify all in-repository generated references against an svMultiPhysics
checkout with:

```bash
python3 verify_reference_data.py --repo /path/to/svMultiPhysics
```

The verifier runs each generator twice to confirm deterministic output, then
compares the result byte-for-byte against the committed canonical CSV.
Bueno--Orovio and Regazzoni are reported as NOT CHECKED (intentional — no
in-repository generator exists for them). The verifier never overwrites
canonical files.
