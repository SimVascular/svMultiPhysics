# Regazzoni active-stress reference

This generator runs Francesco Regazzoni's MIT-licensed
[`cardiac-activation`](https://github.com/FrancescoRegazzoni/cardiac-activation)
implementation at commit
`26f05df28891df7b3c69f16bb136cdced6b63c4d`, using the RDQ20-MF human
body-temperature parameter set. It is not a Python transcription.

The pinned source leaves six structural-zero entries of the cross-bridge
matrix uninitialized. `zero_initialize_xb_matrix.patch` adds
`XB_A.setZero()` before the populated entries are assigned. The patch is
applied only to a temporary source copy.

The canonical protocol uses 600 outer updates with `dt=1 ms`, a normalized
double-exponential calcium transient, and a raised-cosine sarcomere-length
shortening/recovery cycle. The authors' code advances regulatory-unit states
with Forward-Euler substeps and cross-bridge moments with implicit Euler.
Tension is converted from kPa to MPa. The output schema is
`step,Ta,s0,...,s19` at checkpoints `0,30,99,157,299,599`.

Requirements are Python 3.10+, Git, `patch`, a C++17 compiler, Eigen headers,
Boost `property_tree`, and a checkout containing the pinned commit. Eigen and
Boost paths can be supplied explicitly or through `EIGEN3_INCLUDE_DIR` and
`BOOST_INCLUDE_DIR`.

```bash
python3 generate_regazzoni.py \
  --reference-repo /path/to/cardiac-activation \
  --output active_stress_regazzoni_twitch.csv
```

Compiler and math-library differences may affect final printed bits, so the
package verifier compares this reference numerically with strict tolerances.
