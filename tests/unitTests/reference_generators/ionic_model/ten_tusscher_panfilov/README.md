# Ten Tusscher--Panfilov ionic references

The three model cores are adapted from CellML API-generated Python for the
Physiome TP06 [EPI][epi-cellml], [ENDO][endo-cellml], and [M-cell][m-cellml]
models. All are in exposure `a7179d94365ff0c9c0e6eb7c6a787d3d`, from
workspace `tentusscher_panfilov_2006` at pinned revision
[`5dc42395eef6044fe766786f7bff197dea355eb3`][revision]. They are separate
phenotype artifacts in the same revision:

| Phenotype | CellML SHA-256 | Generated Python | Generated-Python SHA-256 |
| --- | --- | --- | --- |
| EPI | `999e4af049776f18a35aa30e09c3a8fa9487ab0b2e9ef6a77489dd371e1d6e6a` | [source][epi-python] | `1ffdfbdd36df8804e0580ca8dc277579c24ff7f104bc89108a9e031a12f94085` |
| ENDO | `b8c49166fa32aa14f6bbde21527fdd9ce6fb29f086a438d56aa479372f70467d` | [source][endo-python] | `cb42893460dc6d96e7aa79d3a8ab935ac43c1979976904a59da1ef7b78fefb0a` |
| M | `f97df16ce337f031dd6019478dfe5393d9e5d1c517ba386016f138ea3ec46074` | [source][m-python] | `af1dbdf3a08eff7611e251cc6587bd6674fda66538ec1b49f05d974e377b16e6` |

The [EPI][epi-metadata], [ENDO][endo-metadata], and [M-cell][m-metadata]
metadata name Penny Noble as the CellML model author and cite ten Tusscher and
Panfilov (2006). The [EPI][epi-license], [ENDO][endo-license], and
[M-cell][m-license] license/citation pages state CC BY 3.0 Unported
([terms][license]); attribution and adaptation notices are retained in each
model core.

## Adaptation and protocol

Each core retains its generated initial values, constants, and `computeRates`
expressions, except that the built-in periodic stimulus is replaced by a
caller-supplied value. The generated algebraic workspace is returned so the
thin driver can use the source steady states and time constants for
Rush--Larsen gate updates. Unused legends, vectorized algebraic output,
NumPy/SciPy helpers, demonstration solver, and plotting are removed; scalar
NumPy operations are replaced by equivalent standard-library operations.

The driver advances generated states `V, Ki, Nai, Cai, Ca_ss, Ca_SR, R_prime`
with Forward Euler and the twelve gates with Rush--Larsen, all from the old
state. It maps the generated CellML ordering to CSV columns containing those
seven states followed by `Xr1, Xr2, Xs, m, h, j, d, f, f2, fCass, s, r`.
EPI, ENDO, and M retain their source initial conditions, conductances, and
`s`-gate kinetics. The source `ICaL` expression is retained without
regularizing its removable singularity at exactly `V=15 mV`.

All profiles use `dt=0.005 ms` for 120000 updates, with
`Istim=-52 pA/pF` for `10 <= t < 11 ms`, zero stimulus otherwise, zero SAC,
and no pre-pacing. Stimulus is evaluated at old-state time. Output contains
only the CSV header and selected checkpoint rows.

```bash
python3 generate_ttp.py --profile epi --output ionic_ttp_epi_trajectory.csv
```

Valid profiles are `epi`, `endo`, and `m`.

[revision]: https://models.physiomeproject.org/workspace/tentusscher_panfilov_2006/@@file/5dc42395eef6044fe766786f7bff197dea355eb3/
[epi-cellml]: https://models.physiomeproject.org/workspace/tentusscher_panfilov_2006/rawfile/5dc42395eef6044fe766786f7bff197dea355eb3/ten_tusscher_model_2006_IK1Ko_epi_units.cellml
[endo-cellml]: https://models.physiomeproject.org/workspace/tentusscher_panfilov_2006/rawfile/5dc42395eef6044fe766786f7bff197dea355eb3/ten_tusscher_model_2006_IK1Ko_endo_units.cellml
[m-cellml]: https://models.physiomeproject.org/workspace/tentusscher_panfilov_2006/rawfile/5dc42395eef6044fe766786f7bff197dea355eb3/ten_tusscher_model_2006_IK1Ko_M_units.cellml
[epi-python]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_epi_units.cellml/@@cellml_codegen/Python/raw
[endo-python]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_endo_units.cellml/@@cellml_codegen/Python/raw
[m-python]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_M_units.cellml/@@cellml_codegen/Python/raw
[epi-metadata]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_epi_units.cellml/cmeta
[endo-metadata]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_endo_units.cellml/cmeta
[m-metadata]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_M_units.cellml/cmeta
[epi-license]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_epi_units.cellml/license_citation
[endo-license]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_endo_units.cellml/license_citation
[m-license]: https://models.physiomeproject.org/exposure/a7179d94365ff0c9c0e6eb7c6a787d3d/ten_tusscher_model_2006_IK1Ko_M_units.cellml/license_citation
[license]: https://creativecommons.org/licenses/by/3.0/legalcode
