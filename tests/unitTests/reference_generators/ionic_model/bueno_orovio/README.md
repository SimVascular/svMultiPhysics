# Bueno--Orovio ionic references

The committed EPI, ENDO, and M-cell references were produced from openCARP's
external [`Bueno.model`](https://git.opencarp.org/openCARP/openCARP/-/blob/ac4e96792db082958fa9a830341f6e0642cc6bb8/physics/limpet/models/Bueno.model)
implementation at commit
`ac4e96792db082958fa9a830341f6e0642cc6bb8`. The model file has SHA-256
`71a70bc3771fc4dce98af1a88051df1f732d0c4dc079be44e17cb7f4c25734bf`;
the generated CPU implementation `physics/limpet/src/imps_src/Bueno.cc` at
that commit has SHA-256
`e0b4886910cb928bc08e0314aa54f22a485e1c69130cc272cef23a5690fd436c`.
The model cites Bueno-Orovio, Cherry, and Fenton (2008),
[doi:10.1016/j.jtbi.2008.03.029](https://doi.org/10.1016/j.jtbi.2008.03.029).

The pinned openCARP source is covered by the
[openCARP Academic Public License v1.1](https://git.opencarp.org/openCARP/openCARP/-/blob/ac4e96792db082958fa9a830341f6e0642cc6bb8/LICENSE.md).
`Bueno.model` has no file-level license override. No openCARP source is copied
into svMultiPhysics.

## Parameters and protocol

The openCARP model defaults are the ENDO parameters. EPI and M-cell runs set
the following phenotype values explicitly:

| Parameter | EPI | ENDO | M |
|---|---:|---:|---:|
| `u_u` | 1.55 | 1.56 | 1.61 |
| `theta_v_minus` | 0.006 | 0.20 | 0.10 |
| `theta_o` | 0.006 | 0.006 | 0.005 |
| `tau_v1_minus` | 60 | 75 | 80 |
| `tau_v2_minus` | 1150 | 10 | 1.4506 |
| `tau_w1_minus` | 60 | 6 | 70 |
| `tau_w2_minus` | 15 | 140 | 8 |
| `k_w_minus` | 65 | 200 | 200 |
| `u_w_minus` | 0.03 | 0.016 | 0.016 |
| `tau_w_plus` | 200 | 280 | 280 |
| `tau_fi` | 0.11 | 0.10 | 0.078 |
| `tau_o1` | 400 | 470 | 410 |
| `tau_o2` | 6 | 6 | 7 |
| `tau_so1` | 30.0181 | 40 | 91 |
| `tau_so2` | 0.9957 | 1.2 | 0.8 |
| `k_so` | 2.0458 | 2 | 2.1 |
| `u_so` | 0.65 | 0.65 | 0.6 |
| `tau_s2` | 16 | 2 | **2** |
| `tau_si` | 1.8875 | 2.9013 | 3.3849 |
| `tau_w_inf` | 0.07 | 0.0273 | 0.01 |
| `w_inf_star` | 0.94 | 0.78 | 0.5 |

All profiles use `u_o=0`, `theta_v=0.3`, `theta_w=0.13`,
`tau_v_plus=1.4506`, `tau_s1=2.7342`, `k_s=2.0994`, `u_s=0.9087`,
`Cm=1`, `V_0=-84 mV`, and `V_mu=85.7 mV`. The five openCARP internal
`theta_w` constants are all `0.13`. The transformed formulation is used, with
`V_mV=V_0+V_mu*u` and initial state `(V_mV,v,w,s)=(-84,1,1,0)`.

Bueno-Orovio et al. Table 1 gives M-cell `tau_s2=4 ms`; svMultiPhysics
currently uses `2 ms`. The canonical M run explicitly uses `2 ms` to test the
shipped parameter set.

Each run applies `Istim=-35.714 pA/pF` for `10 <= t < 12 ms`, evaluated at
the old-state time, and zero otherwise. `Ksac=0`. All four states use
simultaneous Forward Euler with `dt=0.01 ms`; EPI runs for 600 ms (60,000
updates), while ENDO and M run for 1200 ms (120,000 updates). Output columns
are `step,V_mV,v,w,s`, and checkpoint `N` is the state after `N` completed
updates.

EPI checkpoints are `0, 1000, 1001, 1109, 1198, 1200, 2083, 6018, 8589,
10000, 20000, 23059, 27518, 28294, 28403, 30247, 40000, 60000`.
ENDO checkpoints are `0, 1000, 1001, 1107, 1200, 1972, 2500, 4944, 10000,
20000, 24067, 28075, 28557, 28755, 28876, 30720, 40000, 60000, 120000`.
M-cell checkpoints are `0, 1000, 1001, 1106, 1200, 2021, 2500, 4651,
10000, 20000, 30000, 36369, 44724, 46147, 46456, 46640, 48735, 60000,
120000`.

## Reproduction

Check out the pinned openCARP commit and build its generated CPU ionic-model
implementation. Configure `Bueno` in transformed mode with the parameter set
above. Because voltage is an external openCARP state, drive one cell by first
evaluating the generated `Iion` and `diff_v`, `diff_w`, and `diff_s` expressions
from state `n`, then apply the four Forward-Euler updates simultaneously. The
voltage update is

```text
V_mV[n+1] = V_mV[n] - dt * (Iion[n] + Istim[n]).
```

This driver ordering preserves the svMultiPhysics old-state stimulus and
simultaneous-update contract; openCARP's standard `bench` stimulus path applies
the voltage increment before evaluating the ionic model and is therefore not
used for these references. Compare regenerated checkpoint values numerically
with the committed CSVs using the existing trajectory-test tolerance; the
committed CSVs are the authoritative reference data.

The resulting files are:

- `tests/unitTests/reference_data/ionic_bueno_orovio_epi_trajectory.csv`
- `tests/unitTests/reference_data/ionic_bueno_orovio_endo_trajectory.csv`
- `tests/unitTests/reference_data/ionic_bueno_orovio_m_trajectory.csv`

No ready-to-run generator is included because the trusted implementation
remains external and depends on the openCARP build/runtime environment.
