# Ten Tusscher--Panfilov ionic references

The TP06 equations follow the curated Physiome CellML models for EPI, ENDO,
and M cells:

- workspace `604`
- revision `de5a4e600b57c8b9b0bd477695648c82351bc6dd`
- models `ten_tusscher_model_2006_{epi,endo,m}.cellml`

The seven main states are `V, Ki, Nai, Cai, Ca_ss, Ca_SR, R_prime`; the twelve
gates are `Xr1, Xr2, Xs, m, h, j, d, f, f2, fCass, s, r`. Main states use
Forward Euler and gates use Rush--Larsen, with all rates evaluated from the old
state. The CSV stores the seven states followed by the twelve gates.

All profiles use their published CellML initial conditions, `dt=0.005 ms`,
600 ms duration, `Istim=-52 pA/pF` for `10 <= t < 11 ms`, `Ksac=0`, and no
pre-pacing. EPI uses `G_to=0.294` and `G_Ks=0.392`; ENDO uses `G_to=0.073`,
`G_Ks=0.392`, and ENDO-specific `s` kinetics; M uses `G_to=0.294`,
`G_Ks=0.098`, and EPI-type `s` kinetics.

Generate any phenotype with:

```bash
python3 generate_ttp.py --profile epi \
  --output ionic_ttp_epi_trajectory.csv
```

Valid profiles are `epi`, `endo`, and `m`. Checkpoint `N` is the state after
exactly `N` completed updates.
