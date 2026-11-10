# Manual validation scripts

## Expanded-prior warm starts

`expanded_prior_real_data_trial.py` compares an initial narrow-prior fit, an
expanded-prior warm fit, and an independent expanded-prior cold fit on saved
light curves. It retains the observation-specific priors, seeds, and comparison
criteria from the original trial. The light curves are local inputs and are
not included in this repository.

Run it from the repository root with EXOTIC installed in your environment:

```bash
python scripts/expanded_prior_real_data_trial.py \
  --dataset-path 0 "/data/NormalizedFlux_WASP-43 b_2026-03-20.txt" \
  --output /data/expanded-prior-trial.json
```

Use the index corresponding to the original observation:

| Index | Observation preset | Input format |
| --- | --- | --- |
| 0 | WASP-43b full transit, 2026-03-20 | Normalized flux |
| 1 | TrES-5b full transit, 2026-07-16 | Final light curve |
| 2 | KELT-20b partial transit, 2026-07-15 | Normalized flux |

Normalized-flux files are CSV tables with `BJD`, `Norm Flux`, `Norm Err`, and
`AM` headers. Final-light-curve files are CSV tables with time in column 1,
flux in column 3, flux error in column 4, and airmass in column 6; comment lines
start with `#`.

Repeat `--dataset-path INDEX PATH` to supply multiple datasets. By default,
the script runs the supplied indices in argument order. The existing
`--dataset-index INDEX` option can select a subset of those inputs and can
also be repeated. A path is required for every selected dataset. Relative
paths are relative to the directory from which you run the command.

`--live-points` and `--seed` retain their existing defaults of 100 and 4729.
The JSON output retains the fit summaries, comparison diagnostics, and
pass/fail result. Loading or fitting failures are reported directly.
