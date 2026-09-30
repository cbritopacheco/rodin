# Iteration-Resolved Calibration

The C++ WNGIR defaults are 30 outer iterations and 15 inner Newton corrections.
The calibrated coefficient triple and convergence tolerances are unchanged.
Both campaign runners use these caps unless explicitly overridden.

For the next campaign, use a fresh output directory and pass
`--steps=30 --barrier-max-iters=15 --log-iterations` to either runner.
Logging is opt-in so an existing campaign with explicit caps can continue with
its original protocol and binary. Do not rebuild its executable while it runs.
The older broad launcher explicitly requests 20/10 and is not a 30/15 launcher.

Each case saves a lossless trace in `<out-dir>/iteration_logs/`, with its case
identity and command in the first two lines. The terminal CSV schema is unchanged.
Inner rows contain Newton correction and iterate norms, relative correction,
step factor, linear iteration count/error, and convergence status; failed linear
solves and infeasible corrections also produce rows. Geometry rows contain
complete-interface RMS/maximum distance, normal error, Jacobian and distortion,
cumulative Newton corrections, and elapsed solver wall time including setup.
They cover the initial state, each accepted outer state, and the terminal state.
The `outer` index on geometry rows counts accepted outer steps; inner rows use
zero-based outer attempt indices. Inner work does not itself move the accepted mesh.

Tracing evaluates geometry additionally and therefore incurs measurement cost.
It does not change the fitting energy, damping, or stopping conditions. Timings
include this logging/validation cost and must be compared with like-for-like runs.

To compute first-hit counts after selecting independent per-lobe reference
coefficients, supply a JSON map from lobe index to `C_infty`:

```sh
python3 experiments/wngir_calibration/summarize_iteration_logs.py \
  /path/to/campaign/iteration_logs --coefficients /path/to/reference_coefficients.json \
  --order=2 --output=/path/to/campaign/first_hits.csv
```

The P1 target is `D_infty <= C_infty(lobes) * h^2`, with sampled Jacobian greater
than 0.01 and distortion strictly below 10, matching the runners' admissibility
profile. Coefficients are required, not inferred from candidate terminal errors.
The summary reports inner attempts (including failed linear solves), their mean,
median, and maximum per attempted outer solve, damping/failure counts, and first-hit outer
count, cumulative completed Newton corrections, and time. Missing hits remain
missing rather than being replaced by the iteration cap. A sampled target hit at
one resolution does not establish an error order or certify a Hausdorff bound.
