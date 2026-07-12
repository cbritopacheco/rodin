---
name: run-adaptation-examples
description: Run the WNGIR / level-set adaptation and shape-optimization examples with known-good parameters. Use when reproducing, debugging, or benchmarking the moving-interface pipelines.
---

# Run the adaptation / shape-optimization examples

Examples build under `build/examples/<Theme>/<Name>`. They dump
`*.h5`/`*.xdmf`/`*.log` into the CWD — **always run from a scratch
directory**, never the repo root, and never commit these artifacts.

## CLI convention

Flags are `--name=value`. A bare `n=12` is **silently ignored** (the run
falls back to defaults and is just mysteriously slow or wrong). Check flag
spelling first when a parameter seems to have no effect.

## Known-good 2D WNGIR cantilever baseline (user-confirmed)

```sh
cd "$(mktemp -d)" && /path/to/build/examples/ShapeOptimization/LevelSetWNGIRCantilever2D \
  --ell=0.5 --iters=200 --n=25 \
  --dt-max=0.02 --reinit=none --redistance-every=20 --redistance-mode=fmm \
  --wngir-gamma-m=0.25 --wngir-gamma-h=0.25 --wngir-ell=1.5 \
  --wngir-rms-tol=1e-4 --wngir-sup-tol=5e-4 --wngir-steps=120
```

Expected: converges to a proper cantilever (material must NOT flood the
domain). Non-negotiables: objective line search stays ON (never pass
`--objective-linesearch=0`); `--ell=0.5`. The `--wngir-ell` bulk length
scale is regime-split: ~1.5 for shape optimization, ~8 for
reconstruction/3D quality.

## 3D / PETSc-backed runs

- `LevelSetWNGIRCantilever3D`: `--n=12` is the fast repro size (`--n=6` too
  coarse — interface empty; `--n=20` is slow). Run PETSc-backed 3D paths
  with `OMP_NUM_THREADS=1` (OpenMP hazards).
- `-wngir_pc_type lu` only affects the PETSc (3D) path; in 2D the step solve
  is native Eigen CG and the flag is inert.

## WNGIR reconstruction / sweep / advection

`build/examples/Geometry/LevelSetWNGIR{Reconstruction,Reconstruction3D,Sweep,Advection}`.
Advection MUST run with redistancing (`--phi-redistance=fmm-moved`) — the
carried level set drifts catastrophically without it.

## Judging runs

Judge by the per-stage metrics the examples print (fit rms/max, quality
minima, inverted count) — not by eyeballing output. The benchmark suite
(`RodinBenchmarks`; Connectivity/MeshIO/P1/Poisson/UniformGrid cases) is
the oracle for core-path performance claims.
