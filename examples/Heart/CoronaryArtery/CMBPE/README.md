# CMBPE campaign: reduced closures for generalised-Newtonian outlets

These scripts reproduce the numerical study of the revised manuscript on
reduced outlet closures, $\mathcal{R}_\mu$. It runs on one node with a shared
core budget (`CORES`, 60 by default).

| step | script | what it runs |
|---|---|---|
| 0 | `00_build.sh` | configures `$BUILD` if it has no `CMakeCache.txt`, then builds `RmuAxi` and `CoronaryArtery` |
| 1 | `01_validation.sh` | axisymmetric truncation benchmark (`RmuAxi`), one core per run |
| 2 | `02_coronary.sh` | coupled LV-0D / coronary-3D runs, `NP` cores per run |

`queue.py` schedules the jobs of each step within `CORES`. A job starts when
its cores are free and its dependencies have finished. Every job writes
`log.txt` in its directory and `.done` when it succeeds. Relaunching a step
skips the finished jobs.

```sh
cd examples/Heart/CoronaryArtery/CMBPE
BUILD=~/rodin/build ./00_build.sh
CORES=60 OUT_ROOT=/scratch/cmbpe ./01_validation.sh
CORES=60 NP=10 OUT_ROOT=/scratch/cmbpe ./02_coronary.sh
```

## Validation (`01_validation.sh`)

The geometry is a stepped pipe with tubes `1:10, 0.6:30, 0.35:30`, in units of
the inlet radius. The flow is pressure-driven, pulsatile and reversing
($A = 1.2$). The full solution is the reference.

Each of the following cases is solved on the truncated geometry with the
closures `tree`, `single`, `rest` and `plateau`:

- three bloods;
- three inlet radii (`SCALES`, $10^{-4}$, $10^{-3}$ and $3\cdot10^{-3}$ m), which keep the operating shear rates and raise the Womersley number;
- two reference shear rates (`GAMMAS`, 100 and 400 s$^{-1}$).

Two more runs check the convergence of the reference: the time step halved
(`conv_dt`) and the mesh refined (`conv_h`).

There are 100 two-dimensional runs. The outputs are:

- one `summary.md` per root, with the errors against the full solution;
- the $Q_{\mathrm{in}}(t)$ figures;
- `reference_convergence.md`.

## Coronary (`02_coronary.sh`)

The bed resistances are calibrated with the healthy blood at its resting
operating point (`-coronary_calibration reference`). The drainage pressure is
$P_{RA} = 700$ Pa and $\alpha = 0.75$, as in the manuscript.

| group | runs |
|---|---|
| `pra700/<blood>` | `rmu`, `rest` (constant, frozen at rest), `plateau` ($\mu_\infty$) |
| `pra2400/hyperviscous` | the same three closures at the elevated drainage pressure |
| `<group>/rest_matched` | `rest` scaled by secant iteration until its mean $\lvert Q\rvert$ equals that of `rmu` |
| `tree/<blood>` | `rmu`, `rest` with a three-generation arteriolar limb (`TREE`, illustrative) |
| `mesh/{coarse,medium,fine}` | healthy `rmu`, `rest`; the coarse and fine meshes read the medium calibration |
| `dt/dt{1e-3,5e-4}` | hyperviscous at 2400 Pa, `rmu` and `rest` |

`compare_coronary_closures.py` writes `summary.md`. For every group it reports:

- the cycle means of $Q$ and $\lvert Q\rvert$;
- the compressive and relaxed transport;
- the flow-weighted viscosity factor;
- the dissipated energy, the mean $\lvert\tau_w\rvert$ and the driving pressure drop;
- the differences of each closure against `rmu`, split into compressive and relaxed parts.

### Cost

A run of 3 cycles at $\Delta t = 10^{-3}$ s has 2550 steps. Read the time per
step in the first logs (`grep Step */log.txt`) before scaling up.

With `NP=10`, six runs proceed at once. The fine mesh uses `NP_FINE=20`. Set
`DTS="5e-4 2.5e-4"` for a third time-step level, which costs four times the
base run.

## Meshes

`meshes/` holds the three coronary meshes as plain files, not in Git LFS:

| mesh | tetrahedra | outlet area (mm²) | origin |
|---|---|---|---|
| `coronary_medium.mesh` | 83,854 | 12.421 | the medium mesh of the manuscript (`coronaria3d.mesh`, branch `CoronaryPoroelastic`) |
| `coronary_coarse.mesh` | 24,182 | 11.341 | `mmg3d -hmax 1.0 -hausd 0.3` from the medium mesh |
| `coronary_fine.mesh` | 208,807 | 12.399 | `mmg3d -hmax 0.45 -hausd 0.1` from the medium mesh |

The mesh coordinates are in mm. The labels are: wall 2, inlet 4, outlets 7,
8, 9, 10, 14 and 15. `make_meshes.sh` regenerates the coarse and fine meshes.

## Companion benchmarks and analysis

- [Y-bifurcation truncation](../TruncationY/README.md) compares the split
  between two daughter branches against a fully resolved reference in 2D
  and 3D. The archived meshes, calibration and waveform data accompany the
  solver, runner and figure scripts.
- [Coronary outlet comparisons](analysis/README.md) contains a prescribed-pressure
  0D replay and tools for comparing completed coronary runs on a common
  cycle grid and wall mesh.
