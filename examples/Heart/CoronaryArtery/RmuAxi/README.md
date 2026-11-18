# Truncation benchmark of reduced outlet closures

A reduced outlet resistance replaces a part of a vascular bed that is not
resolved. This benchmark resolves that part as well, so each closure can be
measured against a reference solution. The fluid is Carreau–Yasuda blood. The
flow is pulsatile and pressure-driven, and it reverses over part of the cycle.

## Geometry

The domain is a stepped axisymmetric pipe: coaxial tubes of radius $R_k$ and
length $L_k$, joined by abrupt contractions. The meridian half-plane is used,
with $x$ along the axis and $y = r \ge 0$.

A symmetric tree is a series of generations, each made of $N_k$ identical
tubes in parallel. The WRMS closure of a generation depends only on its wall
stress, and identical parallel tubes share that stress. A series of tubes is
therefore the axisymmetric surrogate of a symmetric tree.

## Runs

The **full** run resolves every tube. The outlet pressure is zero, and the
inlet traction is

$$
p_{\mathrm{in}}(t) = \bar P\,[1 + A\sin\omega t].
$$

Each **truncated** run stops at $x_T$ and closes the outlet with the uniform
traction $R(Q)\,Q$. The closures differ only in $R(Q)$:

| closure | $R(Q)$ |
|---|---|
| `tree` | $\sum_k 8\,\mu_{\mathrm{ap}}(4\lvert Q\rvert/\pi R_k^3)\,L_k/(\pi R_k^4)$ over the removed tubes |
| `single` | one compartment of radius $R_s$; its length is set so that it equals `tree` at $Q_0$ (the reduction of the paper) |
| `rest` | constant, `tree` at $Q_0$, multiplied by `-rmu_constant_scale` |
| `plateau` | constant, $\sum_k 8\,\mu_\infty L_k/(\pi R_k^4)$ |

## Resting operating point

The resting operating point is fixed before the solve, as in the paper:

- The narrowest tube carries the mean apparent shear rate `-rmu_gamma_ref`, which gives $Q_0$.
- $\bar P$ is the WRMS budget of the whole series at $Q_0$.

Changing `-rmu_scale` (the inlet radius) leaves the shear rates unchanged and
moves the Womersley number. This tests the quasi-steady assumption of the
closure.

## Numerical method

- **Discretisation:** Taylor–Hood P2/P1 elements and BDF1 time stepping.
- **Picard iteration:** each step iterates on the advecting velocity, the viscosity and the outlet resistance until convergence, so no coefficient is lagged.
- **Outlet traction:** the uniform traction is preconditioned by a local implicit term, $\beta R A\,(u\cdot n)(v\cdot n)$ minus its value at the previous iterate. The two cancel at convergence.

## Running it

```sh
cmake --build build --target RmuAxi -j
NP=4 examples/Heart/CoronaryArtery/RmuAxi/run_rmu_axi.sh
SCALES="1e-4 1e-3 3e-3" GAMMAS="100 400" ./run_rmu_axi.sh   # Womersley and shear-rate sweeps
```

`compare_rmu_axi.py` writes `summary.csv` and one figure of $Q_{\mathrm{in}}(t)$
per case. The table reports the errors of each closure against the full run
in the last cycle:

- the signed transport error and the transport error based on $\lvert Q\rvert$;
- the forward and reversed contributions to the signed error;
- the waveform RMS and maximum errors;
- the dissipation in the common region $x \le x_T$;
- the mean $\lvert\tau_w\rvert$ on the parent and stub walls.

`make_rmu_axi_mesh.py` writes both meshes. Its labels are 1 inlet, 2 outlet,
3 wall, 4 axis, 5 parent wall and 6 stub wall; the stub wall runs from the
first step to $x_T$.
