# Verification against exact pipe-flow solutions

There are two tests:

1. **Steady sPTT pipe flow** (`run_sptt_pipe.sh`). This checks space only, since the solution is steady.
2. **Pulsatile Oldroyd-B pipe flow** (`run_womersley.sh`). This one is time-dependent and exact for any amplitude. It runs a space study with BDF2, mesh by mesh, and a time study comparing BDF1 and BDF2. It works for the 3D driver and for the axisymmetric one.

The derivations are in `sptt_pipe_solution.tex` and `womersley_ob_solution.tex`,
both ready to `\input` into the article; the references are in
`sptt_pipe_solution.bib`.

| file | role |
|---|---|
| `sptt_pipe.py` | Reference solution `SPTTPipe(R, eta_s, eta_p, lam, eps)`: `set_flow_rate(Q)`, `u(r)`, `sigma_xr(r)`, `sigma_xx(r)`, `G`, `wall_shear_stress`. Run `python3 sptt_pipe.py` to self-check it against the Oldroyd-B and Oliveira–Pinho closed forms. |
| `compare_sptt_pipe.py` | Reads the last step of a run (all MPI ranks) and prints the relative L2 errors of u and σ, the error of the pressure gradient, and the wall shear stress error over 5 ≤ x/D ≤ 10. |
| `run_sptt_pipe.sh` | Runs the mesh sequence for Oldroyd-B (ε = 0) and sPTT (ε = 0.25), then prints the table with observed orders. |
| `sptt_pipe_solution.tex/.bib` | The analytical solution and the test set-up for the paper. |
| `womersley_ob.py` | Exact pulsatile Oldroyd-B pipe flow `WomersleyOB(R, rho, eta_s, eta_p, lam, omega, Ubar, A)`: `u(r,t)`, `sigma_xr`, `sigma_xx`. Run `python3 womersley_ob.py` to check the momentum and constitutive residuals. |
| `compare_womersley.py` | Errors of one run (3D or axisymmetric) against the exact solution over the last cycle. With `--ref`, the difference with another run on the same mesh, for self-convergence in time. |
| `run_womersley.sh` | Space study (BDF2, mesh by mesh) and time study (BDF1 and BDF2 at 16 to 128 steps per cycle), with tables, LaTeX tables and a figure. |
| `make_pipe_axi_mesh.py` | Straight axisymmetric pipe meshes for the axisymmetric driver; needs no gmsh. |
| `plot_sptt_pipe.py`, `plot_womersley.py` | Figures. |
| `womersley_ob_solution.tex` | The pulsatile solution, the BDF2 scheme and the test set-up for the paper. |

## Running on pc9

```
cd ~/rodin
bash examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit/verification/run_sptt_pipe.sh
```

Options are set through environment variables:

- `NP` (default 24): processes for the finest mesh. With `PARALLEL=1`, the default, all cases run at once and the coarser meshes get NP(rf/rf_max)³ processes, at least 2, which comes to 76 cores in total. `PARALLEL=0` runs the cases one after another with NP each.
- `RF_LIST` (default `0.7 1.0 1.4 2.0`). The rf 0.5 mesh is also in the repository, but it has only two elements per radius and is pre-asymptotic.
- `CASES` (default `OB sPTT`)
- `OUT` (default `~/Results/axi_jfm/verification`)
- `EXE3D`: path to the driver, if the build directory is not `build/`. On pc9, use `EXE3D=$PWD/build-gcc12/examples/...`.

A case that has finished leaves a `DONE` file in its directory, so a second call
only recomputes the table. The meshes are read from
`resources/examples/viscoelastic_fluids/pipe_verification_rf*.mesh`. If one is
missing, it is generated with gmsh, which needs the gmsh Python module.

## Set-up

- **Geometry.** Straight pipe of length 14D: the inlet is at x = −2D and the outlet at x = 12D. The mesh is the healthy (`H`) case of `make_lesion_mesh_3d.py`, built with `Lu = 2`, `Ld = 12` and every size parameter set to 0.1, scaled by `rf`.
- **Fluid.** η_s = η_p = 1.795 mPa s, so β = 0.5 and ρ = 1060. Wi = λŪ/D = 0.5, which gives λ = 0.236 s. ε is 0 for Oldroyd-B and 0.25 for sPTT; with ε = 0.25 the PTT factor reaches f = 2.4 at the wall.
- **Flow.** Re = 10, so Ū = 8.5 mm/s. The inflow is steady (A = 0), the Wo = 1 period is only used as a clock, and the run lasts 2 cycles of 20 steps. Newton converges to |dψ| ≈ 1e-14 at the last step.
- **Boundary conditions.** The inlet imposes the Newtonian parabola and ψ = 0, so the viscoelastic profile has to develop inside the pipe. The outlet carries zero traction.
- **Mesh sequence.**

| rf | h/D | vertices | tetrahedra |
|---|---|---|---|
| 0.5 | 0.2 | 2.1k | 7.6k |
| 0.7 | 0.143 | 5.0k | 20k |
| 1.0 | 0.1 | 12k | 53k |
| 1.4 | 0.071 | 30k | 141k |
| 2.0 | 0.05 | 79k | 399k |

  `check_3d_mesh.py` reports a volume deficit of the polygonal pipe of −2.0 %, −1.0 %, −0.5 %, −0.25 % and −0.12 % for these meshes. That deficit is O(h²) and is part of the discretisation error being measured.

## How the errors are measured

The reference solution is evaluated at the flow rate the discrete solution
actually carries through the region, Q_h = (1/L)∫u_x dV. This separates the
error of the developed flow from the inflow error, which comes from the
polygonal inlet and the developing region. The difference
`err_Q_region_vs_nominal` is printed for reference.

## Known issue: TAWSS = 0 at partition interfaces

In the 3D output, `TAWSS`, `OSI` and probably `shearStress` are zero at some
wall vertices that lie on MPI partition interfaces. On 2 ranks this affected 27
of 484 wall vertices in the region. The likely cause is that
`GridFunction::project(Region::Boundary, ..., onWall)` writes only the facets
local to each rank, and no ghost-to-owner update follows. The maxima written
to the CSV and the log (`maxTAWSS`, `maxTAWSSLesion`) are not affected. The
fields shown in ParaView are affected, and so is anything averaged over wall
nodes. `compare_sptt_pipe.py` excludes those vertices and counts them in
`n_wall_zero`. Until the driver is fixed, do not average TAWSS or OSI from the
3D XDMF without excluding exact zeros, and do not show the TAWSS field from
runs with many ranks without checking for this.

## Expected behaviour of the wall shear stress

The driver computes the WSS as the tangential part of (2η_s ε(u_rec) + σ)·n
and projects it onto the wall. Here ε(u_rec) is the L2 (nodal) recovery of the
P1 velocity gradient. At wall nodes this recovery averages the gradient of
elements whose centroids lie about h/4 inside the wall, so the solvent part of
the WSS converges at O(h) at best, while u and σ converge at O(h²) in L2. On
the two coarsest meshes, for sPTT, the solvent part is low by 6.1 % at
h/D = 0.2 and by 5.2 % at h/D = 0.143. The polymer part is within 2 %. Use the
finer levels to quote a WSS accuracy. If the WSS error is the limiting one, a
consistent boundary flux would converge at the rate of the solution: compute
the traction from the residual of the momentum equation, the reaction-force
method.

## Time integrator: BDF2 (`-al_time_order 2`)

The 3D and axisymmetric drivers take `-al_time_order 1|2`. The default, 1, is
the BDF1 scheme the production runs used, and it is unchanged bit for bit. With
BDF2:

- **Time derivative.** It is (3x^{n+1} − 4x^n + x^{n−1})/(2Δt) in the momentum equation and in the ψ (and ψ_θθ) equation.
- **Convection and viscosities.** The Galerkin convection advects with u* = 2u^n − u^{n−1}. The Carreau–Yasuda and gsPTT viscosities (axisymmetric driver) are also evaluated at u*.
- **Newton.** It starts from (u*, p*, ψ*). It converges in 2 to 4 iterations, like BDF1.
- **First step.** It takes x^{−1} = x^0. This is exact here because the runs start from rest with the (1 − cos) ramp.

**Stabilisation projections** (`-al_bdf2_extrapolate <mask>`). The split-OSS projections are explicit. With BDF1 they are lagged at step n. Each lagged projection leaves an O(Δt) error times its weight. The mask picks which ones are taken at the extrapolation instead:

| projection | bit | weight | at the extrapolation |
|---|---|---|---|
| Π[div u] | 1 | – | unstable |
| Π[∇p] | 2 | α1 = O(h²) | unstable |
| Π[div σ] | 4 | α1 = O(h²) | unstable |
| Π[ε(u)] | 8 | α3 ≈ η_p/2 | stable |
| Π[(u^n·∇)ψ] | 16 | ≈ λ/4 | stable for Oldroyd-B, unstable for the gsPTT blood model R |
| Π[u^n·∇ψ_θθ] | 32 | (axisymmetric) | – |

- **Default (8).** It extrapolates only Π[ε(u)], whose weight does not vanish with h. Mask 24, which adds Π[(u^n·∇)ψ], was stable for Oldroyd-B but diverged within the first period for R (axisymmetric pipe, Re 100, both at 40 and at 160 steps per cycle). Mask 8 and mask 0 are stable there.
- **Unstable masks.** In the axisymmetric pulsatile test, masks 2, 3, 6 and 7 diverged after 2 to 3 periods. A spurious radial velocity or hoop stress grew near the axis. Recomputing the projections at every Newton iterate instead is consistent but converges linearly with a factor of about 0.7 (12 iterations were not enough), so it was dropped.
- **What stays lagged.** Π[∇p], Π[div σ] and Π[div u] stay at step n, with an O(h²Δt) error. On coarse meshes this error limits the observed order (see the table below). It decreases as h², so the order tends to 2 under refinement.

Axisymmetric time study (pulsatile Oldroyd-B, h/D = 0.1, self-convergence d = ‖x_Δt − x_Δt/2‖ at t = 4T):

| T/Δt | BDF1 d_u | BDF1 d_σ | BDF2 d_u | BDF2 d_σ |
|---|---|---|---|---|
| 16 | 9.4e-3 | 2.3e-2 | 9.7e-3 | 1.1e-2 |
| 32 | 3.5e-3 (p 1.4) | 1.3e-2 (p 0.85) | 3.8e-3 (p 1.3) | 4.7e-3 (p 1.2) |
| 64 | 1.5e-3 (p 1.2) | 6.9e-3 (p 0.9) | | |

On this coarse mesh (h/D = 0.1), BDF2 at a given Δt has a time error in σ two
to three times smaller than BDF1. In u the two schemes are alike, because the
lagged Π[∇p] dominates. The error against the exact solution at T/Δt = 64 is
1.5 % in u and 5 % in σ, and it is spatial for both schemes. The finer meshes,
on pc9, show how the order recovers.

## Formulation (`-al_formulation psi_rb|exp|psi`)

In all three formulations, momentum and continuity carry the split OSS with dynamic subscales, as in Moreno, Codina and Baiges (CMAME 367, 2020):

- **Subscales.** The convective, pressure-gradient and div σ subscales carry the memory ρ/Δt. The pressure subscale is algebraic (grad-div on the orthogonal part of div u).
- **α1.** It uses η0 = η_s + η_p/f, where h is the cell length along u. The constants are c1..c4 = 4, 2, 4, 0.25 (`-al_c1..4`).
- **Storage.** The subscales are kept at the quadrature points and updated when a step is accepted. `-al_dynamic_subscales 0` gives the quasi-static variant.

The formulations differ in the constitutive equation:

- **`psi_rb` (default).** The ψ equation (Fattal and Kupferman), with BDF on ψ:
  R_ψ = ∂_t ψ + (u·∇)ψ + N(ψ, ∇u), N = −Dexp⁻¹[∇u E + E ∇uᵀ − c ε(u)] + (f/λ)(I − e^{−ψ}).
  - **Stabilisation.** It is the OSS of the full residual. A log-conformation subscale ψ̃ = τ′(ψ̃ⁿ/Δt − (R_ψ − Π[R_ψⁿ])), with τ = (λ/2η_p)α3, enters in two places:
    - the ψ equation, as ⟨ψ̃, −(uⁿ·∇)χ⟩;
    - the momentum, as ⟨σ̃, ε(v)⟩ with σ̃ = (η_p/λ0) Dexp[ψ̃]. This term gives the u–ψ compatibility.
  - **Approximation.** The reaction part of the adjoint is left out.
  - **Newton.** It is exact and quadratic.
- **`exp`.** Moreno, Codina, Baiges and Castillo (CMAME 354, 2019). The σ equation is written for E = exp(ψ) and tested with χ:
  R = f(E − I)/(2λ0) − ε(u) + (λ/2λ0)[∂_t E + (u·∇)E − (∇u E + E ∇uᵀ) + 2ε(u)].
  Its stabilisation is ⟨σ̃, A*(χ) + ε(v)⟩, with σ̃ = α3′(λ/(2η_pΔt) σ̃ⁿ − (R − Π[Rⁿ])).
- **`psi`.** The former scheme: the ψ equation with split OSS of (u·∇)ψ and ε(u). This adds α3⟨ε(u) − Π[ε(u)], ε(v)⟩ with α3 ≈ η_p/2, which for blood is several times η_s. It acts whether or not the equation is resolved.

Why `psi_rb`:

- **Consistent and low-dissipation.** Its stabilisation acts only on the part of R that the mesh does not resolve. Start-up in the S75 pipe (Re 600, fluid R) shows the difference. At t = 1 ms the diffusion length √(νt) is 0.1 R, so the core should be a plug flow:
  - `exp` and `psi_rb`: Δp = 2.9 Pa, max|u|/Ū = 4.7. These are the inertial values: the plug through the 75 % stenosis, with ρ L dŪ/dt of a few Pa.
  - `psi`: Δp = 0.29 Pa, max|u|/Ū = 2.2. This is far from that. Its stabilisation differs in two ways that both act on this transient: the quasi-static α1 with lagged Π[∇p], and S3.
- **Well scaled at high Wi.**
  - In the E form, the equation for an eigenvalue l of ψ is scaled by e^l. At blood Weissenberg numbers (λγ̇ ~ 10³), the small eigenvalues are therefore not determined.
  - In the S75 xcoarse test, `exp` had ψ_min = −7.5 at t = 0.042 s and −20 at t = 0.052 s. Newton then stalled at every Δt.
  - With the same data, `psi_rb` ran to t = 0.115 s with no rejected step and 2 Newton iterations per step, which is past the t = 0.108 s where the earlier 3D run had failed. Its ψ_min stayed between −4.5 and −6.6. At that point ψ_max was 9.9 and max|σ| was 205 Pa, the latter levelling off. The test was stopped there.

Steady sPTT pipe, rf 0.5 (h/D = 0.2):

| formulation | e_u | e_σ | e_G | e_τw | Newton iterations / steps |
|---|---|---|---|---|---|
| psi_rb | 6.18e-2 | 1.40e-1 | 1.14e-1 | −4.4 % | 146 / 52 |
| exp | 6.16e-2 | 1.62e-1 | 1.14e-1 | −4.3 % | 101 / 40 |
| psi | 6.38e-2 | 1.39e-1 | 1.15e-1 | −4.4 % | 92 / 40 |

Assembly per Newton iteration costs about 1.5 to 1.7 times that of `psi`, for both `psi_rb` and `exp`.

## Pulsatile Oldroyd-B test

```
DRIVER=3d  bash run_womersley.sh        # pc9: NP=8 MAXJOBS=10 runs ten at a time
DRIVER=axi bash run_womersley.sh        # axisymmetric driver, meshes generated
```

- **Parameters.** Re = 10, Wo = 3, q = 1 + 0.5 sin ωt, β = 0.5, λω = 1. The flow is ramped over half a cycle, then 4 cycles are run, and the last one is compared at 8 phases over 5 ≤ x/D ≤ 10. The periodicity residual is about 1e-8.
- **Space study.** BDF2 with 128 steps per cycle on each mesh. The 3D meshes are rf 0.7, 1.0 and 1.4; the axisymmetric ones are h/D = 0.1, 0.05 and 0.025.
- **Other masks.** `BDF2_EXTRAP=<mask>` passes `-al_bdf2_extrapolate`. Use a different `OUT` for each mask.
- **Time study.** On one mesh, BDF1 and BDF2 with 16, 32, 64 and 128 steps per cycle. The self-convergence differences d = ‖x_Δt − x_Δt/2‖ at t = 4T exclude the spatial error.
- **Results.** `tables.txt`, `womersley_table.tex` (two complete LaTeX tables), `womersley.pdf` and `summary.csv`.
