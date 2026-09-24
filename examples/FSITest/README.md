# FSITest — manufactured-solution check of the Robin–Robin FSI scheme

2D verification of Algorithm 1 of Burman, Durst, Fernández, Guzmán & Ruz,
*Robin–Robin loose coupling for incompressible fluid–structure interaction:
non-linear setting and nearly-optimal error analysis*, J. Sci. Comput. 104
(2025), with the correction iterations of its Appendix A. The meshes are
built inside the code, so no mesh files are needed.

| block    | discretization |
|----------|----------------|
| solid    | compressible neo-Hookean `P = mu (F - F^-T) + lambda ln J F^-T`, P1, mid-point scheme (16), stress at d^{n-1/2} (no Gonzalez correction), SNES |
| geometry | Ωⁿ = harmonic P1 lift of dⁿ⁻¹, wⁿ its time derivative |
| fluid    | Navier–Stokes (17), P1/P1, ALE, form (13) with Temam and interface convective terms, GCL term s_Ω (15), Carreau–Yasuda viscosity at the lagged velocity; projected VMS convection, orthogonal grad-div subscale, PSPG (or pure Brezzi–Pitkäranta) — same integrators and same rheology as `CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq` |
| coupling | Robin terms on the **reference** interface with one P1 trace mass matrix; λ = σ_f n_f updated by (18); `K > 1` = Appendix A corrections; α = γ √(ρ_s E) |

Geometry: fluid (0,1)², solid (0,1)×(1,1+H), interface y = 1.
The manufactured fields and every source/mismatch term come from
`generate_mms.py` (sympy), which writes `ManufacturedSolution.h`.

## Run

```sh
cmake --build build -j --target FSITest
mkdir -p run_fsitest && cd run_fsitest
../build/examples/FSITest/FSITest                  # 4 levels, N = 8..64, dt = 0.04..0.005
python3 ../examples/FSITest/plot_convergence.py    # -> fsitest_convergence.png
```

Options: `-fsitest_levels`, `-fsitest_n0`, `-fsitest_dt0`, `-fsitest_T`,
`-fsitest_gamma`, `-fsitest_coupling_iterations K`,
`-fsitest_lambda update|residual|pointwise`,
`-fsitest_cy_mu0 -fsitest_cy_lambda -fsitest_cy_n -fsitest_cy_a -fsitest_cy_reg`
(Carreau–Yasuda; μ_inf is `-fsitest_mu_f`, and `cy_mu0 = mu_f` is Newtonian),
`-fsitest_pspg_residual 0` (pure Brezzi–Pitkäranta), `-fsitest_vms_scale`,
`-fsitest_graddiv_scale`, `-fsitest_pgp_scale`, `-fsitest_xdmf 1`
(ParaView output in `results_fsitest/`), and the material parameters
`-fsitest_rho_f -fsitest_mu_f -fsitest_rho_s -fsitest_lambda_s -fsitest_mu_s`.

`-fsitest_lambda` (all three give the same fluid stress in the limit; the
first one is what the solvers use):

- `update` (default): eq. (18), λⁿ = λⁿ⁻¹ + α(ḋ^{n−½} − uⁿ). Lives in the
  trace space, costs nothing, evaluates no gradient, and reproduces the fluid
  residual exactly — which is what the energy estimate of Theorem 1 uses;
- `residual`: eq. (19), λⁿ from the fluid variational residual tested with the
  interface basis functions. Prints `max|lambda(19) - lambda(18)|` ≈ 1e-15:
  the two are the same object, as the paper shows;
- `pointwise`: strong trace σ_h n of the discrete stress. Consistent but not
  conservative (its nodal values do not reproduce the fluid residual), and it
  depends on evaluating ∇u_h at the nodes.

## Reference results (T = 0.4, Δt ∝ h)

| N | Δt | ‖u−u_h‖ | ord | ‖p−p_h‖ | ord | ‖d−d_h‖ | ord |
|---|---|---|---|---|---|---|---|
| 8  | 0.04  | 4.10e-3 | –    | 3.34e-2 | –    | 1.04e-3 | –    |
| 16 | 0.02  | 1.77e-3 | 1.21 | 2.34e-2 | 0.51 | 4.25e-4 | 1.30 |
| 32 | 0.01  | 8.24e-4 | 1.10 | 1.34e-2 | 0.81 | 1.85e-4 | 1.20 |
| 64 | 0.005 | 4.02e-4 | 1.04 | 6.90e-3 | 0.96 | 8.68e-5 | 1.09 |

### Non-Newtonian (Carreau–Yasuda)

μ(γ̇) = μ_∞ + (μ_0 − μ_∞)(1 + (λγ̇)^a)^((n−1)/a), γ̇ = √(γ_reg² + 2D:D) — the
same law, the same regularization and the same lagging as the Heart solvers.
The manufactured source terms and the interface data are generated from the
full non-linear stress, so nothing about the law is hidden from the test.

```sh
../build/examples/FSITest/FSITest -fsitest_cy_mu0 1.0 -fsitest_cy_lambda 50
```

| N | Δt | ‖u−u_h‖ | ord | ‖p−p_h‖ | ord | ‖d−d_h‖ | ord |
|---|---|---|---|---|---|---|---|
| 8  | 0.04  | 2.96e-3 | –    | 3.92e-2 | –    | 1.04e-3 | –    |
| 16 | 0.02  | 1.23e-3 | 1.26 | 2.59e-2 | 0.60 | 4.25e-4 | 1.29 |
| 32 | 0.01  | 5.70e-4 | 1.11 | 1.46e-2 | 0.83 | 1.87e-4 | 1.19 |
| 64 | 0.005 | 3.02e-4 | 0.92 | 7.42e-3 | 0.97 | 8.77e-5 | 1.09 |

With those parameters the viscosity the fluid sees spans μ ∈ [0.26, 0.92]
(the run prints the range per level), i.e. a factor of 3.5 across the domain.

First order in time with Δt ∝ h, as expected. `-fsitest_lambda pointwise`
converges at the same rate with a slightly different constant (3.63e-4 at
N = 64); before the H1 gradient fix in `src/Rodin/Variational/H1` it stalled
at ~3.6e-3, because the nodal projection of ∇u_h returned zero at the
collapsed vertex of the Dubiner map.
