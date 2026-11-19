# Atrium: viscoelastic and thixotropic blood on the fixed fluid domain

`AtriumViscoelastic` solves the flow of blood in the fluid region of the Atrium
mesh for four constitutive descriptions that differ by one ingredient each. It
is the log-conformation solver of
`examples/viscoelastic_fluids/ArterialLesion3D_viscous_logarithmic_implicit`
(branch `model/viscoelasticFluid`) with the domain and the openings of
[Atrium/Fluid](../Fluid/README.md). The wall is fixed: the solid and the
fluid–structure coupling of the Atrium example are not part of it.

## Fluids

| Tag | `-av_fluid` | Unknowns | Steady shear viscosity | Memory |
|---|---|---|---|---|
| N | `newtonian` | $u,p$ | $\eta_\infty$ | none |
| CY | `cy` | $u,p$ | $\eta_{CY}(\dot\gamma)$ | none |
| VE | `owens -av_owens_equilibrium 1` | $u,p,\psi$ | $\eta_{CY}(\dot\gamma)$ | elastic |
| TVE | `owens` | $u,p,\psi,n$ | $\eta_{CY}(\dot\gamma)$ | elastic and structural |

TVE is the model of Owens (J. Non-Newtonian Fluid Mech. 140, 2006) as already
implemented in the arterial driver. The extra stress is that of Hookean
dumbbells of constant modulus $G=(\eta_\infty-\eta_{pl})/\lambda_H$ whose
relaxation time depends on the mean aggregate size $n$,

$$
\lambda(n,\dot\gamma)=\frac{n\lambda_H}{1+\left(\tfrac12 b\,n(n-1)+a\right)\lambda_H},
\qquad
\frac{Dn}{Dt}=-\tfrac12\,b(\dot\gamma)\,(n-n_{st})(n+n_{st}-1),
$$

with $n_{st}(\dot\gamma)=\frac{\eta_{CY}(\dot\gamma)-\eta_{pl}}{\eta_\infty-\eta_{pl}}\left(1+\tfrac32 a\lambda_H\right)$
and $b=a/(n_{st}(n_{st}-1))$, so that $\eta_{pl}+G\lambda(n_{st},\dot\gamma)=\eta_{CY}(\dot\gamma)$.

VE is the limit of fast kinetics of the same model: $n=n_{st}(\dot\gamma)$ at
every point and time, hence $\lambda=(\eta_{CY}(\dot\gamma)-\eta_{pl})/G$, a
White–Metzner fluid. It has the modulus and the steady shear viscosity of TVE
and no transported structure. The three pairs therefore isolate one effect
each: CY−N shear thinning, VE−CY elasticity, TVE−VE the memory of the
structure. `gsptt` (constant $\lambda$, steady viscosity $\eta_{CY}$ through
its solvent) and `sptt` remain available.

Both Owens fluids have zero second normal stress difference in viscometric
flow, as any Oldroyd-B type model; `postprocess_atrium.py` reports $N_2$ in a
local frame as a diagnostic of the three-dimensional stress, not as a material
function.

## Constants

| Constant | Value | Source and status |
|---|---|---|
| $\lambda_H$ | 5 ms | Moyers-Gonzalez, Owens & Fang, J. Fluid Mech. 617 (2008), Table 1 |
| $a(\dot\gamma)$ | 1 s⁻¹ at rest, cubic to 2 s⁻¹ at $\dot\gamma_c=5.78$ s⁻¹, cubic to 0 at 900 s⁻¹ | same paper, Appendix A; the coefficients $a_{1,k}$, $a_{2,k}$ expand to these cubics |
| $\eta_{pl}$ | 1.2 mPa s | plasma at 37 °C; the paper uses 1.0 mPa s |
| $\eta_{CY}$ | 72, 4.9 mPa s; 6.1529 s; 0.2081; 1.4173 | the law of Atrium/Fluid; Cho & Kensey (1991) by options |
| $n_{in}$ | 1.5 | assumption: blood ~1 s after the pulmonary capillaries, $n\approx1+at/2$ |

Two departures from the published model follow from the steady viscosity.
The paper fits a Cross law ($\eta_0+\eta_N=33.6$ mPa s, $\eta_\infty+\eta_N=3.34$
mPa s, $m=2.1238$, $\beta=0.7014$ s$^m$) together with $a(\dot\gamma)$. That
Cross law has a decreasing shear stress near $\dot\gamma\approx2$ s⁻¹
($d(\eta\dot\gamma)/d\dot\gamma=-1.2$ mPa s at $\beta\dot\gamma^m=(m+1)/(m-1)$),
which is ill-posed for the generalised Newtonian member of the hierarchy. A
monotone Carreau–Yasuda law is used for all three fluids instead, and the
aggregation rate is kept. The time scale of the structure, $\approx n_{st}/a$,
then changes with $n_{st}$: 19 s at rest here against 14 s with the Cross law.

With these constants ($G=0.74$ Pa):

| $\dot\gamma$ (s⁻¹) | $\eta_{CY}$ (mPa s) | $n_{st}$ | $\lambda_{st}$ (s) | $\lambda_{st}\dot\gamma$ | time of $n$ (s) |
|---|---|---|---|---|---|
| 0 | 72.0 | 19.3 | 0.096 | 0 | 18.8 |
| 1 | 20.2 | 5.2 | 0.026 | 0.03 | 4.3 |
| 10 | 7.5 | 1.7 | 0.0085 | 0.08 | 0.51 |
| 100 | 5.3 | 1.13 | 0.0056 | 0.56 | 0.12 |
| 300 | 5.1 | 1.06 | 0.0052 | 1.6 | 0.07 |

The structure builds up over tens of periods ($T=0.85$ s) and breaks down
within a fraction of one. The ratio of the two times is therefore not a
parameter to be swept at fixed geometry: it is set by $a(\dot\gamma)$ and spans
two decades inside one flow. In steady shear the Weissenberg number
$\lambda_{st}\dot\gamma$ exceeds one only above $\approx200$ s⁻¹; blood that
keeps the structure of rest when it is sheared reaches
$\lambda(n_{st}(0),\dot\gamma)\,\dot\gamma=0.4$ at 10 s⁻¹ and 0.75 at 100 s⁻¹.
The driver prints both tables for the constants of the run.

## Domain and openings

The default mesh is `Coarse.mesh`, the atrial volume of the PTAG atrium (mm,
fluid only), used whole. A combined Atrium mesh (fluid cells 1, solid cells 4)
also works: it is partitioned once and the owned fluid cells are extracted on
each rank, as in Atrium/Fluid. The driver reports the boundary area that
carries none of the labels below, since such faces would be traction-free
openings. The conditions are
those of Atrium/Fluid on a fixed domain: traction $-p_{pv}(t)\,n$ on the veins
(labels 1, 2, 3, 5) with their normal impedance and tangential damping,
traction $-(p_{out}+Z\,u\cdot n)\,n$ on the mitral plane (label 6), no slip on
the wall (label 4), and the backflow term on the entering flux of every
opening. With `-av_lv_model` the lumped ventricle `Rodin::Heart::CCMLC2014`
supplies $p_{out}=p_v$ and the valve law $Z=|\Gamma_m|/K$, and reads the mean
pressure of the mitral plane as its atrial pressure (`-av_lv_feedback`).

Three points differ from the arterial driver.

- The log-conformation $\psi$, the aggregate size $n$ and the age are advected
  fields on boundaries that may reverse. No value is imposed strongly. On every
  opening the term $\int|u\cdot n|_-\,(x-x_{in})\,y$ imposes the entering value:
  $\psi_{in}=0$, $n_{in}$ and age 0 on the veins, the value of the previous step
  on the mitral plane.
- The step control of the solver rejects and repeats steps. The ventricle is
  advanced once to any time the flow has not reached and is never taken back;
  a repeated step reads $p_v$ by linear interpolation in the interval already
  covered. Without rejections this is one ventricular step before each flow
  step, as in Atrium/Fluid.
- A wall that does not move lets blood through only while the valve is open,
  so the reservoir and booster phases of the atrium are absent. The comparison
  between fluids is one of rheology at equal kinematics of the wall, not a
  model of atrial filling.

## Discretisation: what was checked

The following statements were checked against the code of the driver; they are
properties of the scheme, not results of a run on the atrial mesh.

- **Logarithmic formulation.** The unknown is $\psi=\log c$ (Fattal &
  Kupferman 2004), with $\sigma=G(e^{\psi}-I)$ and $D\exp$, $D\exp^{-1}$ by
  Daleckii–Krein. The default `psi_rb` discretises the $\psi$ equation itself.
- **Stabilisation constants.** $c_1=4$, $c_2=2$, $c_3=4$, $c_4=0.25$ in
  $\alpha_1=[c_1\eta_0/h_1^2+c_2\rho|u|/h_2]^{-1}$,
  $\alpha_2=h_1^2/(c_1\alpha_1)$ and
  $\alpha_3=[c_3 f/(2\eta_p)+c_4(\lambda|u|/(2\eta_p h_2)+\lambda|\nabla u|/\eta_p)]^{-1}$:
  the values for linear elements of Castillo, Moreno, Baiges & Codina (Arch.
  Comput. Methods Eng.), which also states that they are not modified by the
  logarithmic change of variable. $\eta_0$ is the total viscosity
  $\eta_s+\eta_p/f$.
- **Owens.** The parameters are written with $(\lambda_H,\eta_{p,\infty})$ and
  the factor $f\lambda_H/\lambda$ by cell. Since $G$ is constant, this is
  $\alpha_3=[c_3/(2G\lambda)+c_4(|u|/(2Gh_2)+|\nabla u|/G)]^{-1}$ with the local
  $\lambda$, and $\tau_\psi=\alpha_3/(2G)$.
- **No stress diffusion.** The $\psi$ equation has no diffusive or
  discontinuity-capturing term. Its stabilisation is the orthogonal subscale
  of its full residual, which vanishes where the residual is resolved. The
  term-by-term variant `psi`, which adds a viscosity of order $\eta_p/2$
  wherever $\varepsilon(u)$ is not in the finite element space, is not the
  default and is not used by the campaign script.
- **Dynamic subscales** are on for velocity and log-conformation.
- **Lag.** $\lambda$ is evaluated by cell at step $n$ for both Owens fluids,
  and $\eta_{CY}$ at the extrapolated velocity for CY: a first-order splitting
  in VE and TVE under BDF2.

## Resolution and the question of turbulence

No statement about a turbulent cascade follows from a run of this driver on
the Atrium meshes without the checks below.

With $U\approx0.5$–$1$ m s⁻¹, $L\approx3$ cm and $\nu=4.6\times10^{-6}$ m² s⁻¹
the Reynolds number is $3\times10^3$–$6\times10^3$ and the Kolmogorov length
$L\,\mathrm{Re}^{-3/4}$ is below 0.1 mm. A mesh of millimetre cells is one order
coarser. The orthogonal subscales then act as an implicit subgrid model, and
the dissipation of the resolved scales is numerical in part. Constants
cannot change this; only the mesh can.

The driver measures it. At every accepted step the CSV holds the solvent and
polymer dissipation and the dissipation of the convective orthogonal subscale,
$\rho^2\int\tau\,|(\nabla u)u-\Pi[(\nabla u)u]|^2$, and the log prints their
ratio. A ratio of order one or above means that the small resolved scales are
set by the mesh. On the test mesh of `make_test_mesh.py` (mean cell size
$(V/N)^{1/3}=2.2$ mm, 40 steps per cycle) the cycle mean of the ratio is 4.7
for N, 4.3 for CY, 6.5 for VE and 5.9 for TVE, and the kinetic energy balance
of each fluid closes to 5 % of the power of the openings once this term is
counted. The driver prints the mean cell size. `psiSubscale/psiNorm` plays the
same role for the log-conformation; it is 0.17 on that mesh.

Before a difference in fluctuation level between fluids is attributed to
rheology, three things remain to be done; none was done here.

1. A second mesh and a second time step for at least N and TVE, with the
   fluctuation statistics of the probes compared across them.
2. The transition case of `examples/viscoelastic_fluids/test_kolmogorov_viscoelastic`
   (Berti & Boffetta 2010, $Wi_0=16$) with formulation `psi_rb`. It is the
   available test of whether the scheme sustains an elastic instability at a
   known threshold; no completed run of it was found in the repository.
3. Enough cycles. The phase average and the cycle-to-cycle fluctuation need
   several cycles after the initial transient; `postprocess_atrium.py` reports
   the difference between the last two cycles and marks every difference
   between fluids that is below it.

The estimate of the Weissenberg number above bears on the same question:
elasto-inertial turbulence is reported for Weissenberg numbers of order ten and
above, and the Owens constants give at most a few in the jets.

## Other limits

- **Extension.** Hookean dumbbells have an unbounded extensional stress for
  $\lambda\dot\varepsilon>1/2$. With $\lambda$ up to 0.1 s in aggregated blood
  this is reached at $\dot\varepsilon=5$ s⁻¹. The logarithm keeps the
  conformation positive, not bounded. `-av_owens_epsilon` adds the linear PTT
  factor and bounds it, but the steady shear viscosity is then
  $\eta_{pl}+G\lambda/f<\eta_{CY}$ and the three fluids no longer share it.
- **White–Metzner.** Models whose relaxation time follows the local shear
  rate are known to lose evolutionarity in some flows. The plasma viscosity
  keeps the momentum equation parabolic; no failure was observed on the test
  mesh, which is not evidence for finer ones.
- **Initial structure.** From a uniform $n$ the structure of a point at rest
  relaxes in $\approx19$ s. `-av_owens_reset_cycle k` sets $n$ at the end of
  cycle $k$ to the state of zero mean rate over that cycle,
  $n(n-1)=\langle b\,n_{st}(n_{st}-1)\rangle/\langle b\rangle$, the periodic
  state of a point whose structure is slow and which is not renewed by
  advection, and the transport continues from it. It is an initial guess: the
  age of the appendage and `nMin`, `nMax` over the following cycles show
  whether it holds.
- **Inlet state.** $\psi_{in}=0$ and a uniform $n_{in}$ ignore the shear history
  in the pulmonary veins. The entry length $U\lambda$ is 2.5 mm at
  $U=0.5$ m s⁻¹ and $\lambda=5$ ms.

## Diagnostics

Written at every accepted step to the CSV (`-av_metrics`), in SI units:
`kinetic`, `enstrophy`, `viscousDissipation` ($\int\eta_s\dot\gamma^2$),
`polymerPower` ($\int\sigma:\varepsilon$), `elasticEnergy`
($\int\mathrm{tr}\,\sigma/2$), `polymerDissipation`
($\int f\,\mathrm{tr}\,\sigma/(2\lambda)$), `ossConvective`, `ossGradDiv`,
`subscaleKinetic`, `psiNorm`, `psiSubscale`, `wiMax`, `wiVolumeFraction`
(volume with $\lambda\dot\gamma>1$), and, inside the ball
`-av_appendage_center x,y,z -av_appendage_radius r` (metres), `wiAppendage`,
`shearAppendage`, `kineticAppendage`, `ageAppendage` and the wall indices.
For Hookean dumbbells of uniform $\lambda$ the elastic energy obeys
$dE/dt=\texttt{polymerPower}-\texttt{polymerDissipation}$ up to its flux
through the openings.

`residenceTime` is the Eulerian age, $\partial_t a+u\cdot\nabla a=1$ with
$a=0$ in the entering blood. It replaces particle tracking for the question of
washout: it is a field, defined in every cell including those no particle
enters, and its appendage mean grows at rate one exactly where no blood is
renewed.

`-av_probes "x,y,z;..."` records the Gaussian-weighted mean of $u$ and
$\mathrm{tr}\,\sigma$ about each point (standard deviation `-av_probe_radius`)
at every step, for the phase average, the cycle-to-cycle fluctuation and its
spectrum in frequency. A probe filters scales below its radius.

`postprocess_atrium.py <campaign>` builds the tables and figures; with
`--fields` it adds the Q criterion, the normal stress differences in the local
shear frame and the energy spectrum of a cube inscribed in the chamber, with
the largest wavenumber the mesh resolves.

## Run

```sh
cmake --build build --target AtriumViscoelastic -j
# Test mesh with the labels of the Atrium mesh (not an atrium):
python3 examples/Heart/Atrium/Viscoelastic/make_test_mesh.py /tmp/la_test.mesh 4
mpirun -n 2 build/examples/Heart/Atrium/Viscoelastic/AtriumViscoelastic \
  -av_mesh /tmp/la_test.mesh -av_fluid owens -av_cycles 2 -av_steps_per_cycle 85 \
  -av_time_order 2 -av_appendage_center 0.034,0,0 -av_appendage_radius 0.012
# The campaign, 60 ranks per case:
MESH=/path/to/Coarse.mesh LAA_CENTER=x,y,z LAA_RADIUS=r \
  bash examples/Heart/Atrium/Viscoelastic/run_atrium_rheology.sh
```

`Coarse.mesh` and the Atrium meshes are Git LFS objects. A pointer file is rejected by the
driver and by the script.

## Status

The driver was compiled against this branch (GCC 13, PETSc 3.19) and run on the
test mesh with two ranks for N, CY, VE and TVE over two cycles with the lumped
ventricle, through `run_atrium_rheology.sh` and `postprocess_atrium.py`; no
step was rejected. It was not run on `Coarse.mesh` or the Atrium meshes, and the arterial driver's
verification cases (steady sPTT and pulsatile Oldroyd-B pipe flow) were not
repeated on this branch. The library changes of `model/viscoelasticFluid` that
cache coefficients during assembly (commit `c3e252e680`) are not on this
branch; the driver does not need them, and assembly is slower without them.
