# Coronary FSI coupled to the reduced poroelastic ventricle

`CoronaryArtery_Explicit_PoroElastic` couples the coronary ALE FSI to
`Rodin::Heart::PoroelasticSphereT` (Ruz et al., *A multiscale formulation for
non-Newtonian coronary FSI coupled to reduced myocardial poromechanics*).

## Reduced ventricle (Sec. 3)

| Paper | Code |
|---|---|
| (20)–(22) shell kinematics, `h = h0 λr` | `PoroelasticSphere/Physics/Kinematics.h` |
| (4), (27)–(30) skeletal stress, `λ̄`, pressure law, `p̃` | `Physics/WallStress.h`, `Physics/Perfusion.h` |
| (33)–(34), (B.3)–(B.6) Usyk energy | `UsykSphericalEnergy` (`PassiveLaw.h`) |
| (35)–(38) exponential storage | `PerfusionEvaluator::storage` |
| (31), (40)–(41) blood-volume balance, drainage, pore shear | `PerfusionEvaluator::evaluate`, `DynamicSystem.h` |

Viscous skeleton: `S^v = η Ė`, i.e. `P^v_RR = η λr² λ̇r`, `P^v_θθ = η λθ² λ̇θ`.
Active: in-plane average of the CCMLC2014 fiber stress, `P^act_θθ = σ1D λθ / 2`.
The Jacobian is exact (checked entrywise against finite differences in
`tests/unit/Rodin/Heart/PoroelasticSphereTest.cpp`).

## Interface (Sec. 4.3, 5.2)

Robin–Robin on the conforming interface, transferred node by node:

```
solid : P n_s + α v_Γ = α u_f + T          v_Γ = (d^{n+1} − d^n)/δt
fluid : J_Γ σ n_f + α u_f = α v_Γ − T      (v_Γ = ALE mesh velocity on Γ)
T ← T + α (u_f − v_Γ)                      (variational fluid reaction)
```

`T` is the traction per reference area and `J_Γ = dS/dS_0`. Both sides use the
same BDF1 interface rate, so the fixed point has `u_f = v_Γ` and
`P n_s = −J_Γ σ n_f`. The recovered velocity gradient is only used for WSS.

## Scenarios

`-coronary_scenario` (shared with `PoroelasticSphere_0D <scenario>`):

| name | change | path to the 3D model |
|---|---|---|
| `healthy` | reference | — |
| `aortic_stenosis` | aortic valve conductance ×0.2 | higher LV pressure → higher `p̃` (pv/2 term) and lower `p_ar` at the ostium |
| `hyperviscosity` | blood viscosity ×1.6 (both plateaus) | 3D stresses/WSS, outlet χ, pore drainage `ν_eff`, Windkessel |
| `lvh` | `h0` ×1.35, `c0` ×2.5, `C0` ×0.6, `K` ×0.5 | stiffer, thicker wall: `p̃`, storage and attachment motion |

The distal network (`R_a`, `R_v`, `C`) is calibrated on the healthy hand-off
for every scenario, so only the scenario differs between runs.

## Geometries

`meshes/` holds the two FSI pairs as plain gzip files (no Git LFS):

| key | fluid / solid | notes |
|---|---|---|
| `nostenosis` | `NoStenosisFluid` / `NoStenosisSolid` | built by `meshes/make_nostenosis_mesh.py` |
| `stenosis80` | `Stenosis80Fluid` / `Stenosis80Solid` | = `Stenosis80.mesh` / `Prismatic80.mesh` |

The stenotic pair is the healthy `FSIFluid`/`FSISolid` pair scaled by 25/66 and
deformed at the lesion. The stenosis-free solid is the stenotic solid with only
its 832 lesion vertices replaced (same wedges, labels, wall thickness and
contact patch); away from the lesion both pairs agree to 1e-7 mm.

```sh
CORES_PER_CASE=50 CYCLES=4 examples/Heart/CoronaryArtery/run_poroelastic_scenarios.sh
# batch 1: nostenosis x {healthy, aortic_stenosis, hyperviscosity}, 3 x 50 cores
# batch 2: stenosis80 x {...} when batch 1 ends; results/coronary_poroelastic/<geometry>/<scenario>
```

## Options added

`-coronary_scenario`, `-coronary_results_dir`, `-coronary_fluid_mesh`,
`-coronary_solid_mesh`, `-coronary_lumped_supply` (0D lumped supply
conductance, default 0), `-coronary_heart_disp_epicardial` (attachment target
with `y + (h − h0)/2` instead of `y`).
