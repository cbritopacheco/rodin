# Ectatic tube thrombosis

An ectatic tube introduces regions of slow recirculating flow in which washout
and coagulation can be studied. A sequential flow solve supplies the velocity
for one-way transport of thrombin, fibrinogen, fibrin, and a washout tracer.
The outlet is traction-free; no lumped outlet circulation is imposed.

For tracer concentration $c$ and fluid velocity $u$, the washout model is

$$
\partial_t c+u\cdot\nabla c=0,
\qquad c(X,0)=1,\qquad c|_{\Gamma_{\mathrm{in}}}=0.
$$

The volume average of $c$ measures the remaining fraction of the initial bolus.
Coagulation species use the reaction and wall-source terms defined by the
[thrombosis model](../CoronaryArtery/ThrombosisModel.h). Species transport does
not alter the flow or deform the wall.

Requires `RODIN_USE_PETSC=ON`; MPI is not required.

## Build and run

From the repository root, using a configured build tree:

```sh
cmake --build build --target EctaticTubeThrombosis -j4
cd build
./examples/Heart/EctaticTubeThrombosis/EctaticTubeThrombosis
```

Run from the build root so the default relative resource paths resolve.
Output files are written relative to the working directory. Hydrate mesh
resources with `git lfs pull` before running.

## Input and configuration

The default mesh is `resources/examples/Heart/EctaticTubeThrombosis.mesh`. Generate it from the repository root with:

```sh
python3 examples/Heart/EctaticTubeThrombosis/make_test_mesh.py resources/examples/Heart/EctaticTubeThrombosis.mesh
```

The mesh can be selected with `-pt_mesh`. The `-pt_flow_cycles` option controls
flow warm-up, `-pt_species_cycles` controls the subsequent transport interval,
and `-pt_dt` selects the time step in seconds.
