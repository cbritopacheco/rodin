# Branching truncation benchmark

The full Y resolves the parent, junction, and both straight daughter
extensions. The cut Y retains the same cells up to the cuts. The two outlets
independently replace the removed extensions using the signed pressure law
`p_cut = R(Q) Q`. A steady solution on the full Y supplies the reference
flow and cross-sectional mean pressure at each cut. Both outlet laws match
these values, and the runner uses the same calibration at both pressure heads.

The parent radius is 100 µm in 3D, and the daughter radii are 72 and 46 µm.
The daughters turn by +30° and −38°. The fully resolved distal ends have
zero applied pressure. The inlet receives a sinusoidal pressure of period
0.85 s and amplitude 45% of its mean. The default mean heads are 8.2 and
1.64 Pa. The reported narrow-branch share is its signed cycle volume divided
by the sum of both branch volumes. The comparison also reports the RMS error
of the instantaneous split in percentage points and the RMS branch-flow
errors normalized by the full model's mean absolute branch flows. Since the
extension walls are rigid, the integrated full-model flow at each distal
end equals the flow at its cut for the incompressible reference.

## Build and run

Use PETSc with MUMPS and MPI, and install NumPy and Matplotlib.
Mesh regeneration also requires Gmsh. From the Rodin repository root:

```bash
git submodule update --init third-party/mmg third-party/termcolor
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DRODIN_USE_MPI=ON -DRODIN_USE_PETSC=ON -DRODIN_BUILD_EXAMPLES=ON \
  -DRODIN_INSTALL_RESOURCES=OFF
cmake --build build --target TruncationY -j 6
solver="$(pwd)/build/examples/Heart/CoronaryArtery/TruncationY/TruncationY"
cd examples/Heart/CoronaryArtery/TruncationY
python3 run_benchmark.py --solver "$solver" --mesh data/meshes/3d-coarse \
  --output /tmp/results-y-3d --steps 40
```

The archived waveforms were computed with Rodin revision
`e2674b4ffa2e5b26b00d7643c0f1278f1e9caf20`. Record the Rodin revision for
any new runs. The archived meshes and time series are in `data/`.
Use `python3 report_refinement.py --data data/truncation_y_3d.csv` to print
the 3D refinement results without plotting. Pass the planar CSV to the
same command to report the planar results.

The campaign writes one CSV per case, a calibration file, `summary.csv`,
and `flow_split.png`. Each case uses one local MPI process. Double `--steps`
for temporal refinement. For mesh refinement, regenerate the mesh at a
smaller `--mesh-size`. The runner's `--calibration` option carries forward a
given calibration file without fitting the outlets again.

## 3D refinement and manuscript figure

The manuscript uses Gmsh 4.15.2 with mesh sizes 0.46 and 0.34 in units of
the parent radius. The underlying full meshes have 3031 and 5231 tetrahedra.
On each grid the truncated mesh is a subset of the full mesh. Use the
archived coarse-grid calibration for both fine-grid campaigns. Each campaign
runs the full Y and both outlet laws at each pressure head for two cycles.

The frozen wide and narrow reference flows are, respectively,
$1.8652251\times10^{-11}$ and $2.2258514\times10^{-12}\,\mathrm{m^3/s}$.
The corresponding mean cut pressures are $3.8175956$ and $4.5499727\,\mathrm{Pa}$.

```bash
python3 ./run_benchmark.py --solver "$solver" --mesh data/meshes/3d-coarse --output baseline40 --steps 40
python3 ./run_benchmark.py --solver "$solver" --mesh data/truncation_y_refined --output fine40 --steps 40 --calibration data/truncation_y_refined/calibration.txt
python3 ./run_benchmark.py --solver "$solver" --mesh data/truncation_y_refined --output fine80 --steps 80 --calibration data/truncation_y_refined/calibration.txt
pvpython --no-mpi --disable-registry ./sample_velocity.py --field fine80/reduced_head/full.csv.vtk --mesh data/truncation_y_refined/full.mesh --trace fine80/reduced_head/full.csv --volume peak-volume.npz --vtu peak-recomputed.vtu
python3 ./report_refinement.py --runs baseline40 fine40 fine80 --data reproduced_3d.csv --mesh data/truncation_y_refined/full.mesh --volume peak-volume.npz --diagnostics fig_y_diagnostics.pdf
```

To regenerate Figure 8 from its three peak-pressure fields and the
planar and circular time series, run the following from this directory:

```bash
pvpython --no-mpi --disable-registry ./sample_comparison.py --mesh-dir data/truncation_y_refined --trace data/truncation_y_3d.csv --full data/truncation_y_low_head_peak.vtu --variable data/truncation_y_low_head_peak_variable.vtu --constant data/truncation_y_low_head_peak_constant.vtu --output peak-comparison.npz
python3 ./report_refinement.py --data data/truncation_y_3d.csv --planar-data data/truncation_y_2d.csv --comparison peak-comparison.npz --figure fig_y_flow_split.pdf
```

All three fields are sampled at maximum inlet pressure in the second cycle
under the lower driving pressure. The sampler checks the three cap fluxes
against each case's trace. It then transfers the coincident quadratic
velocity nodes to the exact same retained tetrahedra and interpolates the
three fields at identical volume positions. The views share the camera,
speed scale and speed-dependent opacity. The narrow-branch flow below
each rendering is its instantaneous value at this phase.

The lower panels show the change in the narrow daughter's share of second-cycle
volume on the finest grid for each geometry. Both use the same scale in
percentage points. A negative value means the branch receives a smaller share
at the lower pressure. Compare each outlet with the fully resolved row.

To regenerate the mesh, resolved velocity field and instantaneous
waveforms in Figure 9, use `--field data/truncation_y_low_head_peak.vtu --trace data/truncation_y_3d.csv`
in the sampling command, omit `--vtu`, and pass the resulting volume and
finest mesh with `--diagnostics fig_y_diagnostics.pdf` to the
plotting command. The VTK bindings used for sampling are available in
ParaView's `pvpython`. The mesh SHA-256 digest must match the mesh recorded
with the waveforms. The refinement table reports each truncated model's
share minus the resolved share at the same pressure and discretisation.

All three low-load 3D runs can export the quadratic velocity field at peak
inlet pressure in the second cycle (step 100 at 80 steps per cycle). The
reference field is archived in `data/truncation_y_low_head_peak.vtu` and
the truncated fields in `data/truncation_y_low_head_peak_constant.vtu` and
`data/truncation_y_low_head_peak_variable.vtu`. Each file contains the
velocity at all ten quadratic nodes of every tetrahedron in SI units.
The sampler checks the tetrahedra and integrates normal velocity over the
inlet and both terminal faces of each mesh against its recorded flows.
The renderer interpolates the velocity vector inside the tetrahedra, takes
its magnitude, and composites colour and opacity along each viewing ray.
Surface edges and the silhouette come from the refined mesh.

## Planar version

The reported planar runs are archived in `data/truncation_y_2d.csv`. Each
row records one accepted step, including the two outlet fluxes per unit
depth in square metres per second. The full meshes contain 407 and 739
triangles, with 40 and 80 steps per cycle, respectively. Both pressure
heads and both discretisations use the planar coarse-mesh calibration.
Compute the narrow-branch share by summing its flux over the second cycle
and dividing by the summed fluxes from both daughters over the same steps.

Run `run_benchmark.py --dim 2` with `data/meshes/2d-coarse` and 40 steps
per cycle, then with `data/meshes/2d-fine` and 80 steps per cycle. Use
separate output directories and pass the coarse calibration file to the
fine run with `--calibration`. To regenerate a planar mesh, pass `--dim 2`
to `make_mesh.py` as well. The same layout gives channels with half-widths equal to the
3D radii. The planar solver uses flow per unit depth and slit hydraulic
viscosity. For half-width $H$ and flow per unit depth $Q'$, the apparent
shear rate is $3|Q'|/(2H^2)$. The slit integration uses the same
Carreau–Yasuda material parameters as the circular-tube table. Compare the
2D and 3D flow-split errors against their respective fully resolved Y cases.
