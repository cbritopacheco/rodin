# PDE tutorials

The standalone `Poisson`, `Elasticity`, `Helmholtz`, `Periodic`, `SurfaceMesh`,
`Stokes`, and `NavierStokes` targets generate their own meshes and write XDMF
files with accompanying HDF5 data in the working directory. Enable examples
with `-DRODIN_BUILD_EXAMPLES=ON`. Run numerical examples one at a time.

To regenerate the tutorial figures and cavity video, run each executable from
its own directory named after the target, under a common output directory:

```sh
mkdir -p output/Helmholtz
cd output/Helmholtz
/path/to/build/examples/PDEs/Helmholtz
```

Repeat for `Periodic`, `SurfaceMesh`, `Stokes`, and `NavierStokes`. Then, from
this repository's root, run:

```sh
python3 examples/PDEs/visualize.py /path/to/output /path/to/figures
```

The plotting script requires NumPy, h5py, Matplotlib, and FFmpeg with an H.264
encoder. It reads mesh coordinates and evaluated point values from the XDMF
references, rather than interpreting finite element coefficient ordering.
Streamlines use linear interpolation of those exported samples. They are
instantaneous streamlines, not transient particle trajectories.

The cavity uses a corner-compatible polynomial lid at Reynolds number 100,
Taylor–Hood elements, and one semi-implicit Oseen solve per time step. Its
animation covers simulation time 0–10; it is not a steady-state or mesh-converged
benchmark. These small tutorial solves do not change the benchmark workloads.
