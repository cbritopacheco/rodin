# Coronary outlet comparisons

These scripts supplement the [CMBPE campaign](../README.md). The two
post-processors take output from a resolved coronary calculation. The 0D
replay carries its own pressure input and runs independently of that solver.

## 0D pressure replay

Run from this directory with Python, NumPy, SciPy (1.12 or later) and
Matplotlib:

```sh
python3 compare_terminal_0d.py --output-dir /tmp/cmbpe-terminal-0d
```

`replay_pressures.csv` contains the time, aortic pressure and intramyocardial
pressure columns of the simulated healthy shear-dependent coronary waveform
`data/New data/Nor_700_variable.csv` from manuscript revision
`71c5bc319504516ef0268129dff7c873b000f6d2`. The script repeats the third cycle
of these pressures, then reduces the driving head while holding the outlet
calibration fixed. It reports the cycle mean and waveform for a
shear-dependent outlet and two ways of setting a constant resistance.
The 0D calculation isolates the effect of the outlet law under prescribed
pressures.

## Resolved coronary comparisons

Pass two completed run directories to compare cycle-3 inflow, phase
volumes, mean wall shear stress and dissipation. The script reads
`CoronaryArtery.csv` and `CoronaryArtery_fields.csv` from each run:

```sh
python3 postprocess.py /path/to/variable-run /path/to/constant-run \
  --out comparison.json --waveform comparison-waveform.csv
```

To compare wall shear stress at corresponding nodes, use two XDMF/HDF5
outputs with the same mesh and snapshot times. Export the nodal
`shearStress` vector and use the MEDIT mesh supplied to the simulations.
The script requires NumPy, Matplotlib, and either h5py or HDF5's `h5dump`.

```sh
python3 wall_map.py --reference /path/to/variable.xdmf \
  --matched /path/to/constant.xdmf --mesh ../meshes/coronary_medium.mesh \
  --out wall-shear-comparison
```

The wall map gives the difference in cycle-averaged wall shear stress
magnitude. The constant-minus-variable convention matches the scalar
comparison. The scripts require the run outputs as inputs.
