"""Plot the actual native triangle fields exported by Eikonal_ReactionEikonal."""
from pathlib import Path
import argparse
import csv
import sys
import xml.etree.ElementTree as ET

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.collections import LineCollection, PolyCollection
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from plot_style import apply_style, turbo_colors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, nargs='?', default=Path('.'))
    parser.add_argument('--time', type=float, default=20, help='Voltage snapshot time in ms')
    parser.add_argument('--fibers', action='store_true', help='Preview the prescribed fiber axes')
    args = parser.parse_args()
    directory = args.directory
    data = np.genfromtxt(directory / 'ReactionEikonalActivation.csv', delimiter=',', names=True)
    with h5py.File(directory / 'ReactionEikonalVoltage.mesh.h5') as source:
        topology = source['/Mesh/XDMF/Topology'][:].reshape(-1, 4)
    if not np.all(topology[:, 0] == 4):
        raise ValueError('Expected native triangles')
    triangulation = mtri.Triangulation(data['x_mm'], data['y_mm'], topology[:, 1:])
    if args.fibers:
        with (directory / 'ReactionEikonalParameters.csv').open() as source:
            parameters = {row['parameter']: float(row['value']) for row in csv.DictReader(source)}
        if parameters['axis_aligned_geometry'] != 1:
            raise ValueError('Expected the axis-aligned uniform-grid example')
        fiber = np.array([parameters['fiber_x'], parameters['fiber_y']])
        points = np.column_stack([data['x_mm'], data['y_mm']])
        x_indices = np.searchsorted(np.unique(points[:, 0]), points[:, 0])
        y_indices = np.searchsorted(np.unique(points[:, 1]), points[:, 1])
        n = int(parameters['intervals'])
        mask = ((x_indices % 4 == 0) & (y_indices % 4 == 0)
                & (x_indices > 0) & (x_indices < n) & (y_indices > 0) & (y_indices < n))
        centers = points[mask]
        segments = np.stack([centers - 0.3 * fiber, centers + 0.3 * fiber], axis=1)
        apply_style()
        figure, axis = plt.subplots(figsize=(7, 7), layout='constrained')
        axis.add_collection(PolyCollection(points[topology[:, 1:]], facecolors='#eef2f7',
                                           edgecolors='#c5d0dc', linewidths=0.2))
        axis.add_collection(LineCollection(segments, colors=[turbo_colors(2)[1]],
                                           linewidths=2.5, capstyle='round'))
        axis.scatter(parameters['source_x_mm'], parameters['source_y_mm'], marker='*',
                     s=100, facecolor='white', edgecolor='#303b49', zorder=5)
        angle = np.degrees(np.arctan2(fiber[1], fiber[0]))
        axis.set(xlim=(-0.3, 10.3), ylim=(-0.3, 10.3), aspect='equal',
                 xlabel=r'$x\;[\mathrm{mm}]$', ylabel=r'$y\;[\mathrm{mm}]$',
                 title='Prescribed fibers on the axis-aligned uniform grid')
        axis.text(0.035, 0.97,
                  rf'$\mathbf{{f}}=({fiber[0]:g},\,{fiber[1]:g})^T,\quad\theta={angle:.2f}^\circ$',
                  transform=axis.transAxes, va='top', fontsize=13,
                  bbox=dict(facecolor='white', edgecolor='none', alpha=0.9))
        figure.savefig(directory / 'ReactionEikonalFibers.png', dpi=180, transparent=True)
        figure.savefig(directory / 'ReactionEikonalFibers.svg', transparent=True)
        plt.close(figure)
        return
    frames = ET.parse(directory / 'ReactionEikonalVoltage.xdmf').findall('.//Grid[@GridType="Uniform"]')
    frame = min(frames, key=lambda grid: abs(float(grid.find('Time').get('Value')) - args.time))
    time = float(frame.find('Time').get('Value'))
    reference = frame.find('Attribute[@Name="voltage_mV"]/DataItem').text.strip()
    filename, dataset = reference.split(':', 1)
    with h5py.File(directory / filename) as source:
        voltage = source[dataset][:]
    apply_style()
    figure, axes = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')
    maximum = max(data['fim_foot_onset_ms'].max(), data['replus_lat_ms'].max())
    panels = [
        (data['fim_foot_onset_ms'], r'FIM foot-current onset $T$', r'$T\;[\mathrm{ms}]$', 0, maximum),
        (data['replus_lat_ms'], 'RE+ steepest-upstroke activation', r'$t_{\mathrm{LAT}}\;[\mathrm{ms}]$', 0, maximum),
        (voltage, rf'RE+ voltage at $t={time:g}\,\mathrm{{ms}}$', r'$V\;[\mathrm{mV}]$', -85, 35),
        (data['lat_minus_foot_onset_ms'], 'Onset-to-upstroke lag', r'$t_{\mathrm{LAT}}-T\;[\mathrm{ms}]$', None, None),
    ]
    for axis, (values, title, label, lower, upper) in zip(axes.flat, panels):
        image = axis.tripcolor(triangulation, values, shading='gouraud', cmap='turbo', vmin=lower, vmax=upper)
        axis.set(title=title, xlabel=r'$x\;[\mathrm{mm}]$', ylabel=r'$y\;[\mathrm{mm}]$', aspect='equal')
        figure.colorbar(image, ax=axis, label=label, shrink=0.8)
    figure.savefig(directory / 'ReactionEikonal.png', dpi=180)
    figure.savefig(directory / 'ReactionEikonal.svg')
    plt.close(figure)


if __name__ == '__main__':
    main()
