"""Plot the actual native triangle fields exported by Eikonal_ReactionEikonal."""
from pathlib import Path
import argparse
import sys
import xml.etree.ElementTree as ET

import h5py
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from plot_style import apply_style


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path, nargs='?', default=Path('.'))
    parser.add_argument('--time', type=float, default=20, help='Voltage snapshot time in ms')
    args = parser.parse_args()
    directory = args.directory
    data = np.genfromtxt(directory / 'ReactionEikonalActivation.csv', delimiter=',', names=True)
    with h5py.File(directory / 'ReactionEikonalVoltage.mesh.h5') as source:
        topology = source['/Mesh/XDMF/Topology'][:].reshape(-1, 4)
    if not np.all(topology[:, 0] == 4):
        raise ValueError('Expected native triangles')
    triangulation = mtri.Triangulation(data['x_mm'], data['y_mm'], topology[:, 1:])
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
