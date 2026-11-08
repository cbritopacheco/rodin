#!/usr/bin/env python3
"""Render finite-element and topology figures from the C++ examples' CSV output."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path as _PlotPath
import sys as _plot_sys
_plot_examples = next(parent for parent in _PlotPath(__file__).resolve().parents
                      if (parent / 'rodin.mplstyle').is_file())
_plot_sys.path.insert(0, str(_plot_examples))
from plot_style import apply_style
apply_style()
from matplotlib.collections import PolyCollection
import numpy as np

PALETTE = plt.get_cmap('turbo')


def sample_triangles(samples):
    """Connect the documented per-cell triangular sampling lattice."""
    result = []
    for cell in np.unique(samples['cell']):
        rows = np.flatnonzero(samples['cell'] == cell)
        indices = {(int(samples['i'][row]), int(samples['j'][row])): row for row in rows}
        n = max(i for i, _ in indices)
        for i in range(n):
            for j in range(n - i):
                result.append([indices[i, j], indices[i + 1, j], indices[i, j + 1]])
                if i + j < n - 1:
                    result.append([indices[i + 1, j], indices[i + 1, j + 1], indices[i, j + 1]])
    return np.asarray(result)


def render(spaces, connectivity, output):
    (output / 'Images').mkdir(parents=True, exist_ok=True)
    (output / 'Examples/Geometry').mkdir(parents=True, exist_ok=True)
    basis = np.genfromtxt(spaces / 'ReferenceP1.csv', delimiter=',', names=True)
    for index in range(3):
        fig, ax = plt.subplots(figsize=(4.8, 4.3), dpi=180, layout='constrained')
        field = ax.tripcolor(basis['r'], basis['s'], basis[f'phi{index}'],
                             shading='gouraud', cmap=PALETTE, vmin=0, vmax=1)
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color='#718096', linewidth=.8)
        ax.scatter([0, 1, 0], [0, 0, 1], s=25, color='#17679A', edgecolor='#071F46', zorder=3)
        for r, t, label, offset in [(0, 0, '0', (-10, -12)), (1, 0, '1', (8, -5)), (0, 1, '2', (-10, 8))]:
            ax.annotate(label, (r, t), xytext=offset, textcoords='offset points', fontsize=9)
        ax.set(xlabel='$r$', ylabel='$s$', title=f'$\\widehat{{\\phi}}_{index}$',
               xlim=(-.05, 1.15), ylim=(-.05, 1.15), aspect='equal')
        ax.spines[['top', 'right']].set_visible(False)
        fig.colorbar(field, ax=ax, fraction=.04, pad=.03, ticks=[0, .5, 1])
        fig.savefig(output / f'Images/P1Basis-{index}.png')
        plt.close(fig)

    samples = np.genfromtxt(spaces / 'FiniteElementSpaces.csv', delimiter=',', names=True)
    triangles = sample_triangles(samples)
    cells, values = [], []
    for cell in np.unique(samples['cell']):
        rows = samples[samples['cell'] == cell]
        n = int(rows['i'].max())
        corners = [rows[(rows['i'] == i) & (rows['j'] == j)][0]
                   for i, j in [(0, 0), (n, 0), (0, n)]]
        cells.append([(r['x'], r['y']) for r in corners])
        if not np.allclose(rows['p0'], rows['p0'][0], rtol=0, atol=1e-14):
            raise ValueError('P0 export must contain one constant value per cell')
        values.append(rows['p0'][0])
    for name in ['p0', 'p1', 'p2']:
        fig, ax = plt.subplots(figsize=(4.8, 4.3), dpi=180, layout='constrained')
        if name == 'p0':
            # Original FE cells each receive their own constant, without nodal averaging.
            field = PolyCollection(cells, array=np.asarray(values), cmap=PALETTE,
                                   edgecolors='#718096', linewidths=.3)
            field.set_clim(0, 2)
            ax.add_collection(field)
        else:
            field = ax.tripcolor(samples['x'], samples['y'], triangles, samples[name],
                                 shading='gouraud', cmap=PALETTE, vmin=0, vmax=2)
        ax.set(xlabel='$x$', ylabel='$y$',
               title={'p0': r'$P_0$: cell values', 'p1': r'$P_1$', 'p2': r'$P_2$'}[name],
               xlim=(0, 1), ylim=(0, 1), aspect='equal')
        ax.spines[['top', 'right']].set_visible(False)
        fig.colorbar(field, ax=ax, fraction=.045, pad=.03).set_label(r'$I_h f$')
        fig.savefig(output / f'Images/SpaceInterpolation-{name}.png')
        plt.close(fig)

    vertices = np.genfromtxt(connectivity / 'Connectivity.vertices.csv', delimiter=',', names=True)
    edges = np.genfromtxt(connectivity / 'Connectivity.edges.csv', delimiter=',', names=True)
    positions = {int(v['id']): np.array([v['x'], v['y']]) for v in vertices}
    fig, ax = plt.subplots(figsize=(5.8, 5.1), dpi=180)
    for edge in edges:
        a, b = positions[int(edge['v0'])], positions[int(edge['v1'])]
        color = '#17679A' if edge['cells'] == 1 else '#13B6D0'
        ax.plot([a[0], b[0]], [a[1], b[1]], color=color, linewidth=2)
        middle = (a + b) / 2
        normal = np.array([-(b - a)[1], (b - a)[0]])
        normal /= np.linalg.norm(normal)
        # Place boundary labels away from the center; label the interior separately.
        if edge['cells'] == 1 and np.dot(normal, middle - .5) < 0:
            normal *= -1
        label = middle + .12 * normal
        count = int(edge['cells'])
        ax.text(*label, f"Edge {int(edge['id'])}\n{count} {'cell' if count == 1 else 'cells'}",
                ha='center', va='center', color='#17679A', fontsize=10)
    for index, point in positions.items():
        ax.scatter(*point, color='#071F46', s=35, zorder=3)
        offset = .04 * np.sign(point - .5)
        ax.text(*(point + offset), f'v{index}', ha='center', va='center')
    ax.set(xlim=(-.3, 1.3), ylim=(-.3, 1.3), aspect='equal')
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(output / 'Examples/Geometry/Connectivity.png')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spaces', type=Path, required=True)
    parser.add_argument('--connectivity', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    render(args.spaces, args.connectivity, args.output)
