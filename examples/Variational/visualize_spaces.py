#!/usr/bin/env python3
"""Render finite-element and topology figures from the C++ examples' CSV output."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path as _PlotPath
import sys as _plot_sys
from pathlib import Path as _PlotPath
_plot_examples = next(parent for parent in _PlotPath(__file__).resolve().parents
                      if (parent / 'rodin.mplstyle').is_file())
_plot_sys.path.insert(0, str(_plot_examples))
from plot_style import apply_style
apply_style()
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

PALETTE = LinearSegmentedColormap.from_list(
    'Rodin', ['#071F46', '#17679A', '#13B6D0', '#F4F8FB'])


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
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.5), dpi=180)
    for index, ax in enumerate(axes):
        field = ax.tripcolor(basis['r'], basis['s'], basis[f'phi{index}'],
                             shading='gouraud', cmap=PALETTE, vmin=0, vmax=1)
        ax.plot([0, 1, 0, 0], [0, 0, 1, 0], color='#718096', linewidth=.8)
        ax.scatter([0, 1, 0], [0, 0, 1], s=25, color='#17679A', edgecolor='#071F46', zorder=3)
        for r, s, label, offset in [(0, 0, '0', (-10, -12)), (1, 0, '1', (8, -5)), (0, 1, '2', (-10, 8))]:
            ax.annotate(label, (r, s), xytext=offset, textcoords='offset points', fontsize=9)
        ax.set(xlabel='$r$', ylabel='$s$', title=f'$\\widehat{{\\phi}}_{index}$',
               xlim=(-.05, 1.15), ylim=(-.05, 1.15), aspect='equal')
        ax.spines[['top', 'right']].set_visible(False)
        fig.colorbar(field, ax=ax, fraction=.04, pad=.03, ticks=[0, .5, 1])
    fig.tight_layout()
    fig.savefig(output / 'Images/P1Basis.png')
    plt.close(fig)

    samples = np.genfromtxt(spaces / 'FiniteElementSpaces.csv', delimiter=',', names=True)
    triangles = sample_triangles(samples)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.8), dpi=180)
    for ax, name in zip(axes, ['p0', 'p1', 'p2']):
        error = np.abs(samples[name] - samples['exact'])
        field = ax.tripcolor(samples['x'], samples['y'], triangles, error,
                             shading='gouraud', cmap=PALETTE,
                             vmin=0, vmax=max(error.max(), 3e-15))
        ax.set(xlabel='$x$', ylabel='$y$', title={'p0': 'P0', 'p1': 'P1', 'p2': 'H1 degree 2'}[name],
               aspect='equal')
        bar = fig.colorbar(field, ax=ax, fraction=.045, pad=.03)
        bar.ax.ticklabel_format(style='sci', scilimits=(-2, 2))
    fig.suptitle('Absolute interpolation error · same mesh and quadratic data', y=.99)
    fig.tight_layout(rect=(0, 0, 1, .91))
    fig.savefig(output / 'Images/SpaceInterpolation.png')
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
