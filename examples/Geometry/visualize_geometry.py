#!/usr/bin/env python3
"""Plot reference polytopes and mapped quadrature from Quadrature.cpp CSV output."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from mpl_toolkits.mplot3d import proj3d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    import sys as _plot_sys
    from pathlib import Path as _PlotPath
    _plot_examples = next(parent for parent in _PlotPath(__file__).resolve().parents
                          if (parent / 'rodin.mplstyle').is_file())
    _plot_sys.path.insert(0, str(_plot_examples))
    from plot_style import apply_style, turbo_colors
    apply_style()
    args.output.mkdir(parents=True, exist_ok=True)
    data = np.genfromtxt(args.root / 'Quadrature.csv', names=True, delimiter=',')
    triangles = [[[0, 0], [1, 0], [0, 1]], [[1, 1], [3, 1], [1, 4]]]
    for name, coordinates, triangle, title in zip(['reference', 'physical'], [('r', 's'), ('x', 'y')],
            triangles, ['Reference triangle', 'Physical triangle']):
        fig, ax = plt.subplots(figsize=(4.8, 4.5), layout='constrained')
        ax.add_patch(Polygon(triangle, closed=True, facecolor='#166699',
                             alpha=.09, edgecolor='#718096', linewidth=1.5))
        ax.scatter(data[coordinates[0]], data[coordinates[1]],
                   s=data['w'] / max(data['w']) * 200, color=turbo_colors(1)[0])
        ax.set_aspect('equal')
        ax.autoscale_view()
        ax.margins(.15)
        ax.set_xlabel('$'+coordinates[0]+'$')
        ax.set_ylabel('$'+coordinates[1]+'$', rotation=0)
        ax.set_title(title, pad=14)
        ax.spines[['top', 'right']].set_visible(False)
        fig.savefig(args.output / f'Quadrature-{name}.png', dpi=180, transparent=True)
        plt.close(fig)

    # A zero-dimensional reference point has no coordinates. Its schematic
    # places the unique vertex at the center of the panel and labels its index.
    fig, ax = plt.subplots(figsize=(4.8, 4.3), layout='constrained')
    ax.scatter([0], [0], color='#13b6d0', edgecolors='white', linewidths=1,
               s=48, zorder=5)
    ax.annotate('0', (0, 0), xytext=(0, 25), textcoords='offset points',
                ha='center', va='center', fontsize=10, color='white',
                bbox=dict(boxstyle='circle,pad=.22', facecolor='#123657',
                          edgecolor='white', linewidth=.7))
    ax.set(xlim=(-.7, .7), ylim=(-.7, .7), aspect='equal')
    ax.set_axis_off()
    ax.set_title('Point', fontsize=13, pad=3)
    fig.savefig(args.output / 'ReferencePolytopes-point.png', dpi=180, transparent=True)
    plt.close(fig)

    vertices = np.genfromtxt(args.root / 'ReferencePolytopes.csv', names=True, delimiter=',')
    names = ['Segment', 'Triangle', 'Quadrilateral', 'Tetrahedron', 'Pyramid', 'Hexahedron', 'Wedge']
    edges = [ [(0,1)], [(0,1),(1,2),(2,0)], [(0,1),(1,2),(2,3),(3,0)],
              [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],
              [(0,1),(1,2),(2,3),(3,0),(0,4),(1,4),(2,4),(3,4)],
              [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)],
              [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)] ]
    for i, name in enumerate(names):
        selected = vertices[vertices['type'] == i+1]
        xyz = np.stack([selected['x'], selected['y'], selected['z']], axis=1)
        dim = int(selected['dimension'][0])
        fig = plt.figure(figsize=(4.8, 4.3), layout='constrained')
        ax = fig.add_subplot(projection='3d' if dim == 3 else None)
        if dim == 3:
            ax.add_collection3d(Line3DCollection([xyz[list(edge)] for edge in edges[i]],
                                               colors='#166699', linewidths=1.7))
            ax.scatter(*xyz.T, color='#13b6d0', edgecolors='white', linewidths=1,
                       s=48, depthshade=False)
            ax.set(xlim=(-.1,1.1), ylim=(-.1,1.1), zlim=(-.1,1.1))
            ax.set_box_aspect((1,1,1)); ax.view_init(22, -58)
        else:
            for edge in edges[i]:
                ax.plot(*xyz[list(edge), :2].T, color='#166699', linewidth=1.7)
            ax.scatter(*xyz[:,:2].T, color='#13b6d0', edgecolors='white', linewidths=1,
                       s=48, zorder=5)
            ax.set(xlim=(-.2,1.2), ylim=(-.2,1.2)); ax.set_aspect('equal')
        if dim == 1:
            ax.set_ylim(-.7, .7)
        ax.set_axis_off()
        ax.set_title(name, fontsize=13, pad=3)
        # Project the reference vertices before placing their local indices.
        # Outward offsets keep the badges clear of the incident edges.
        if dim == 3:
            px, py, _ = proj3d.proj_transform(*xyz.T, ax.get_proj())
            projected = np.column_stack((px, py))
        else:
            projected = xyz[:, :2]
        center = projected.mean(axis=0)
        for vertex, position in zip(selected['vertex'], projected):
            direction = position - center
            length = np.linalg.norm(direction)
            offset = 25 * direction / length if length > 1e-12 else np.array([0, 25])
            # The two interior projected cube vertices have edges close to
            # their radial directions; use the clear sectors above/below them.
            if name == 'Hexahedron' and int(vertex) in (3, 5):
                offset = np.array([0, -30 if int(vertex) == 3 else 25])
            if name == 'Wedge' and int(vertex) == 2:
                offset = np.array([-25, 15])
            ax.annotate(str(int(vertex)), position, xytext=offset,
                        textcoords='offset points', ha='center', va='center',
                        fontsize=10, color='white', clip_on=False,
                        bbox=dict(boxstyle='circle,pad=.22', facecolor='#123657',
                                  edgecolor='white', linewidth=.7))
        fig.savefig(args.output / f'ReferencePolytopes-{name.lower()}.png', dpi=180, transparent=True)
        plt.close(fig)


if __name__ == '__main__':
    main()
