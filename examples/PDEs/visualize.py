#!/usr/bin/env python3
"""Render tutorial figures from the examples' XDMF/HDF5 output.

Requires NumPy, h5py, Matplotlib, and FFmpeg for the cavity video.
Usage: python visualize.py OUTPUT_ROOT ASSET_DIRECTORY
Each executable is run in its own OUTPUT_ROOT/<executable> directory.
"""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

import h5py
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
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.tri import LinearTriInterpolator, Triangulation
import numpy as np

SEQUENTIAL = LinearSegmentedColormap.from_list(
    'Rodin', ['#071F46', '#17679A', '#13B6D0', '#F4F8FB'])
DIVERGING = LinearSegmentedColormap.from_list(
    'Rodin signed', ['#071F46', '#F4F8FB', '#13B6D0'])


def snapshot(root, name, index=-1):
    """Read exported point values; never infer a finite element DOF layout."""
    path = root / name / (name + '.xdmf')
    grids = ET.parse(path).findall('.//Grid[@GridType="Uniform"]')
    grid = grids[index]

    def data(element):
        filename, dataset = element.text.strip().split(':', 1)
        with h5py.File(path.parent / filename) as handle:
            return handle[dataset][:]

    vertices = data(grid.find('Geometry/DataItem'))
    topology = data(grid.find('Topology/DataItem')).reshape(-1, 4)
    if not np.all(topology[:, 0] == 4):
        raise ValueError('These tutorial figures require triangular output meshes')
    fields = {a.get('Name'): data(a.find('DataItem'))
              for a in grid.findall('Attribute') if a.get('Center') == 'Node'}
    if not all(np.isfinite(field).all() for field in fields.values()):
        raise ValueError('Non-finite exported field')
    return vertices, topology[:, 1:].astype(int), fields, float(grid.find('Time').get('Value'))


def planar(ax, xy, triangles, values, title, limit=1):
    field = ax.tripcolor(*xy.T, triangles, values, shading='gouraud',
                         cmap=DIVERGING, vmin=-limit, vmax=limit)
    ax.set_aspect('equal')
    ax.set(xlabel='$x$', ylabel='$y$', title=title)
    return field


def flow_figure(root, name, index=-1):
    xy, triangles, fields, time = snapshot(root, name, index)
    velocity = fields['velocity']
    tri = Triangulation(*xy.T, triangles)
    coordinates = np.linspace(0, 1, 121)
    xx, yy = np.meshgrid(coordinates, coordinates)
    ux = LinearTriInterpolator(tri, velocity[:, 0])(xx, yy)
    uy = LinearTriInterpolator(tri, velocity[:, 1])(xx, yy)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.8), dpi=140)
    speed = np.linalg.norm(velocity, axis=1)
    image = axes[0].tripcolor(tri, speed, shading='gouraud', cmap=SEQUENTIAL, vmin=0,
                                vmax=2 if name == 'NavierStokes' else 1)
    axes[0].set_title('Speed', pad=36)
    axes[1].streamplot(coordinates, coordinates, ux, uy, density=1.8,
                       color='#17679A', linewidth=.8, arrowsize=.8)
    axes[1].set_title('Instantaneous streamlines', pad=36)
    for ax in axes:
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='$x$', ylabel='$y$', aspect='equal')
        if name == 'NavierStokes':
            # The trace is the computed boundary velocity, not an illustrative lid.
            lid = np.isclose(xy[:, 1], 1)
            order = np.argsort(xy[lid, 0])
            ax.plot(xy[lid, 0][order], 1.025 + .035 * velocity[lid, 0][order],
                    transform=ax.get_xaxis_transform(), color='#13B6D0',
                    linewidth=2, clip_on=False)
            center = .5 - .22 * np.cos(2 * np.pi * time / 8)
            ax.annotate('', xy=(center + .08, 1.12), xytext=(center - .08, 1.12),
                        xycoords='axes fraction',
                        arrowprops={'arrowstyle': '->', 'color': '#17679A', 'lw': 1.4})
        else:
            ax.annotate('', xy=(.82, 1.035), xytext=(.18, 1.035), xycoords='axes fraction',
                        arrowprops={'arrowstyle': '->', 'color': '#17679A', 'lw': 1.4})
    fig.colorbar(image, ax=axes[0], fraction=.045, pad=.03)
    if name == 'NavierStokes':
        fig.suptitle(rf'Traveling lid forcing · $\mathrm{{Re}}=200$ · $t={time:.2f}$', y=.99)
    fig.tight_layout(pad=1.2, rect=(0, 0, 1, .93 if name == 'NavierStokes' else 1))
    return fig


def render_mixed_boundary(root, output):
    """Plot the computed mixed-boundary solution and identify the boundary data."""
    output.mkdir(parents=True, exist_ok=True)
    xy, triangles, fields, _ = snapshot(root, 'MixedBoundary')
    fig, ax = plt.subplots(figsize=(8, 5), dpi=180)
    image = ax.tripcolor(*xy.T, triangles, fields['u'], shading='gouraud',
                         cmap=SEQUENTIAL, vmin=0, vmax=1)
    ax.plot([0, 0], [0, 1], color='#071F46', linewidth=3)
    ax.plot([1, 1], [0, 1], color='#13B6D0', linewidth=3)
    ax.plot([0, 1], [0, 0], color='#718096', linewidth=1)
    ax.plot([0, 1], [1, 1], color='#718096', linewidth=1)
    for y in [.25, .5, .75]:
        ax.annotate('', xy=(1.18, y), xytext=(1.02, y),
                    arrowprops={'arrowstyle': '->', 'color': '#17679A'})
    ax.text(-.08, .5, 'Fixed value\n$u=0$', ha='right', va='center')
    ax.text(1.22, .5, 'Normal derivative\n$\\partial_n u=1$', ha='left', va='center')
    ax.text(.5, 1.07, 'Zero normal derivative', ha='center')
    ax.text(.5, -.12, 'Zero normal derivative', ha='center')
    ax.set_aspect('equal')
    ax.set(xlim=(-.45, 1.85), ylim=(-.2, 1.2))
    ax.set_axis_off()
    fig.colorbar(image, ax=ax, fraction=.03, pad=.01, shrink=.7, label='$u$')
    fig.tight_layout(pad=1)
    fig.savefig(output / 'MixedBoundary.png')
    plt.close(fig)


def render_poisson(root, output):
    xy, triangles, fields, _ = snapshot(root, 'Poisson')
    u = fields['u'].reshape(-1)
    fig = plt.figure(figsize=(6, 4.8), dpi=220)
    ax = fig.add_subplot(projection='3d')
    surface = ax.plot_trisurf(*xy.T, u, triangles=triangles, shade=False,
                             linewidth=.28, edgecolor=(7/255, 31/255, 70/255, .55))
    surface.set_facecolors(SEQUENTIAL(plt.Normalize(0, u.max())(u[triangles].mean(axis=1))))
    ax.set_box_aspect((1, 1, .7))
    ax.view_init(elev=27, azim=-55)
    ax.set_axis_off()
    ax.set_position([-.08, -.08, 1.16, 1.16])
    fig.savefig(output / 'Poisson.png', bbox_inches='tight', pad_inches=.02)
    plt.close(fig)


def render_elasticity(root, output):
    xy, triangles, fields, _ = snapshot(root, 'Elasticity')
    displacement = fields['displacement']
    magnitude = np.linalg.norm(displacement, axis=1)
    warped = xy + 10 * displacement
    fig, ax = plt.subplots(figsize=(9, 3.9), dpi=200)
    field = ax.tripcolor(*warped.T, triangles, magnitude, shading='gouraud', cmap=SEQUENTIAL)
    ax.triplot(*warped.T, triangles, color='#17679A', linewidth=.26, alpha=.65)
    ax.plot([0, 4, 4, 0, 0], [0, 0, 1, 1, 0], color='#718096',
            linestyle=(0, (4, 3)), linewidth=1, label='Undeformed beam')
    ax.plot([0, 0], [0, 1], color='#718096', linewidth=3)
    ax.annotate('Fixed', xy=(0, .5), xytext=(-.25, .5), ha='right', va='center')
    ax.set(aspect='equal', xlim=(-.65, 4.15), ylim=(-.55, 1.2), xlabel='$x$', ylabel='$y$')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(loc='lower left', labelcolor='#718096', fontsize=9)
    fig.colorbar(field, ax=ax, pad=.025, fraction=.035, shrink=.8,
                 label='Displacement magnitude')
    fig.tight_layout(pad=1.3)
    fig.savefig(output / 'Elasticity.png')
    plt.close(fig)


def render(root, output):
    output.mkdir(parents=True, exist_ok=True)
    render_poisson(root, output)
    render_elasticity(root, output)
    if (root / 'MixedBoundary/MixedBoundary.xdmf').exists():
        render_mixed_boundary(root, output)
    xy, triangles, fields, _ = snapshot(root, 'Helmholtz')
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), dpi=180)
    for ax, key, title in zip(axes, ['real', 'imaginary'], ['Real part', 'Imaginary part']):
        image = planar(ax, xy, triangles, fields[key], title)
        fig.colorbar(image, ax=ax, fraction=.045, pad=.03)
    fig.tight_layout()
    fig.savefig(output / 'Helmholtz.png')
    plt.close(fig)

    xy, triangles, fields, _ = snapshot(root, 'Periodic')
    fig, ax = plt.subplots(figsize=(5.8, 5.2), dpi=180)
    image = planar(ax, xy, triangles, fields['u'], 'Periodic Poisson solution')
    ticks = [0, np.pi, 2 * np.pi]
    ax.set_xticks(ticks, ['0', '$\\pi$', '$2\\pi$'])
    ax.set_yticks(ticks, ['0', '$\\pi$', '$2\\pi$'])
    fig.colorbar(image, ax=ax, fraction=.045, pad=.04)
    fig.tight_layout()
    fig.savefig(output / 'Periodic.png')
    plt.close(fig)

    xyz, triangles, fields, _ = snapshot(root, 'SurfaceMesh')
    u = fields['u']
    fig = plt.figure(figsize=(6, 5.2), dpi=180)
    ax = fig.add_subplot(projection='3d')
    surface = ax.plot_trisurf(*xyz.T, triangles=triangles, shade=False,
                             linewidth=.2, edgecolor=(.1, .3, .45, .35))
    surface.set_facecolors(DIVERGING((u[triangles].mean(axis=1) + .5)))
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=24, azim=38)
    ax.set_axis_off()
    fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-.5, .5), cmap=DIVERGING),
                 ax=ax, shrink=.65, pad=0, label='$u$')
    fig.tight_layout()
    fig.savefig(output / 'SurfaceMesh.png')
    plt.close(fig)

    for name in ['Stokes', 'NavierStokes']:
        fig = flow_figure(root, name)
        fig.savefig(output / (name + '.png'), bbox_inches='tight', pad_inches=.1)
        plt.close(fig)

    grids = ET.parse(root / 'NavierStokes/NavierStokes.xdmf').findall('.//Grid[@GridType="Uniform"]')
    # Every second dt=0.05 snapshot at 10 fps preserves simulation-time playback.
    # VP9/WebM retains alpha; the MP4/H.264 format does not support it.
    writer = FFMpegWriter(fps=10, codec='libvpx-vp9', bitrate=-1,
                          extra_args=['-pix_fmt', 'yuva420p', '-crf', '30',
                                      '-b:v', '0', '-auto-alt-ref', '0'])
    first = flow_figure(root, 'NavierStokes', 0)
    with writer.saving(first, str(output / 'NavierStokes.webm'), dpi=140):
        for index in range(0, len(grids), 2):
            fig = flow_figure(root, 'NavierStokes', index)
            writer.fig = fig
            writer.grab_frame(transparent=True, facecolor='none')
            plt.close(fig)
    plt.close(first)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--example', choices=['all', 'mixed-boundary'], default='all')
    args = parser.parse_args()
    if args.example == 'mixed-boundary':
        render_mixed_boundary(args.root, args.output)
    else:
        render(args.root, args.output)
