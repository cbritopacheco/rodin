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
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.tri import LinearTriInterpolator, Triangulation
import numpy as np

SEQUENTIAL = LinearSegmentedColormap.from_list(
    'Rodin', ['#071F46', '#17679A', '#13B6D0', '#F4F8FB'])
DIVERGING = LinearSegmentedColormap.from_list(
    'Rodin signed', ['#071F46', '#F4F8FB', '#13B6D0'])
plt.rcParams.update({'font.size': 10, 'text.color': '#27323f',
                     'axes.labelcolor': '#27323f', 'xtick.color': '#627083',
                     'ytick.color': '#627083', 'axes.edgecolor': '#d1d8e0',
                     'figure.facecolor': '#F4F8FB', 'axes.facecolor': '#F4F8FB'})


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
    image = axes[0].tripcolor(tri, speed, shading='gouraud', cmap=SEQUENTIAL, vmin=0, vmax=1)
    axes[0].set_title('Speed', pad=22)
    axes[1].streamplot(coordinates, coordinates, ux, uy, density=1.8,
                       color='#17679A', linewidth=.8, arrowsize=.8)
    axes[1].set_title('Instantaneous streamlines', pad=22)
    for ax in axes:
        ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='$x$', ylabel='$y$', aspect='equal')
        ax.annotate('', xy=(.82, 1.035), xytext=(.18, 1.035), xycoords='axes fraction',
                    arrowprops={'arrowstyle': '->', 'color': '#17679A', 'lw': 1.4})
    fig.colorbar(image, ax=axes[0], fraction=.045, pad=.03)
    if name == 'NavierStokes':
        fig.suptitle(f'Lid-driven cavity · Re = 100 · t = {time:.2f}', y=.99)
    fig.tight_layout(pad=1.2, rect=(0, 0, 1, .93 if name == 'NavierStokes' else 1))
    return fig


def render(root, output):
    output.mkdir(parents=True, exist_ok=True)
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
    writer = FFMpegWriter(fps=10, codec='libx264', bitrate=1800,
                          extra_args=['-pix_fmt', 'yuv420p', '-movflags', '+faststart'])
    first = flow_figure(root, 'NavierStokes', 0)
    with writer.saving(first, str(output / 'NavierStokes.mp4'), dpi=140):
        for index in range(0, len(grids), 2):
            fig = flow_figure(root, 'NavierStokes', index)
            writer.fig = fig
            writer.grab_frame(facecolor=fig.get_facecolor())
            plt.close(fig)
    plt.close(first)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    render(args.root, args.output)
