#!/usr/bin/env pvpython
"""Render one computed cantilever update from its XDMF history.

Run with ParaView's pvpython. Temporary single-frame XDMF files reference the
original HDF5 arrays; field values and mesh coordinates are retained verbatim.
"""
import argparse
import copy
import importlib.util
from pathlib import Path
import tempfile
import xml.etree.ElementTree as ET

from paraview.simple import *

spec = importlib.util.spec_from_file_location(
    'rodin_figures', Path(__file__).resolve().parents[1] / 'PDEs/visualize_paraview.py')
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)


def snapshot(path, name, index, temporary):
    document = ET.parse(path)
    collection = next(g for g in document.findall('./Domain/Grid') if g.get('Name') == name)
    grids = collection.findall('Grid')
    grid = copy.deepcopy(grids[index])
    for item in grid.findall('.//DataItem'):
        if item.get('Format') == 'HDF':
            filename, dataset = item.text.strip().split(':', 1)
            item.text = str((path.parent / filename).resolve()) + ':' + dataset
    root = ET.Element('Xdmf', Version='3.0')
    ET.SubElement(root, 'Domain').append(grid)
    filename = temporary / f'{name}-{index}.xdmf'
    ET.ElementTree(root).write(filename, encoding='utf-8', xml_declaration=True)
    source = XDMFReader(FileNames=[str(filename)])
    source.UpdatePipeline()
    return source


def mesh(source, target, fill=None):
    display = Show(source, target)
    display.ColorArrayName = ['POINTS', '']
    display.Representation = 'Surface With Edges'
    display.Ambient = 1
    display.Diffuse = 0
    display.AmbientColor = fill or [.91, .96, .97]
    display.EdgeColor = style.EDGE
    display.LineWidth = .5
    return display


def contour(source, field, target, color, width=2):
    curve = Contour(Input=source)
    curve.ContourBy = ['POINTS', field]
    curve.Isosurfaces = [0]
    curve.UpdatePipeline()
    display = Show(curve, target)
    display.ColorArrayName = ['POINTS', '']
    display.AmbientColor = color
    display.DiffuseColor = color
    display.LineWidth = width
    return curve


def render(args):
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='rodin-cantilever-stage-') as folder:
        temporary = Path(folder)
        if args.mode == 'simple':
            current = snapshot(args.history, 'mesh', args.iteration, temporary)
            following = snapshot(args.history, 'mesh', args.iteration + 1, temporary)
            state = current
            vector = 'descent'
            prefix = 'SimpleCantileverStages'
        else:
            current = snapshot(args.history, 'domain', args.iteration, temporary)
            following = snapshot(args.history, 'state', args.iteration + 1, temporary)
            state = snapshot(args.history, 'state', args.iteration, temporary)
            vector = 'dJ'
            prefix = 'LevelSetStages'
        bounds = current.GetDataInformation().GetBounds()
        span = max(bounds[1] - bounds[0], bounds[3] - bounds[2])
        point_count = current.GetDataInformation().GetNumberOfPoints()
        stride = (max(1, (point_count + 99) // 100)
                  if args.glyph_stride is None else args.glyph_stride)
        if stride < 1:
            raise ValueError('Glyph stride must be positive')
        maximum = current.PointData[vector].GetRange(-1)[1]
        if maximum <= 0:
            raise ValueError('A nonzero computed design direction is required')

        def view(title, note, legend=False):
            target = style.view()
            target.CameraFocalPoint = [(bounds[0] + bounds[1]) / 2 + (.1 * span if legend else 0),
                                      (bounds[2] + bounds[3]) / 2, 0]
            target.CameraPosition = [target.CameraFocalPoint[0], target.CameraFocalPoint[1], 5]
            target.CameraParallelScale = .72 * span
            style.text(target, title, (.5, .95), 24)
            style.text(target, note, (.5, .025), 18)
            return target

        def save(target, name):
            style.save(args.output / f'{prefix}-{name}.png', target, (720, 640))
            Delete(target)

        def descent():
            normalized = Calculator(Input=current)
            normalized.ResultArrayName = 'normalized_descent'
            normalized.Function = f'mag({vector})/{maximum:.17g}'
            normalized.UpdatePipeline()
            target = view('Regularized design direction',
                          f'Iteration {args.iteration}; node stride {stride}', legend=True)
            style.color(normalized, target, 'normalized_descent', (0, 1),
                        title=r'$\Vert\theta_h\Vert/\max\Vert\theta_h\Vert$')
            arrows = Glyph(Input=current, GlyphType='Arrow')
            arrows.OrientationArray = ['POINTS', vector]
            arrows.ScaleArray = ['POINTS', 'No scale array']
            arrows.ScaleFactor = .055 * span
            arrows.GlyphMode = 'Every Nth Point'
            arrows.Stride = stride
            display = Show(arrows, target)
            display.ColorArrayName = ['POINTS', '']
            display.AmbientColor = style.CYAN
            display.DiffuseColor = style.CYAN
            display.Opacity = .5
            if args.mode == 'levelset':
                contour(current, 'dist', target, style.BLUE)
            save(target, 'descent')

        if args.mode == 'simple':
            target = view('Elastic state on the current design', f'Iteration {args.iteration}', legend=True)
            limit = state.PointData['displacement'].GetRange(-1)[1]
            style.color(state, target, 'displacement', (0, limit), title=r'$\Vert u_h\Vert$')
            save(target, 'state')
            descent()
            moved = WarpByVector(Input=current, Vectors=['POINTS', vector],
                                 ScaleFactor=args.step / maximum)
            moved.UpdatePipeline()
            target = view('Applied design displacement', rf'$\max\Vert\delta x\Vert={args.step:g}$')
            mesh(moved, target)
            surface = ExtractSurface(Input=current)
            boundary = FeatureEdges(Input=surface)
            boundary.FeatureEdges = 0
            boundary.ManifoldEdges = 0
            boundary.NonManifoldEdges = 0
            boundary.BoundaryEdges = 1
            display = Show(boundary, target)
            display.ColorArrayName = ['POINTS', '']
            display.AmbientColor = style.BLUE
            display.DiffuseColor = style.BLUE
            display.LineWidth = 2
            save(target, 'displaced')
        else:
            target = view('Reconstructed signed distance',
                          rf'Iteration {args.iteration}; colors clipped at $\pm0.05$', legend=True)
            style.color(current, target, 'dist', (-.05, .05), signed=True, title=r'$d_h$')
            contour(current, 'dist', target, style.BLUE)
            save(target, 'distance')
            descent()
            target = view('Transported level set', 'Navy: current interface; cyan: proposed zero contour', legend=True)
            style.color(current, target, 'advect', (-.05, .05), signed=True, title=r'$\phi_h^{\mathrm{new}}$')
            contour(current, 'dist', target, style.BLUE, 2)
            contour(current, 'advect', target, style.CYAN, 3)
            save(target, 'transport')
        target = view('Mesh used by the next elastic state',
                      f'Iteration {args.iteration + 1}; after remeshing and preparation')
        mesh(following, target)
        save(target, 'remeshed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('history', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--mode', choices=['simple', 'levelset'], required=True)
    parser.add_argument('--iteration', type=int, default=0)
    parser.add_argument('--glyph-stride', type=int,
                        help='Sample every Nth node; default targets at most 100 glyphs')
    parser.add_argument('--step', type=float, default=.01,
                        help='SimpleCantilever2D displacement scale; match the run hmin')
    render(parser.parse_args())
