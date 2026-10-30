#!/usr/bin/env pvpython
"""Render the three actual BackgroundFitting mesh blocks with ParaView."""
import argparse
import importlib.util
from pathlib import Path

from paraview.simple import *

style_path = Path(__file__).resolve().parents[2] / 'PDEs/visualize_paraview.py'
spec = importlib.util.spec_from_file_location('rodin_figures', style_path)
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--reconstruction', help='XDMF stem for a two-block Reconstruction run')
    parser.add_argument('--lobes', type=int, default=0, help='Target lobe frequency')
    args = parser.parse_args()
    stem = args.reconstruction or 'BackgroundFitting'
    source = XDMFReader(FileNames=[str((args.root / (stem + '.xdmf')).resolve())])
    source.UpdatePipeline()
    layout = CreateLayout(name='Background adaptation and interface fitting')
    names = ['background', 'moved'] if args.reconstruction else ['uniform', 'classified', 'fitted']
    titles = ['Reference mesh', 'Quadratic fitted mesh'] if args.reconstruction else [
        'Uniform background', 'Adapted and classified', 'Fitted interface']
    views = [style.view() for _ in names]
    AssignViewToLayout(views[0], layout, hint=0)
    layout.SplitHorizontal(0, 1 / len(names))
    AssignViewToLayout(views[1], layout, hint=2)
    if len(names) == 3:
        layout.SplitHorizontal(2, .5)
        AssignViewToLayout(views[2], layout, hint=6)
    for target, name, title in zip(views, names, titles):
        block = ExtractBlock(Input=source)
        block.Selectors = ['/Root/' + name]
        block.UpdatePipeline()
        display = Show(block, target)
        ColorBy(display, None)
        display.Representation = 'Surface With Edges'
        display.Ambient = 1
        display.Diffuse = 0
        display.AmbientColor = [.91, .94, .96]
        display.EdgeColor = [.44, .50, .56]
        display.LineWidth = .6
        if name != 'uniform':
            inside = Threshold(Input=block)
            inside.Scalars = ['CELLS', 'Attribute']
            inside.LowerThreshold = 1
            inside.UpperThreshold = 1
            inside_display = Show(inside, target)
            ColorBy(inside_display, None)
            inside_display.Representation = 'Surface With Edges'
            inside_display.Ambient = 1
            inside_display.Diffuse = 0
            inside_display.AmbientColor = [.075, .714, .816]
            inside_display.EdgeColor = [.09, .40, .60]
        circle = ProgrammableSource()
        circle.OutputDataSetType = 'vtkPolyData'
        circle.Script = f'''
import math
import vtk
points = vtk.vtkPoints()
line = vtk.vtkPolyLine()
line.GetPointIds().SetNumberOfIds(241)
for i in range(241):
    angle = 2 * math.pi * i / 240
    radius = .25 + (.05 * math.cos({args.lobes} * angle) if {args.lobes} else 0)
    points.InsertNextPoint(.5 + radius * math.cos(angle), .5 + radius * math.sin(angle), .002)
    line.GetPointIds().SetId(i, i)
cells = vtk.vtkCellArray()
cells.InsertNextCell(line)
output = self.GetPolyDataOutput()
output.SetPoints(points)
output.SetLines(cells)
'''
        reference = Show(circle, target)
        reference.AmbientColor = [.027, .122, .275]
        reference.DiffuseColor = [.027, .122, .275]
        reference.LineWidth = 2
        target.CameraFocalPoint = [.5, .5, 0]
        target.CameraPosition = [.5, .5, 5]
        target.CameraParallelScale = .66
        style.text(target, title, size=19)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    style.save(args.output, layout, (600 * len(names), 660))
    # Individual panels allow the documentation to stack stages on phones.
    for target, name in zip(views, names):
        style.save(args.output.with_name(args.output.stem + '-' + name + '.png'),
                   target, (600, 660))


if __name__ == '__main__':
    main()
