#!/usr/bin/env pvpython
"""Render the MFEM round trip and XDMF time series from their actual output.

Usage: pvpython examples/IO/visualize_paraview.py OUTPUT_ROOT IMAGE_DIRECTORY
OUTPUT_ROOT contains MFEM/Square.xdmf and XDMF/TimeEvolution.xdmf.
"""
import argparse
import importlib.util
from pathlib import Path

from paraview.simple import *

style_path = Path(__file__).resolve().parents[1] / 'PDEs/visualize_paraview.py'
spec = importlib.util.spec_from_file_location('rodin_figures', style_path)
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    source = XDMFReader(FileNames=[str((args.root / 'MFEM/Square.xdmf').resolve())])
    source.UpdatePipeline()
    target = style.view()
    target.CameraParallelScale = .9
    display = style.color(source, target, 'u', (0, 2))
    display.Representation = 'Surface With Edges'
    display.EdgeColor = style.TEXT
    bar = GetScalarBar(GetColorTransferFunction('u'), target)
    bar.Position = [.83, .19]
    style.text(target, 'Restored scalar field')
    style.save(args.output / 'MFEM.png', target, (720, 640))
    Delete(target)
    Delete(source)

    source = XDMFReader(FileNames=[str((args.root / 'XDMF/TimeEvolution.xdmf').resolve())])
    layout = CreateLayout(name='Time snapshots')
    left, right = style.view(), style.view()
    AssignViewToLayout(left, layout, hint=0)
    layout.SplitHorizontal(0, .5)
    AssignViewToLayout(right, layout, hint=2)
    for target, time in [(left, 0), (right, .25)]:
        target.CameraFocalPoint = [.5, .5, 0]
        target.CameraPosition = [.5, .5, 5]
        target.CameraParallelScale = .95
        target.ViewTime = time
        source.UpdatePipeline(time)
        style.color(source, target, 'u', (0, 2 ** .5), title='Magnitude')
        bar = GetScalarBar(GetColorTransferFunction('u'), target)
        bar.Orientation = 'Horizontal'
        bar.Position = [.22, .025]
        bar.ScalarBarLength = .56
        bar.TitleFontSize = 17
        bar.LabelFontSize = 15
        bar.CustomLabels = [0, 2 ** .5 / 2, 2 ** .5]
        arrows = Glyph(Input=source, GlyphType='Arrow')
        arrows.OrientationArray = ['POINTS', 'u']
        arrows.ScaleArray = ['POINTS', 'u']
        arrows.ScaleFactor = .065
        arrows.GlyphMode = 'Every Nth Point'
        arrows.Stride = 40
        display = Show(arrows, target)
        ColorBy(display, None)
        display.DiffuseColor = [.06, .13, .21]
        display.AmbientColor = [.06, .13, .21]
        style.text(target, rf'$t = {time:g}$')
    style.save(args.output / 'XDMF.png', layout, (1440, 640))


if __name__ == '__main__':
    main()
