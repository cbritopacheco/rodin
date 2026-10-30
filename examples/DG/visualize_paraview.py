#!/usr/bin/env pvpython
"""Render the cell-local DG output without merging duplicated vertices."""
import argparse
import importlib.util
from pathlib import Path
from paraview.simple import *

spec = importlib.util.spec_from_file_location('rodin_figures',
    Path(__file__).resolve().parents[1] / 'PDEs/visualize_paraview.py')
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('root', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
layout = CreateLayout(name='DG solution')
left, right = style.view(), style.view()
AssignViewToLayout(left, layout, hint=0)
layout.SplitHorizontal(0, .5)
AssignViewToLayout(right, layout, hint=2)
for target, stem, title in [(left,'DGCoarse','4 intervals per axis'),
                            (right,'DGPoisson','32 intervals per axis')]:
    source = XDMFReader(FileNames=[str((args.root / (stem+'.xdmf')).resolve())])
    source.UpdatePipeline()
    target.CameraFocalPoint=[.5,.42,0]
    target.CameraPosition=[.5,.42,5]
    target.CameraParallelScale = .92
    display = style.color(source, target, 'u', (1,4), title=r'$u_h$')
    if stem == 'DGCoarse':
        display.Representation = 'Surface With Edges'
        display.EdgeColor = style.TEXT
    bar = GetScalarBar(GetColorTransferFunction('u'), target)
    bar.Orientation = 'Horizontal'
    bar.Position = [.22,.025]
    bar.ScalarBarLength = .56
    bar.CustomLabels = [1,2,3,4]
    bar.TitleFontSize = 17
    bar.LabelFontSize = 15
    style.text(target,title)
style.save(args.output,layout,(1440,640))

style.save_panels(args.output, [left, right], ['coarse', 'fine'])
