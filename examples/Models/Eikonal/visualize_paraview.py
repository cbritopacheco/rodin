#!/usr/bin/env pvpython
"""Render the checked FMM arrival-time field and isolines with ParaView."""
import argparse
import importlib.util
from pathlib import Path

from paraview.simple import *

style_path = Path(__file__).resolve().parents[2] / 'PDEs/visualize_paraview.py'
spec = importlib.util.spec_from_file_location('rodin_figures', style_path)
style = importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('root', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
source = XDMFReader(FileNames=[str((args.root / 'FMMArrival.xdmf').resolve())])
source.UpdatePipeline()
target = style.view()
target.CameraParallelScale = .72
target.CameraFocalPoint = [.55, .5, 0]
target.CameraPosition = [.55, .5, 5]
style.color(source, target, 'arrival', (0, .75), title=r'$T_h$')
lines = Contour(Input=source)
lines.ContourBy = ['POINTS', 'arrival']
lines.Isosurfaces = [.1, .2, .3, .4, .5, .6]
display = Show(lines, target)
ColorBy(display, None)
display.AmbientColor = [.027, .122, .275]
display.DiffuseColor = [.027, .122, .275]
display.LineWidth = 1.5
style.text(target, 'Arrival time from the center')
args.output.parent.mkdir(parents=True, exist_ok=True)
style.save(args.output, target, (850, 740))
