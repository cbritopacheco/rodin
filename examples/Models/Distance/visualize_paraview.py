#!/usr/bin/env pvpython
"""Render the three-circle distance normalization and its reference with ParaView."""
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
args.output.mkdir(parents=True, exist_ok=True)
source = XDMFReader(FileNames=[str((args.root / 'DistanceNormalization.xdmf').resolve())])
source.UpdatePipeline()
error = Calculator(Input=source)
error.ResultArrayName = 'error'
error.Function = 'abs(normalized-reference)'
error.UpdatePipeline()
for name, title, label in [('reference', 'Reference distance', r'$d$'),
                           ('normalized', 'Normalized potential', r'$d_h$'),
                           ('error', 'Absolute nodal error', r'$\mathrm{abs}(d_h-d)$')]:
    target = style.view()
    target.CameraParallelScale = .74
    target.CameraFocalPoint = [.55, .5, 0]
    target.CameraPosition = [.55, .5, 5]
    is_error = name == 'error'
    style.color(error if is_error else source, target, name,
                (0, .12) if is_error else (-.55, .55), signed=not is_error, title=label)
    contour = Contour(Input=source)
    contour.ContourBy = ['POINTS', 'potential']
    contour.Isosurfaces = [0]
    line = Show(contour, target)
    ColorBy(line, None)
    line.AmbientColor = [.027, .122, .275]
    line.DiffuseColor = [.027, .122, .275]
    line.LineWidth = 1.6
    style.text(target, title)
    style.save(args.output / ('DistanceNormalization-' + name + '.png'), target, (850, 740))
    Delete(target)
