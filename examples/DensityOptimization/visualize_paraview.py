#!/usr/bin/env pvpython
"""Render the last consistent state of the density optimization example."""
import argparse
import importlib.util
from pathlib import Path
import xml.etree.ElementTree as ET
from paraview.simple import *

spec=importlib.util.spec_from_file_location('rodin_figures',
    Path(__file__).resolve().parents[1]/'PDEs/visualize_paraview.py')
style=importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('root',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
# Retain the final metadata entry; HDF5 paths stay relative to the same directory.
tree=ET.parse(args.root/'TemperatureMinimization.xdmf')
collection=tree.find('.//Grid[@CollectionType="Temporal"]')
for grid in list(collection)[:-1]:
    collection.remove(grid)
selected=args.root/'TemperatureMinimization.last.xdmf'
tree.write(selected,encoding='utf-8',xml_declaration=True)
source=XDMFReader(FileNames=[str(selected.resolve())])
source.UpdatePipeline()
layout=CreateLayout(name='Density optimization')
left,right=style.view(),style.view()
AssignViewToLayout(left,layout,hint=0)
layout.SplitHorizontal(0,.5)
AssignViewToLayout(right,layout,hint=2)
for target,field,title,symbol in [(left,'gamma','Material density',r'$\gamma$'),
                                 (right,'u','Temperature',r'$u$')]:
    target.CameraFocalPoint=[.5,.42,0]
    target.CameraPosition=[.5,.42,5]
    target.CameraParallelScale=.95
    maximum=1 if field=='gamma' else source.PointData[field].GetRange()[1]
    style.color(source,target,field,(0,maximum),title=symbol)
    if field=='gamma':
        lut=GetColorTransferFunction(field)
        lut.RGBPoints=[0,.957,.973,.984, .5,.075,.714,.816, 1,.027,.122,.275]
    bar=GetScalarBar(GetColorTransferFunction(field),target)
    bar.Orientation='Horizontal';bar.Position=[.22,.025];bar.ScalarBarLength=.56
    bar.CustomLabels=[0,maximum/2,maximum];bar.TitleFontSize=17;bar.LabelFontSize=15
    style.text(target,title)
style.save(args.output,layout,(1440,640))

style.save_panels(args.output, [left, right], ['density', 'temperature'])
