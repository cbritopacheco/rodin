#!/usr/bin/env pvpython
"""Render the first and last solved cantilever states and its objective history."""
import argparse
import importlib.util
from pathlib import Path
from paraview.simple import *

spec=importlib.util.spec_from_file_location('rodin_figures',
    Path(__file__).resolve().parents[1]/'PDEs/visualize_paraview.py')
style=importlib.util.module_from_spec(spec)
spec.loader.exec_module(style)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('root',type=Path)
parser.add_argument('output',type=Path)
parser.add_argument('--stem',default='SimpleCantilever')
parser.add_argument('--field',default='displacement')
args=parser.parse_args()
source=XDMFReader(FileNames=[str((args.root/(args.stem+'.xdmf')).resolve())])
source.UpdatePipeline()
times=list(source.TimestepValues)
bounds=[]; maximum=0
for time in [times[0],times[-1]]:
    source.UpdatePipeline(time)
    bounds.append(source.GetDataInformation().GetBounds())
    maximum=max(maximum,source.PointData[args.field].GetRange()[1])
x0=min(b[0] for b in bounds);x1=max(b[1] for b in bounds)
y0=min(b[2] for b in bounds);y1=max(b[3] for b in bounds)
layout=CreateLayout(name='Cantilever descent')
left,right=style.view(),style.view()
AssignViewToLayout(left,layout,hint=0)
layout.SplitHorizontal(0,.5)
AssignViewToLayout(right,layout,hint=2)
for target,time in [(left,times[0]),(right,times[-1])]:
    target.ViewTime=time
    target.CameraFocalPoint=[(x0+x1)/2,(y0+y1)/2,0]
    target.CameraPosition=[(x0+x1)/2,(y0+y1)/2,5]
    target.CameraParallelScale=max(x1-x0,y1-y0)*.85
    source.UpdatePipeline(time)
    display=style.color(source,target,args.field,(0,maximum),title=r'$\Vert\mathbf{u}_h\Vert$')
    display.Representation='Surface With Edges';display.EdgeColor=style.TEXT
    bar=GetScalarBar(GetColorTransferFunction(args.field),target)
    bar.Orientation='Horizontal';bar.Position=[.22,.025];bar.ScalarBarLength=.56
    bar.CustomLabels=[0,maximum/2,maximum];bar.TitleFontSize=17;bar.LabelFontSize=15
    style.text(target,f'Iteration {int(time)}')
style.save(args.output,layout,(1440,640))
