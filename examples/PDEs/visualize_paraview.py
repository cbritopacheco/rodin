#!/usr/bin/env pvpython
"""Render the PDE tutorials from XDMF output with ParaView.

Run with pvpython, not ordinary Python. PNG frames retain alpha; FFmpeg encodes
those frames into the cavity's alpha-capable VP9/WebM animation.
"""
import argparse
from pathlib import Path
import subprocess
import os
import math
import matplotlib
from PIL import Image

matplotlib.rcParams['mathtext.fontset'] = 'cm'

from paraview.simple import *
from paraview.simple import _DisableFirstRenderCameraReset

FONT = str(Path(__file__).resolve().parents[2] / 'brand/fonts/static/Manrope-Regular.ttf')
TEXT = [value / 255 for value in ([222, 227, 235] if os.environ.get('RODIN_PLOT_THEME') == 'dark' else [48, 59, 73])]
EDGE = [.44, .50, .56]
BLUE = [.09, .40, .60]
CYAN = [.075, .714, .816]
_DisableFirstRenderCameraReset()


def reader(root, name):
    source = XDMFReader(FileNames=[str((root / name / (name + '.xdmf')).resolve())])
    source.UpdatePipeline()
    return source


def view():
    result = CreateView('RenderView')
    result.ViewSize = [640, 640]
    result.OrientationAxesVisibility = 0
    result.UseColorPaletteForBackground = 0
    result.Background = [1, 1, 1]
    result.CameraParallelProjection = 1
    result.CameraPosition = [.5, .54, 5]
    result.CameraFocalPoint = [.5, .54, 0]
    result.CameraParallelScale = .76
    return result


def text(target, value, position=(.5, .94), size=22):
    source = Text(Text=value)
    display = Show(source, target)
    display.WindowLocation = 'Lower Center' if position[1] < .1 else 'Upper Center'
    display.FontFamily = 'File'
    display.FontFile = FONT
    display.FontSize = size
    display.Color = TEXT
    display.Justification = 'Center'
    return source


def color(source, target, array, limits, signed=False, title=None, association=None):
    """Plot the exported association; cell fields remain constant on each cell."""
    if association is None:
        point_field = array in source.PointData.keys()
        cell_field = array in source.CellData.keys()
        if point_field == cell_field:
            raise ValueError(f"Select POINTS or CELLS explicitly for {array}")
        association = 'POINTS' if point_field else 'CELLS'
    display = Show(source, target)
    display.Representation = 'Surface'
    display.Ambient = 1
    display.Diffuse = 0
    ColorBy(display, (association, array))
    if association == 'CELLS':
        display.InterpolateScalarsBeforeMapping = 0
    lut = GetColorTransferFunction(array)
    # Use the same continuous palette across the tutorial fields.
    lut.ApplyPreset('Turbo', True)
    lut.RescaleTransferFunction(*limits)
    display.SetScalarBarVisibility(target, True)
    bar = GetScalarBar(lut, target)
    bar.Title = title or rf'${array}$'
    bar.ComponentTitle = ''
    bar.TitleFontFamily = 'File'
    bar.TitleFontFile = FONT
    bar.LabelFontFamily = 'File'
    bar.LabelFontFile = FONT
    bar.TitleColor = TEXT
    bar.LabelColor = TEXT
    bar.TitleFontSize = 20
    bar.LabelFontSize = 18
    bar.WindowLocation = 'Any Location'
    bar.Position = [.84, .19]
    bar.AutomaticLabelFormat = 0
    bar.LabelFormat = '{:.2g}'
    bar.RangeLabelFormat = '{:.2g}'
    bar.AllowOverlappingLabels = 0
    bar.UseCustomLabels = 1
    bar.CustomLabels = [limits[0] + i * (limits[1] - limits[0]) / 4 for i in range(5)]
    bar.AddRangeLabels = 0
    bar.ScalarBarThickness = 12
    bar.ScalarBarLength = .58
    return display


def save(path, target, size=(1280, 640), crop=True):
    if target.GetXMLName() == 'RenderView':
        target.ViewSize = list(size)
    SaveScreenshot(str(path), target, ImageResolution=list(size), TransparentBackground=1)
    if crop and path.suffix.lower() == '.png' and target.GetXMLName() == 'RenderView':
        # Trim unused transparent margins from still figures.
        with Image.open(path) as image:
            bounds = image.getchannel('A').getbbox()
            if bounds:
                image.crop((max(0, bounds[0] - 12), max(0, bounds[1] - 12),
                            min(image.width, bounds[2] + 12),
                            min(image.height, bounds[3] + 12))).save(path)


def save_panels(path, views, names, size=(720, 640)):
    """Export each scientific panel so HTML controls rows and responsive layout."""
    for target, name in zip(views, names):
        save(path.with_name(path.stem + '-' + name + path.suffix), target, size)


def flow(root, name, output, frames=None, fps=30):
    source = reader(root, name)
    times = list(source.TimestepValues)
    if frames:
        source = TemporalInterpolator(Input=source)
        source.UpdatePipeline(times[0])
    left, right = view(), view()
    layout = CreateLayout(name=name)
    AssignViewToLayout(left, layout, hint=0)
    layout.SplitHorizontal(0, .5)
    AssignViewToLayout(right, layout, hint=2)
    layout.SetSize(1280, 680)
    peak_speed = 2 if name == 'NavierStokes' else 1
    color(source, left, 'velocity', (0, peak_speed), title='Speed')
    speed_bar = GetScalarBar(GetColorTransferFunction('velocity'), left)
    speed_bar.Position = [.90, .19]
    speed_bar.Title = ''
    outline = Outline(Input=source)
    for target in [left, right]:
        display = Show(outline, target)
        display.DiffuseColor = TEXT
        display.LineWidth = 1
        text(target, 'Speed' if target == left else 'Instantaneous streamlines', (.5, .96))
    seeds = Line(Point1=[.5, .06, 0], Point2=[.5, .94, 0], Resolution=14)
    streams = StreamTracerWithCustomSource(Input=source, SeedSource=seeds)
    streams.Vectors = ['POINTS', 'velocity']
    streams.IntegrationDirection = 'BOTH'
    streams.SurfaceStreamlines = 1
    streams.IntegrationStepUnit = 'Length'
    streams.InitialStepLength = .005
    streams.MaximumStreamlineLength = 1.5
    streams.MaximumSteps = 4000
    streams.TerminalSpeed = 1e-8
    display = Show(streams, right)
    ColorBy(display, None)
    display.DiffuseColor = BLUE
    display.AmbientColor = BLUE
    display.LineWidth = 1.5
    arrows = Glyph(Input=streams, GlyphType='Arrow')
    arrows.OrientationArray = ['POINTS', 'velocity']
    arrows.ScaleArray = ['POINTS', 'No scale array']
    arrows.ScaleFactor = .018
    arrows.GlyphMode = 'Every Nth Point'
    arrows.Stride = 250
    display = Show(arrows, right)
    ColorBy(display, None)
    display.DiffuseColor = BLUE
    display.AmbientColor = BLUE
    lid = PlotOverLine(Input=source, Point1=[0, 1, 0], Point2=[1, 1, 0], Resolution=128)
    profile = Calculator(Input=lid)
    profile.Function = f'coordsX*iHat+(1.04+{0.08 / peak_speed}*mag(velocity))*jHat'
    profile.CoordinateResults = 1
    for target in [left, right]:
        display = Show(profile, target)
        ColorBy(display, None)
        display.DiffuseColor = CYAN
        display.AmbientColor = CYAN
        display.LineWidth = 3
    label = text(left, '', (.5, .02), 18)
    selected = ([(index, min(times[0] + index / fps, times[-1]))
                 for index in range(math.ceil((times[-1] - times[0]) * fps) + 1)]
                if frames else [(len(times)-1, times[-1])])
    if frames:
        frames.mkdir(parents=True, exist_ok=True)
    for frame, (_, time) in enumerate(selected):
        source.UpdatePipeline(time)
        left.ViewTime = right.ViewTime = time
        label.Text = rf'$\mathrm{{Re}}=200,\quad t={time:.2f}$' if name == 'NavierStokes' else 'Steady Stokes flow'
        streams.UpdatePipeline(time)
        if frames:
            save(frames / f'{frame:04d}.png', layout, (1280, 680))
        if time == times[-1]:
            save(output / (name + '.png'), layout, (1280, 680))
            save_panels(output / (name + '.png'), [left, right], ['speed', 'streamlines'])
        print(name, time, flush=True)
    Delete(layout)
    Delete(left)
    Delete(right)


def scalar(root, name, output):
    source = reader(root, name)
    source.UpdatePipeline()
    target = view()
    dataset = source
    if name == 'Elasticity':
        dataset = WarpByVector(Input=source, Vectors=['POINTS', 'displacement'], ScaleFactor=10)
        dataset.UpdatePipeline()
        magnitude = source.PointData['displacement'].GetRange(-1)[1]
        display = color(dataset, target, 'displacement', (0, magnitude), title='Displacement')
        display.Representation = 'Surface With Edges'
        display.EdgeColor = [.22, .34, .43]
        display.LineWidth = .4
        original = Show(Outline(Input=source), target)
        original.DiffuseColor = TEXT
        target.CameraPosition = [2, .4, 8]
        target.CameraFocalPoint = [2, .4, 0]
        target.CameraParallelScale = 1.9
        bar = GetScalarBar(GetColorTransferFunction('displacement'), target)
        bar.Orientation = 'Horizontal'
        bar.Position = [.22, .08]
        bar.ScalarBarLength = .55
        text(target, 'Displacement amplified 10 times', (.5, .94))
    elif name == 'Poisson':
        dataset = WarpByScalar(Input=source, Scalars=['POINTS', 'u'], ScaleFactor=8)
        dataset.UpdatePipeline()
        display = color(dataset, target, 'u', source.PointData['u'].GetRange(), title=r'$u$')
        display.Representation = 'Surface With Edges'
        display.EdgeColor = [.18, .30, .40]
        display.LineWidth = .4
        target.CameraPosition = [2.2, -2.8, 2.3]
        target.CameraFocalPoint = [.5, .5, .2]
        target.CameraViewUp = [0, 0, 1]
        target.CameraParallelScale = .85
        bar = GetScalarBar(GetColorTransferFunction('u'), target)
        bar.Position = [.87, .19]
    elif name == 'SurfaceMesh':
        display = color(source, target, 'u', (-.5, .5), signed=True)
        display.Representation = 'Surface With Edges'
        display.EdgeColor = [.18, .30, .40]
        display.LineWidth = .35
        target.CameraPosition = [3, -4, 2.5]
        target.CameraFocalPoint = [0, 0, 0]
        target.CameraViewUp = [0, 0, 1]
        target.CameraParallelScale = 1.5
    elif name == 'Periodic':
        color(source, target, 'u', (-1, 1), signed=True)
        target.CameraPosition = [3.14, 3.14, 8]
        target.CameraFocalPoint = [3.14, 3.14, 0]
        target.CameraParallelScale = 4.6
    elif name == 'MixedBoundary':
        color(source, target, 'u', (0, 1))
        text(target, 'Mixed boundary conditions', (.5, .94))
        text(target, r'Left: $u=0$     Right: $\partial u/\partial n=1$' + '\n' + r'Top/bottom: $\partial u/\partial n=0$',
             (.5, .02), 18)
        for x in [0, 1]:
            edge = Show(Line(Point1=[x, 0, 0], Point2=[x, 1, 0]), target)
            edge.DiffuseColor = BLUE if x == 0 else CYAN
            edge.LineWidth = 3
    else:
        raise ValueError(name)
    size = (800, 640) if name == 'Poisson' else (720, 640)
    save(output / (name + '.png'), target, size)
    Delete(target)


def poisson_readme(root, output, frames=None, fps=30):
    """Render the exported P1 nodal field with interpolated point colors."""
    source = reader(root, 'Poisson')
    target = view()
    surface = WarpByScalar(Input=source, Scalars=['POINTS', 'u'], ScaleFactor=9)
    surface.UpdatePipeline()
    display = color(surface, target, 'u', source.PointData['u'].GetRange(),
                    title=r'$u$', association='POINTS')
    display.InterpolateScalarsBeforeMapping = 1
    display.SetScalarBarVisibility(target, False)
    display.Representation = 'Surface With Edges'
    display.EdgeColor = [0, 0, 0]
    display.LineWidth = .4
    target.CameraPosition = [2.2, -2.8, 2.3]
    target.CameraFocalPoint = [.5, .5, .2]
    target.CameraViewUp = [0, 0, 1]
    target.CameraParallelScale = .85
    save(output / 'Poisson.png', target, (1000, 800))
    if frames is not None:
        frames.mkdir(parents=True, exist_ok=True)
        count = 6 * fps
        center = list(target.CameraFocalPoint)
        offset = [target.CameraPosition[i] - center[i] for i in range(3)]
        bounds = None
        for index in range(count):
            angle = 2 * math.pi * index / count
            target.CameraPosition = [
                center[0] + offset[0] * math.cos(angle) - offset[1] * math.sin(angle),
                center[1] + offset[0] * math.sin(angle) + offset[1] * math.cos(angle),
                center[2] + offset[2]]
            path = frames / f'{index:04d}.png'
            save(path, target, (720, 576), crop=False)
            with Image.open(path) as image:
                box = image.getchannel('A').getbbox()
            bounds = box if bounds is None else (
                min(bounds[0], box[0]), min(bounds[1], box[1]),
                max(bounds[2], box[2]), max(bounds[3], box[3]))
            if index % fps == 0:
                print('README rotation:', index, '/', count, flush=True)
        # One crop for the whole orbit prevents the image from shifting or resizing.
        bounds = (max(0, bounds[0] - 12), max(0, bounds[1] - 12),
                  min(720, bounds[2] + 12), min(576, bounds[3] + 12))
        for index in range(count):
            path = frames / f'{index:04d}.png'
            with Image.open(path) as image:
                image.crop(bounds).save(path)
        pattern = str(frames / '%04d.png')
        palette = str(frames / 'palette.png')
        subprocess.run(['ffmpeg', '-y', '-framerate', str(fps), '-t', str(count / fps), '-i', pattern,
                        '-vf', 'palettegen=stats_mode=full:reserve_transparent=1',
                        '-frames:v', '1', '-update', '1', palette], check=True)
        subprocess.run(['ffmpeg', '-y', '-framerate', str(fps), '-t', str(count / fps), '-i', pattern,
                        '-i', palette, '-lavfi',
                        'paletteuse=dither=sierra2_4a:alpha_threshold=128',
                        '-loop', '0', '-final_delay', str(round(100 / fps)),
                        str(output / 'Poisson.gif')], check=True)

    print('README Poisson: point-data P1 interpolation,',
          source.GetDataInformation().GetNumberOfPoints(), 'vertices,',
          source.GetDataInformation().GetNumberOfCells(), 'cells,',
          'nodal range', source.PointData['u'].GetRange(), flush=True)
    Delete(target)


def helmholtz(root, output):
    source = reader(root, 'Helmholtz')
    left, right = view(), view()
    layout = CreateLayout(name='Helmholtz')
    AssignViewToLayout(left, layout, hint=0)
    layout.SplitHorizontal(0, .5)
    AssignViewToLayout(right, layout, hint=2)
    layout.SetSize(1280, 640)
    for target, array in [(left, 'real'), (right, 'imaginary')]:
        color(source, target, array, (-1, 1), signed=True,
              title=r'$\mathrm{Re}(u_h)$' if array == 'real' else r'$\mathrm{Im}(u_h)$')
        text(target, 'Real part' if array == 'real' else 'Imaginary part')
    save(output / 'Helmholtz.png', layout)
    save_panels(output / 'Helmholtz.png', [left, right], ['real', 'imaginary'])
    Delete(layout)
    Delete(left)
    Delete(right)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--frames', type=Path)
    parser.add_argument('--fps', type=int, choices=[30, 60], default=30)
    cases = ['Poisson', 'Elasticity', 'SurfaceMesh', 'Periodic', 'MixedBoundary']
    parser.add_argument('--example', choices=['all', 'cavity', 'Helmholtz', 'Stokes', 'NavierStokes', 'READMEPoisson', 'READMEPoissonRotation'] + cases,
                        default='all')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.example in ['all', 'cavity'] and args.frames is None:
        parser.error('--frames is required when rendering the cavity animation')
    if args.example == 'all':
        for name in cases:
            scalar(args.root, name, args.output)
        helmholtz(args.root, args.output)
        flow(args.root, 'Stokes', args.output)
    elif args.example in cases:
        scalar(args.root, args.example, args.output)
    elif args.example in ['READMEPoisson', 'READMEPoissonRotation']:
        if args.example == 'READMEPoissonRotation' and args.frames is None:
            parser.error('--frames is required for the README rotation')
        poisson_readme(args.root, args.output,
                       args.frames if args.example == 'READMEPoissonRotation' else None, args.fps)
    elif args.example == 'Helmholtz':
        helmholtz(args.root, args.output)
    elif args.example == 'NavierStokes':
        flow(args.root, 'NavierStokes', args.output)
    elif args.example == 'Stokes':
        flow(args.root, 'Stokes', args.output)
    if args.example in ['all', 'cavity']:
        flow(args.root, 'NavierStokes', args.output, args.frames, args.fps)
        subprocess.run(['ffmpeg', '-y', '-framerate', str(args.fps), '-i', str(args.frames / '%04d.png'),
                        '-c:v', 'libvpx-vp9', '-pix_fmt', 'yuva420p', '-crf', '30', '-b:v', '0',
                        '-auto-alt-ref', '0', str(args.output / 'NavierStokes.webm')], check=True)
