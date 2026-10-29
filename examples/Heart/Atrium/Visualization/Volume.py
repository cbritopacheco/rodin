#!/usr/bin/env pvpython
# Copyright Carlos BRITO PACHECO 2026.
# Distributed under the Boost Software License, Version 1.0.
"""Render solved atrial velocity states as a volume animation using ParaView."""
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import paraview.simple as pv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='Fluid XDMF time series.')
    parser.add_argument('output', type=Path, help='Output MP4 movie or PVSM scene.')
    parser.add_argument('--solid', type=Path, help='Solid XDMF; otherwise use sibling Solid.xdmf when present.')
    parser.add_argument('--state-only', action='store_true',
                        help='Save a ParaView scene without rendering or encoding a movie.')
    parser.add_argument('--units', choices=['m/s', 'cm/s'], default='cm/s')
    parser.add_argument('--fps', type=float, default=12)
    parser.add_argument('--resolution', type=int, default=128,
                        help='Samples along each axis of the visualization grid.')
    parser.add_argument('--size', type=int, nargs=2, default=[1920, 900])
    parser.add_argument('--max-frames', type=int, default=0,
                        help='Subsample solved states; zero renders every state.')
    parser.add_argument('--vmax', type=float, help='Fixed legend maximum in selected units.')
    parser.add_argument('--title', default='Blood velocity (volume rendering)')
    args = parser.parse_args()
    if args.fps <= 0 or args.resolution < 8 or args.max_frames < 0:
        parser.error('fps must be positive, resolution at least 8, and max-frames nonnegative')
    if any(n <= 0 or n % 2 for n in args.size):
        parser.error('image width and height must be positive even integers for MP4')
    if args.vmax is not None and (not math.isfinite(args.vmax) or args.vmax <= 0):
        parser.error('vmax must be finite and positive')
    encoder = shutil.which('ffmpeg')
    if not encoder and not args.state_only:
        parser.error('ffmpeg is required to encode MP4')
    source, output = args.input.resolve(), args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    pv._DisableFirstRenderCameraReset()
    reader = pv.XDMFReader(FileNames=[str(source)])
    reader.UpdatePipelineInformation()
    solid_path = args.solid.resolve() if args.solid else source.with_name('Solid.xdmf')
    if args.solid and not solid_path.exists():
        parser.error(f'Solid series does not exist: {solid_path}')
    solid_reader = None
    solid_surface = None
    if solid_path.exists():
        solid_reader = pv.XDMFReader(registrationName='Solid reference mesh (mm)', FileNames=[str(solid_path)])
        solid_reader.UpdatePipelineInformation()
    merged = pv.MergeBlocks(Input=reader)
    merged.MergePoints = 1
    times = list(reader.TimestepValues) or [0.0]
    if solid_reader is not None:
        solid_times = list(solid_reader.TimestepValues) or [0.0]
        times = [t for t in times if any(math.isclose(t, u, rel_tol=1e-10, abs_tol=1e-12) for u in solid_times)]
        if not times:
            parser.error('Fluid and solid series have no common saved times')
        solid_reader.UpdatePipeline(times[0])
        if 'Displacement' not in solid_reader.PointData.keys():
            parser.error('Solid series has no nodal Displacement field')
        solid_merged = pv.MergeBlocks(Input=solid_reader)
        solid_merged.MergePoints = 1
        warped = pv.WarpByVector(registrationName='Solid physical displacement (scale 1)', Input=solid_merged)
        warped.Vectors = ['POINTS', 'Displacement']
        warped.ScaleFactor = 1.0
        displacement = pv.Calculator(Input=warped)
        displacement.ResultArrayName = 'DisplacementMagnitude'
        displacement.Function = 'mag(Displacement)'
        solid_metres = pv.Transform(registrationName='Deformed solid coordinates (m)', Input=displacement)
        solid_metres.Transform.Scale = [.001, .001, .001]
        # Scaling geometry must not rescale the displacement array, whose legend is mm.
        solid_metres.TransformAllInputVectors = 0
        solid_surface = pv.ExtractSurface(Input=solid_metres)
    if args.max_frames and len(times) > args.max_frames:
        indices = [round(i * (len(times)-1) / max(args.max_frames-1, 1))
                   for i in range(args.max_frames)]
        times = [times[i] for i in indices]
    reader.UpdatePipeline(times[0])
    velocity = 'Velocity' if 'Velocity' in reader.PointData.keys() else 'velocity'
    speed = pv.Calculator(Input=merged)
    speed.ResultArrayName = 'Speed'
    factor = 100 if args.units == 'cm/s' else 1
    speed.Function = f'{factor}*mag({velocity})'
    maximum = 0.0
    displacement_maximum = 0.0
    bounds = None
    for t in times:
        speed.UpdatePipeline(t)
        maximum = max(maximum, speed.PointData['Speed'].GetRange()[1])
        pipelines = [speed]
        if solid_surface is not None:
            solid_surface.UpdatePipeline(t)
            displacement_maximum = max(displacement_maximum, solid_surface.PointData['DisplacementMagnitude'].GetRange()[1])
            pipelines.append(solid_surface)
        for pipeline in pipelines:
            current = pipeline.GetDataInformation().GetBounds()
            bounds = list(current) if bounds is None else [min(bounds[i], current[i]) if i % 2 == 0 else max(bounds[i], current[i]) for i in range(6)]
    vmax = args.vmax if args.vmax is not None else maximum
    if not math.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    volume = pv.ResampleToImage(Input=speed)
    volume.SamplingDimensions = [args.resolution] * 3
    volume.UseInputBounds = 1
    surface = pv.ExtractSurface(Input=merged)
    lut = pv.GetColorTransferFunction('Speed')
    speed.UpdatePipeline(times[0])
    if solid_surface is not None:
        displacement_lut = pv.GetColorTransferFunction('DisplacementMagnitude')
        solid_surface.UpdatePipeline(times[0])
    diagonal = math.sqrt(sum((bounds[2*i+1]-bounds[2*i])**2 for i in range(3)))
    center = [(bounds[2*i]+bounds[2*i+1])/2 for i in range(3)]
    views = [pv.CreateView('RenderView') for _ in range(2)]
    layout = pv.CreateLayout(name='Atrium Flow')
    layout.SplitHorizontal(0, .5)
    for location, view in enumerate(views, start=1):
        layout.AssignView(location, view)
    layout.SetSize(*args.size)
    time = pv.AnnotateTimeFilter(Input=reader)
    time.Format = 't = {time:.3f} s'
    for view, label, direction in zip(views, ['Back', 'Front'], [1, -1]):
        view.UseColorPaletteForBackground = 0
        view.Background = [.015, .025, .055]
        view.OrientationAxesVisibility = 0
        display = pv.Show(volume, view)
        display.Representation = 'Volume'
        pv.ColorBy(display, ('POINTS', 'Speed'))
        display.ScalarOpacityUnitDistance = diagonal / 60
        skin = pv.Show(surface, view)
        skin.Representation = 'Surface'
        pv.ColorBy(skin, None)
        skin.DiffuseColor = [.85, .92, 1.]
        skin.Opacity = .06 if solid_surface is not None else .16
        skin.Ambient = .3
        skin.Diffuse = .7
        skin.Specular = .3
        skin.SpecularPower = 25
        if solid_surface is not None:
            wall = pv.Show(solid_surface, view)
            wall.Representation = 'Surface'
            wall.Opacity = .45
            wall.Specular = .3
            pv.ColorBy(wall, ('POINTS', 'DisplacementMagnitude'))
            wall.SetScalarBarVisibility(view, True)
            displacement_bar = pv.GetScalarBar(displacement_lut, view)
            displacement_bar.Title = 'Solid displacement'
            displacement_bar.ComponentTitle = '(mm)'
            displacement_bar.WindowLocation = 'Lower Left Corner'
            displacement_bar.TitleColor = [.95, .95, .98]
            displacement_bar.LabelColor = [.95, .95, .98]
        display.SetScalarBarVisibility(view, True)
        bar = pv.GetScalarBar(lut, view)
        bar.Title = 'Speed'
        bar.ComponentTitle = f'({args.units})'
        bar.TitleColor = [.95, .95, .98]
        bar.LabelColor = [.95, .95, .98]
        bar.WindowLocation = 'Lower Right Corner'
        title = pv.Text(Text=f'{args.title}\n{label} view')
        title_display = pv.Show(title, view)
        title_display.WindowLocation = 'Upper Left Corner'
        title_display.FontSize = 16
        title_display.Color = [.95, .95, .98]
        time_display = pv.Show(time, view)
        time_display.WindowLocation = 'Upper Right Corner'
        time_display.FontSize = 14
        time_display.Color = [.95, .95, .98]
        view.CameraFocalPoint = center
        view.CenterOfRotation = center
        view.CameraPosition = [center[0]+direction*diagonal,
                               center[1]+direction*diagonal, center[2]+diagonal]
        view.CameraViewUp = [0, 0, 1]
        view.CameraParallelProjection = 1
        view.ViewTime = times[0]
        pv.ResetCamera(view)
    if solid_surface is not None:
        dmax = max(displacement_maximum, 1e-12)
        displacement_lut.ColorSpace = 'RGB'
        displacement_lut.RGBPoints = [0, .267, .005, .329,
            .25*dmax, .230, .322, .546, .5*dmax, .128, .567, .551,
            .75*dmax, .369, .789, .383, dmax, .993, .906, .144]
        displacement_lut.RescaleTransferFunction(0, max(displacement_maximum, 1e-12))
        displacement_lut.AutomaticRescaleRangeMode = 'Never'
    lut.ColorSpace = 'RGB'
    lut.RGBPoints = [0, .10, .35, 1., .25*vmax, .05, 1., 1.,
                     .5*vmax, 1., 1., .55, .75*vmax, 1., .70, .12,
                     vmax, 1., .28, .12]
    lut.RescaleTransferFunction(0, vmax)
    opacity = pv.GetOpacityTransferFunction('Speed')
    opacity.Points = [0, 0, .5, 0, .02*vmax, .002, .5, 0,
                      .25*vmax, .04, .5, 0, .5*vmax, .16, .5, 0,
                      vmax, .6, .5, 0]
    opacity.RescaleTransferFunction(0, vmax)
    lut.AutomaticRescaleRangeMode = 'Never'
    pv.SetActiveView(views[0])
    animation = pv.GetAnimationScene()
    animation.UpdateAnimationUsingDataTimeSteps()
    animation.AnimationTime = times[0]
    state = output.with_suffix('.pvsm')
    pv.SaveState(str(state))
    if args.state_only:
        print('Scene:', state)
        return
    with tempfile.TemporaryDirectory(prefix='atrium-volume-') as tmp:
        for i, t in enumerate(times):
            animation.AnimationTime = t
            volume.UpdatePipeline(t)
            surface.UpdatePipeline(t)
            if solid_surface is not None:
                solid_surface.UpdatePipeline(t)
            pv.RenderAllViews()
            frame = Path(tmp)/f'frame-{i:05d}.png'
            pv.SaveScreenshot(str(frame), layout, ImageResolution=args.size)
            if i == len(times)-1:
                shutil.copyfile(frame, output.with_suffix('.png'))
            print(f'Rendered {i+1}/{len(times)}: t={t:.6g} s', flush=True)
        subprocess.run([encoder, '-y', '-v', 'error', '-framerate', str(args.fps),
                        '-i', str(Path(tmp)/'frame-%05d.png'), '-c:v', 'libx264',
                        '-pix_fmt', 'yuv420p', '-crf', '18', str(output)], check=True)
    metadata = dict(solid_source=str(solid_path) if solid_reader is not None else None, displacement_scale=1.0, displacement_units='mm', source=str(source), units=args.units, velocity_factor=factor,
                    legend_maximum=vmax, selected_times_s=times, fps=args.fps,
                    movie_duration_s=len(times)/args.fps,
                    sampling_dimensions=[args.resolution]*3,
                    views=['Back', 'Front'])
    output.with_suffix('.json').write_text(json.dumps(metadata, indent=2))
    print('Movie:', output)


if __name__ == '__main__':
    main()
