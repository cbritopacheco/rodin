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
    parser.add_argument('output', type=Path, help='Output MP4 movie.')
    parser.add_argument('--units', choices=['m/s', 'cm/s'], default='cm/s')
    parser.add_argument('--fps', type=float, default=12)
    parser.add_argument('--resolution', type=int, default=128,
                        help='Samples along each axis of the visualization grid.')
    parser.add_argument('--size', type=int, nargs=2, default=[1280, 720])
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
    if not encoder:
        parser.error('ffmpeg is required to encode MP4')
    source, output = args.input.resolve(), args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    pv._DisableFirstRenderCameraReset()
    reader = pv.XDMFReader(FileNames=[str(source)])
    reader.UpdatePipelineInformation()
    merged = pv.MergeBlocks(Input=reader)
    merged.MergePoints = 1
    times = list(reader.TimestepValues) or [0.0]
    if args.max_frames and len(times) > args.max_frames:
        indices = [round(i * (len(times)-1) / max(args.max_frames-1, 1))
                   for i in range(args.max_frames)]
        times = [times[i] for i in indices]
    speed = pv.Calculator(Input=merged)
    speed.ResultArrayName = 'speed'
    factor = 100 if args.units == 'cm/s' else 1
    speed.Function = f'{factor}*mag(velocity)'
    maximum = 0.0
    for t in times:
        speed.UpdatePipeline(t)
        maximum = max(maximum, speed.PointData['speed'].GetRange()[1])
    vmax = args.vmax if args.vmax is not None else maximum
    if not math.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    volume = pv.ResampleToImage(Input=speed)
    volume.SamplingDimensions = [args.resolution] * 3
    volume.UseInputBounds = 1
    surface = pv.ExtractSurface(Input=merged)
    view = pv.CreateView('RenderView')
    view.ViewSize = args.size
    view.UseColorPaletteForBackground = 0
    view.Background = [.015, .025, .055]
    view.OrientationAxesVisibility = 0
    display = pv.Show(volume, view)
    display.Representation = 'Volume'
    pv.ColorBy(display, ('POINTS', 'speed'))
    lut = pv.GetColorTransferFunction('speed')
    lut.ColorSpace = 'RGB'
    lut.RGBPoints = [0, .02, .05, .18, .25*vmax, .08, .65, .85,
                     .5*vmax, .96, .96, .30, .75*vmax, 1., .5, .02,
                     vmax, .95, .12, .015]
    lut.RescaleTransferFunction(0, vmax)
    opacity = pv.GetOpacityTransferFunction('speed')
    opacity.Points = [0, 0, .5, 0, .02*vmax, .002, .5, 0,
                      .25*vmax, .04, .5, 0, .5*vmax, .16, .5, 0,
                      vmax, .6, .5, 0]
    opacity.RescaleTransferFunction(0, vmax)
    speed.UpdatePipeline(times[0])
    bounds = speed.GetDataInformation().GetBounds()
    diagonal = math.sqrt(sum((bounds[2*i+1]-bounds[2*i])**2 for i in range(3)))
    display.ScalarOpacityUnitDistance = diagonal / 60
    skin = pv.Show(surface, view)
    skin.Representation = 'Surface'
    pv.ColorBy(skin, None)
    skin.DiffuseColor = [.8, .85, .95]
    skin.Opacity = .06
    skin.Specular = .1
    display.SetScalarBarVisibility(view, True)
    bar = pv.GetScalarBar(lut, view)
    bar.Title = 'Speed'
    bar.ComponentTitle = f'({args.units})'
    bar.TitleColor = [.95, .95, .98]
    bar.LabelColor = [.95, .95, .98]
    bar.WindowLocation = 'Lower Right Corner'
    title = pv.Text(Text=args.title)
    title_display = pv.Show(title, view)
    title_display.WindowLocation = 'Upper Left Corner'
    title_display.FontSize = 18
    title_display.Color = [.95, .95, .98]
    time = pv.AnnotateTimeFilter(Input=reader)
    time.Format = 't = {time:.3f} s'
    time_display = pv.Show(time, view)
    time_display.WindowLocation = 'Upper Right Corner'
    time_display.FontSize = 14
    time_display.Color = [.95, .95, .98]
    center = [(bounds[2*i]+bounds[2*i+1])/2 for i in range(3)]
    view.CameraFocalPoint = center
    view.CameraPosition = [center[0]+diagonal, center[1]+.65*diagonal,
                           center[2]+.8*diagonal]
    view.CameraViewUp = [0, 0, 1]
    view.ViewTime = times[0]
    pv.ResetCamera(view)
    view.GetActiveCamera().Zoom(1.3)
    pv.SaveState(str(output.with_suffix('.pvsm')))
    with tempfile.TemporaryDirectory(prefix='atrium-volume-') as tmp:
        for i, t in enumerate(times):
            view.ViewTime = t
            volume.UpdatePipeline(t)
            surface.UpdatePipeline(t)
            pv.Render(view)
            frame = Path(tmp)/f'frame-{i:05d}.png'
            pv.SaveScreenshot(str(frame), view, ImageResolution=args.size)
            if i == len(times)-1:
                shutil.copyfile(frame, output.with_suffix('.png'))
            print(f'Rendered {i+1}/{len(times)}: t={t:.6g} s', flush=True)
        subprocess.run([encoder, '-y', '-v', 'error', '-framerate', str(args.fps),
                        '-i', str(Path(tmp)/'frame-%05d.png'), '-c:v', 'libx264',
                        '-pix_fmt', 'yuv420p', '-crf', '18', str(output)], check=True)
    metadata = dict(source=str(source), units=args.units, velocity_factor=factor,
                    legend_maximum=vmax, selected_times_s=times, fps=args.fps,
                    movie_duration_s=len(times)/args.fps,
                    sampling_dimensions=[args.resolution]*3)
    output.with_suffix('.json').write_text(json.dumps(metadata, indent=2))
    print('Movie:', output)


if __name__ == '__main__':
    main()
