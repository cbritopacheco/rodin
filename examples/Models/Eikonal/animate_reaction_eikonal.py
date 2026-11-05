"""Compare RE+ voltage propagation on identical isotropic and anisotropic meshes."""
from pathlib import Path
import argparse
import csv
import json
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from plot_style import apply_style, turbo_colors


def read_case(directory):
    """Read native mesh, exported nodal snapshots, activation times and run inputs."""
    with (directory / 'ReactionEikonalParameters.csv').open() as source:
        parameters = {row['parameter']: float(row['value'])
                      for row in csv.DictReader(source)}
    table = np.genfromtxt(directory / 'ReactionEikonalActivation.csv',
                          delimiter=',', names=True)
    with h5py.File(directory / 'ReactionEikonalVoltage.mesh.h5') as source:
        vertices = source['/Mesh/Geometry/Vertices'][:]
        topology = source['/Mesh/XDMF/Topology'][:].reshape(-1, 4)
    if not np.all(topology[:, 0] == 4):
        raise ValueError('Expected native triangles')
    if not np.allclose(vertices, np.column_stack([table['x_mm'], table['y_mm']]),
                       rtol=0, atol=1e-12):
        raise ValueError('CSV vertex order differs from the exported mesh')
    frames = ET.parse(directory / 'ReactionEikonalVoltage.xdmf').findall(
        './/Grid[@GridType="Uniform"]')
    times, values = [], []
    for frame in frames:
        times.append(float(frame.find('Time').get('Value')))
        reference = frame.find('Attribute[@Name="voltage_mV"]/DataItem').text.strip()
        filename, dataset = reference.split(':', 1)
        with h5py.File(directory / filename) as source:
            values.append(source[dataset][:])
    times, values = np.asarray(times), np.asarray(values)
    if (len(times) < 2 or not np.all(np.diff(times) > 0)
            or values.shape != (len(times), len(vertices))
            or not np.isfinite(values).all()
            or not np.isfinite(table['replus_lat_ms']).all()):
        raise ValueError('Invalid or incomplete voltage/activation data')
    return parameters, table, vertices, topology[:, 1:], times, values


def repolarization_times(case):
    """First post-peak 90% recovery, relative to each node's initial voltage."""
    values, times = case[5], case[4]
    peaks = values.argmax(axis=0)
    columns = np.arange(values.shape[1])
    threshold = values[0] + 0.1 * (values[peaks, columns] - values[0])
    crossings = ((np.arange(len(times))[:, None] > peaks)
                 & (values <= threshold))
    valid = crossings.any(axis=0)
    columns = columns[valid]
    upper = crossings.argmax(axis=0)[valid]
    lower = upper - 1
    result = np.full(values.shape[1], np.inf)
    fraction = ((threshold[valid] - values[lower, columns])
                / (values[upper, columns] - values[lower, columns]))
    result[valid] = times[lower] + fraction * (times[upper] - times[lower])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('isotropic', type=Path)
    parser.add_argument('anisotropic', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--end-ms', type=float, default=50)
    parser.add_argument('--duration', type=float, default=10,
                        help='Playback duration in seconds; physical time remains labelled')
    parser.add_argument('--repolarization', action='store_true',
                        help='Show 90% repolarization and export matched-probe traces')
    args = parser.parse_args()
    ffmpeg = shutil.which('ffmpeg')
    if ffmpeg is None:
        raise RuntimeError('ffmpeg with libvpx-vp9 and libx264 is required')
    cases = [read_case(args.isotropic), read_case(args.anisotropic)]
    first, second = cases
    for key in ['intervals', 'dt_ms', 'conditioning_beats', 'axis_aligned_geometry', 'fiber_x', 'fiber_y', 'source_x_mm', 'source_y_mm', 'source_radius_mm']:
        if first[0][key] != second[0][key]:
            raise ValueError(f'Comparison differs in {key}')
    if (first[0]['axis_aligned_geometry'] != 1
            or first[0]['nominal_velocity_ratio'] != 1
            or second[0]['nominal_velocity_ratio'] <= 1):
        raise ValueError('Use axis-aligned uniform-grid runs with isotropic and anisotropic tensors')
    if (first[2].shape != second[2].shape
            or not np.allclose(first[2], second[2], rtol=0, atol=1e-12)
            or not np.array_equal(first[3], second[3])
            or not np.array_equal(first[1]['stimulated'], second[1]['stimulated'])
            or not np.array_equal(first[4], second[4])):
        raise ValueError('Comparison requires the same mesh, sources and output times')
    if not (0 < args.end_ms <= min(case[4][-1] for case in cases)
            and 0 < args.duration):
        raise ValueError('Requested duration/end time is outside the exported series')
    repolarization = [repolarization_times(case) for case in cases] if args.repolarization else []
    coordinates = first[2]
    nx = int(first[0]['intervals']) + 1
    if (len(coordinates) != nx**2
            or not np.allclose(np.unique(coordinates[:, 0]), np.linspace(0, 10, nx))
            or not np.allclose(np.unique(coordinates[:, 1]), np.linspace(0, 10, nx))):
        raise ValueError('Expected an axis-aligned uniform 10 mm square grid')
    args.output.mkdir(parents=True, exist_ok=True)
    apply_style()
    plt.rcParams['animation.ffmpeg_path'] = ffmpeg
    figure, axes = plt.subplots(1, 2, figsize=(10, 7), dpi=150)
    figure.subplots_adjust(left=0.07, right=0.85, bottom=0.16, top=0.79, wspace=0.18)
    title = ('RE+ depolarization and repolarization' if args.repolarization
             else 'RE+ propagation on the same tissue sheet')
    figure.suptitle(title, fontsize=17, y=0.98)
    source_x, source_y = first[0]['source_x_mm'], first[0]['source_y_mm']
    source_point = np.array([source_x, source_y])
    figure.text(0.5, 0.905,
                rf'Source at $({source_x:g},{source_y:g})\,\mathrm{{mm}}$ • same mesh, kinetics and time step',
                ha='center', fontsize=11)
    surfaces, percentages = [], []
    for axis, case, title in zip(axes, cases, ['Isotropic', 'Anisotropic']):
        parameters, table, vertices, triangles, times, values = case
        ratio = parameters['nominal_velocity_ratio']
        axis.set_title(title + rf'  ($\sigma_l/\sigma_t={ratio**2:g}$)', fontsize=14, pad=26)
        axis.text(0.5, 1.025,
                  rf'$c_l={parameters["longitudinal_velocity_mm_per_ms"]:.3f},\;'
                  rf'c_t={parameters["transverse_velocity_mm_per_ms"]:.3f}'
                  r'\;[\mathrm{mm/ms}]$', transform=axis.transAxes,
                  ha='center', fontsize=11)
        triangulation = mtri.Triangulation(vertices[:, 0], vertices[:, 1], triangles)
        surface = axis.tripcolor(triangulation, values[0], shading='gouraud',
                                 cmap='turbo', vmin=-85, vmax=35, rasterized=True)
        surfaces.append(surface)
        axis.scatter(*source_point, marker='*', s=75, facecolor='white',
                     edgecolor='#303b49', linewidth=1, zorder=5)
        axis.set(xlabel=r'$x\;[\mathrm{mm}]$', ylabel=r'$y\;[\mathrm{mm}]$', aspect='equal')
        axis.set_xlim(vertices[:, 0].min() - 0.5, vertices[:, 0].max() + 0.5)
        axis.set_ylim(vertices[:, 1].min() - 0.5, vertices[:, 1].max() + 0.5)
        box = axis.get_position()
        percentages.append(figure.text(box.x0 + box.width / 2, 0.12, '',
                                       ha='center', fontsize=10))
    fiber = np.array([first[0]['fiber_x'], first[0]['fiber_y']])
    fiber_origin = np.array([0.8, 1.0])
    axes[1].annotate(r'$\mathbf{f}$', xy=fiber_origin + 1.5 * fiber, xytext=fiber_origin,
                     arrowprops=dict(arrowstyle='->', color='#303b49', lw=1.8),
                     fontsize=14, color='#303b49')
    colorbar = figure.colorbar(surfaces[0], cax=figure.add_axes([0.9, 0.25, 0.022, 0.48]))
    colorbar.set_label(r'$V\;[\mathrm{mV}]$', fontsize=12)
    clock = figure.text(0.5, 0.045 if args.repolarization else 0.075, '', ha='center', fontsize=15)
    figure.text(0.5, 0.008 if args.repolarization else 0.025,
                rf'Playback: ${args.end_ms / args.duration:g}\,\mathrm{{ms}}$ per second'
                '  |  Linear interpolation between exported snapshots',
                ha='center', fontsize=10)

    def update(time):
        for index, (case, surface, percentage) in enumerate(zip(cases, surfaces, percentages)):
            times, values = case[4:]
            upper = min(max(1, np.searchsorted(times, time)), len(times) - 1)
            weight = (time - times[upper - 1]) / (times[upper] - times[upper - 1])
            surface.set_array((1 - weight) * values[upper - 1] + weight * values[upper])
            fraction = np.count_nonzero(case[1]['replus_lat_ms'] <= time) / len(case[1])
            label = f'Activated: {100 * fraction:5.1f}%' if args.repolarization else f'Nodes activated: {100 * fraction:5.1f}%'
            if args.repolarization:
                recovered = np.count_nonzero(repolarization[index] <= time) / len(case[1])
                label += f' | 90% repolarized: {100 * recovered:5.1f}%'
            percentage.set_text(label)
        clock.set_text(rf'$t={time:05.2f}\;\mathrm{{ms}}$')

    poster_time = min(300 if args.repolarization else 20, args.end_ms)
    update(poster_time)
    figure.savefig(args.output / 'REPlusComparison.png', transparent=True, dpi=150)
    fps = 30
    frame_times = np.linspace(0, args.end_ms, round(args.duration * fps) + 1)
    writer = animation.FFMpegWriter(fps=fps, codec='libvpx-vp9',
                                   extra_args=['-pix_fmt', 'yuva420p', '-auto-alt-ref', '0',
                                               '-crf', '28', '-b:v', '0'])
    webm = args.output / 'REPlusComparison.webm'
    with writer.saving(figure, str(webm), dpi=150):
        for time in frame_times:
            update(time)
            writer.grab_frame(transparent=True)
    plt.close(figure)
    # Flatten the alpha-capable WebM onto white for broad MP4 player compatibility.
    subprocess.run([ffmpeg, '-hide_banner', '-loglevel', 'error', '-y',
                    '-c:v', 'libvpx-vp9', '-i', str(webm), '-filter_complex',
                    f'color=c=white:s=1500x1050:r={fps}[bg];'
                    '[bg][0:v]overlay=shortest=1:format=auto,format=yuv420p',
                    '-c:v', 'libx264', '-crf', '18', '-movflags', '+faststart',
                    str(args.output / 'REPlusComparison.mp4')], check=True)
    summary = dict(fps=fps, frames=len(frame_times), physical_end_ms=args.end_ms,
                   playback_ms_per_second=args.end_ms / args.duration,
                   identical_mesh=True, identical_sources=True,
                   spatial_display='native P1 triangles', temporal_display='linear interpolation',
                   cases=[dict(parameters=case[0], latest_activation_ms=float(case[1]['replus_lat_ms'].max()))
                          for case in cases])
    if args.repolarization:
        summary['repolarization90_definition'] = 'First post-peak crossing of V_initial + 0.1*(V_peak_exported - V_initial)'
        summary['apd90_definition'] = 'Repolarization90 time minus steepest-upstroke LAT; unweighted nodal statistics'
        summary['repolarization_sampling'] = 'Peaks and linearly interpolated crossings from exported snapshots, approximately 1 ms cadence'
        probe = np.linalg.norm(first[2] - first[2].mean(axis=0), axis=1).argmin()
        traces, trace_axes = plt.subplots(1, 2, figsize=(11, 4), layout='constrained')
        for case, recovery, color, label, item in zip(
                cases, repolarization, turbo_colors(2), ['Isotropic', 'Anisotropic'], summary['cases']):
            lat = case[1]['replus_lat_ms']
            apd = recovery - lat
            finite = np.isfinite(recovery)
            item['repolarized_nodes'] = int(finite.sum())
            item['latest_repolarization90_ms'] = float(recovery[finite].max()) if finite.any() else None
            item['mean_nodal_apd90_ms'] = float(apd[finite].mean()) if finite.any() else None
            item['apd90_range_ms'] = [float(apd[finite].min()), float(apd[finite].max())] if finite.any() else None
            item['probe'] = dict(vertex=int(probe), coordinates_mm=case[2][probe].tolist(),
                                 lat_ms=float(lat[probe]),
                                 apd90_ms=float(apd[probe]) if np.isfinite(apd[probe]) else None)
            for axis, shift in zip(trace_axes, [0, lat[probe]]):
                axis.plot(case[4] - shift, case[5][:, probe], label=label, color=color)
                if np.isfinite(recovery[probe]):
                    voltage90 = np.interp(recovery[probe], case[4], case[5][:, probe])
                    axis.plot(recovery[probe] - shift, voltage90, 'o', color=color, markersize=4)
                axis.set(ylabel=r'$V\;[\mathrm{mV}]$', ylim=(-85, 35))
                axis.legend()
        trace_axes[0].set(title='Same central probe: absolute time', xlabel=r'$t\;[\mathrm{ms}]$')
        trace_axes[1].set(title='Aligned by local activation time', xlabel=r'$t-t_{\mathrm{LAT}}\;[\mathrm{ms}]$')
        traces.savefig(args.output / 'REPlusRepolarizationTraces.png', dpi=180, transparent=True)
        traces.savefig(args.output / 'REPlusRepolarizationTraces.svg', transparent=True)
        plt.close(traces)
    (args.output / 'REPlusComparison.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
