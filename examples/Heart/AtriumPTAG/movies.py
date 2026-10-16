#!/usr/bin/env pvpython
# -*- coding: utf-8 -*-
"""
atrium_ptag_movies.py — videos y snapshots del caso AtriumPTAG (Rodin) para
presentar a público clínico. Se ejecuta con el pvpython de ParaView (5.10+):

    pvpython atrium_ptag_movies.py AtriumPTAG.xdmf --out movies

    # sólo algunas escenas, más rápido (menos frames, menor resolución):
    pvpython atrium_ptag_movies.py AtriumPTAG.xdmf --scenes streamlines,thrombin_wall \
        --res 1280x720 --every 2 --volume-dims 120

Escenas (--scenes, separadas por coma; "all" = todas):
    streamlines           líneas de corriente (tubos, color = velocidad) dentro
                          de una pared "fantasma" translúcida
    velocity_volume       render volumétrico de la velocidad
    thrombin_wall         trombina sobre el endocardio (pared coloreada) con un
                          resplandor volumétrico de trombina en el interior
    fibrin_volume         fibrina (F = G0 - G) como volumen, "el coágulo"
    prothrombin_depletion consumo de protrombina, (P0 - P)/P0 en %

Salida (en --out):
    <escena>/frame_000000.png ...   fotogramas
    <escena>.mp4                    video (si hay ffmpeg en el PATH)
    snapshots/<escena>_{first,mid,last}.png  tres instantáneas en alta resolución

En macOS el pvpython está dentro de la app, p. ej.
    /Applications/ParaView-5.12.0.app/Contents/bin/pvpython
"""
import argparse
import math
import os
import shutil
import subprocess
import sys

from paraview.simple import *  # noqa: F401,F403

paraview.simple._DisableFirstRenderCameraReset()

# --------------------------------------------------------------------------
# Parámetros de las escenas: nombre del campo que se pinta, unidad mostrada,
# factor de conversión, preset de color, texto del título.
# --------------------------------------------------------------------------
# Paletas propias (fracción del rango -> RGB). Los campos de especies parten
# en un color visible, no en negro, porque se dibujan sobre fondo oscuro.
LUT_TISSUE_HOT = [(0.0, [0.88, 0.84, 0.82]), (0.25, [0.97, 0.75, 0.45]),
                  (0.55, [0.92, 0.25, 0.10]), (0.8, [0.75, 0.02, 0.05]),
                  (1.0, [1.0, 0.95, 0.5])]
LUT_EMBER = [(0.0, [0.35, 0.0, 0.12]), (0.3, [0.85, 0.05, 0.12]),
             (0.65, [1.0, 0.5, 0.1]), (1.0, [1.0, 0.97, 0.7])]
LUT_GOLD = [(0.0, [0.4, 0.15, 0.0]), (0.35, [0.95, 0.55, 0.05]),
            (0.7, [1.0, 0.85, 0.3]), (1.0, [1.0, 1.0, 0.85])]
LUT_ICE = [(0.0, [0.05, 0.15, 0.45]), (0.35, [0.1, 0.55, 0.9]),
           (0.7, [0.55, 0.9, 1.0]), (1.0, [1.0, 1.0, 1.0])]
LUT_FLOW = [(0.0, [0.05, 0.1, 0.5]), (0.3, [0.1, 0.8, 0.9]),
            (0.6, [0.95, 0.9, 0.2]), (1.0, [1.0, 0.3, 0.1])]

SCENES = {
    "streamlines": dict(
        title="Flujo sanguíneo en la aurícula izquierda",
        field="speed", label="Velocidad", unit="cm/s", source="velocity",
        preset="Turbo", lut=None, species=False, range_scale=0.8),
    "velocity_volume": dict(
        title="Velocidad de la sangre (render volumétrico)",
        field="speed", label="Velocidad", unit="cm/s", source="velocity",
        preset=None, lut=LUT_FLOW, species=False, range_scale=0.5,
        opacity=[(0.0, 0.0), (0.06, 0.0), (0.25, 0.2), (0.6, 0.7), (1.0, 1.0)]),
    "thrombin_wall": dict(
        title="Trombina sobre el endocardio",
        field="thrombin_nM", label="Trombina", unit="nM", source="thrombin",
        preset=None, lut=LUT_TISSUE_HOT, species=True),
    "thrombin_volume": dict(
        title="Trombina en la cavidad auricular",
        field="thrombin_nM", label="Trombina", unit="nM", source="thrombin",
        preset=None, lut=LUT_GOLD, species=True,
        opacity=[(0.0, 0.0), (0.05, 0.0), (0.25, 0.4), (0.6, 0.85), (1.0, 1.0)]),
    "fibrin_volume": dict(
        title="Formación de fibrina (coágulo)",
        field="fibrin_uM", label="Fibrina", unit="µM", source="fibrin",
        preset=None, lut=LUT_EMBER, species=True,
        opacity=[(0.0, 0.0), (0.05, 0.0), (0.25, 0.5), (0.6, 0.9), (1.0, 1.0)]),
    "prothrombin_depletion": dict(
        title="Consumo de protrombina",
        field="prothrombin_depletion", label="Protrombina consumida", unit="%",
        source="prothrombin", preset=None, lut=LUT_ICE, species=True,
        opacity=[(0.0, 0.0), (0.08, 0.0), (0.3, 0.3), (0.65, 0.75), (1.0, 1.0)]),
}

DARK_TOP = [0.02, 0.03, 0.07]
DARK_BOTTOM = [0.10, 0.12, 0.20]
WALL_COLOR = [0.86, 0.80, 0.78]       # tejido pálido
TEXT_COLOR = [0.95, 0.95, 0.97]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xdmf", help="AtriumPTAG.xdmf escrito por Rodin")
    p.add_argument("--out", default="movies", help="carpeta de salida")
    p.add_argument("--scenes", default="all",
                   help="lista separada por comas: " + ",".join(SCENES) + " | all")
    p.add_argument("--res", default="1920x1080", help="resolución WxH del video")
    p.add_argument("--snap-res", default="3840x2160", help="resolución de los snapshots")
    p.add_argument("--fps", type=int, default=12)
    p.add_argument("--every", type=int, default=1, help="usar 1 de cada N pasos de tiempo")
    p.add_argument("--from-time", type=float, default=None)
    p.add_argument("--to-time", type=float, default=None)
    p.add_argument("--max-frames", type=int, default=0, help="0 = sin límite")
    p.add_argument("--orbit", type=float, default=90.0,
                   help="grados de giro de la cámara a lo largo del video (0 = fija)")
    p.add_argument("--view-dir", default="1,-0.7,0.9",
                   help="dirección desde la que mira la cámara (x,y,z)")
    p.add_argument("--view-up", default="0,0,1", help="vector 'arriba' de la cámara")
    p.add_argument("--zoom", type=float, default=1.0)
    p.add_argument("--wall-opacity", type=float, default=0.12,
                   help="opacidad de la pared fantasma en las escenas de volumen/líneas")
    p.add_argument("--seeds", type=int, default=500, help="semillas de las líneas de corriente")
    p.add_argument("--tube", type=float, default=0.0045,
                   help="radio de los tubos, como fracción de la diagonal del dominio")
    p.add_argument("--volume-dims", type=int, default=160,
                   help="resolución del remuestreo para el render volumétrico (n^3)")
    p.add_argument("--range-samples", type=int, default=4,
                   help="pasos de tiempo usados para fijar el rango de color (0 = todos)")
    p.add_argument("--no-video", action="store_true", help="sólo fotogramas y snapshots")
    p.add_argument("--light", action="store_true", help="fondo claro en vez de oscuro")
    return p.parse_args()


def vec(s):
    v = [float(x) for x in s.split(",")]
    n = math.sqrt(sum(x * x for x in v))
    return [x / n for x in v]


def log(*a):
    print("[movies]", *a)
    sys.stdout.flush()


# --------------------------------------------------------------------------
# Lectura
# --------------------------------------------------------------------------
def load(path, arrays):
    # El lector Xdmf3 resuelve los .h5 relativos a la carpeta del .xdmf; con una
    # ruta sin directorio concatena mal, así que siempre se pasa la absoluta.
    reader = Xdmf3ReaderS(FileName=[os.path.abspath(path)])
    reader.UpdatePipelineInformation()
    available = list(reader.PointArrays)
    wanted = [a for a in arrays if a in available]
    missing = [a for a in arrays if a not in available]
    if missing:
        log("campos no encontrados en el XDMF:", missing, "— disponibles:", available)
    reader.PointArrays = wanted
    merged = MergeBlocks(Input=reader)          # los bloques por rango -> una malla
    clean = CleantoGrid(Input=merged)           # une los nodos duplicados de interfaz
    return reader, clean


def add_fields(src, needed):
    """Campos derivados con unidades clínicas; devuelve el último filtro."""
    cur = src
    if "speed" in needed:
        cur = Calculator(Input=cur, ResultArrayName="speed", Function="mag(velocity)*100")
    if "thrombin_nM" in needed:
        cur = Calculator(Input=cur, ResultArrayName="thrombin_nM", Function="thrombin*1e6")
    if "fibrin_uM" in needed:
        cur = Calculator(Input=cur, ResultArrayName="fibrin_uM", Function="fibrin*1e3")
    if "prothrombin_depletion" in needed:
        p0 = needed["prothrombin_depletion"]
        cur = Calculator(Input=cur, ResultArrayName="prothrombin_depletion",
                         Function="100*(%.9g-prothrombin)/%.9g" % (p0, p0))
    return cur


def field_range(src, field, times, samples, component=None):
    """(min, max) del campo sobre varios pasos de tiempo (escala de color fija)."""
    pick = times if samples <= 0 or samples >= len(times) else \
        [times[int(round(i * (len(times) - 1) / (samples - 1)))] for i in range(samples)]
    lo, hi = float("inf"), float("-inf")
    for t in pick:
        src.UpdatePipeline(t)
        arr = src.PointData[field]
        r = arr.GetRange(-1) if component == -1 else arr.GetRange()
        lo, hi = min(lo, r[0]), max(hi, r[1])
    return lo, hi


def first_time_with_signal(src, field, times, threshold=0.0):
    for t in times:
        src.UpdatePipeline(t)
        if src.PointData[field].GetRange()[1] > threshold:
            return t
    return None


# --------------------------------------------------------------------------
# Vista
# --------------------------------------------------------------------------
def make_view(res, light):
    view = CreateView("RenderView")
    view.ViewSize = res
    view.OrientationAxesVisibility = 0
    view.UseColorPaletteForBackground = 0
    if light:
        view.BackgroundColorMode = "Single Color"
        view.Background = [1, 1, 1]
    else:
        view.BackgroundColorMode = "Gradient"
        view.Background = DARK_BOTTOM
        view.Background2 = DARK_TOP
    try:
        view.EnableRayTracing = 0
    except Exception:
        pass
    return view


def frame_camera(view, bounds, view_dir, view_up, zoom):
    cx = 0.5 * (bounds[0] + bounds[1])
    cy = 0.5 * (bounds[2] + bounds[3])
    cz = 0.5 * (bounds[4] + bounds[5])
    diag = math.sqrt((bounds[1] - bounds[0]) ** 2 + (bounds[3] - bounds[2]) ** 2 +
                     (bounds[5] - bounds[4]) ** 2)
    dist = 1.35 * diag / zoom
    view.CameraFocalPoint = [cx, cy, cz]
    view.CameraPosition = [cx + view_dir[0] * dist, cy + view_dir[1] * dist,
                           cz + view_dir[2] * dist]
    view.CameraViewUp = view_up
    view.CameraViewAngle = 30
    view.CameraParallelProjection = 0
    return diag


def text_overlay(view, text, location, size):
    src = Text(Text=text)
    disp = Show(src, view)
    disp.WindowLocation = location
    disp.FontSize = size
    disp.Color = TEXT_COLOR if view.BackgroundColorMode == "Gradient" else [0.1, 0.1, 0.1]
    disp.FontFamily = "Arial"
    disp.Bold = 1 if size > 20 else 0
    return src, disp


def style_scalar_bar(rep, view, lut, label, unit, dark):
    rep.SetScalarBarVisibility(view, True)
    bar = GetScalarBar(lut, view)
    bar.Title = label
    bar.ComponentTitle = "[%s]" % unit
    color = TEXT_COLOR if dark else [0.1, 0.1, 0.1]
    bar.TitleColor = color
    bar.LabelColor = color
    bar.TitleFontSize = 18
    bar.LabelFontSize = 16
    bar.WindowLocation = "Lower Right Corner"
    bar.ScalarBarLength = 0.35
    bar.ScalarBarThickness = 14
    bar.RangeLabelFormat = "%-#6.3g"
    bar.AutomaticLabelFormat = 0
    bar.LabelFormat = "%-#6.3g"
    bar.DrawTickMarks = 1


def ghost_wall(surface, view, opacity):
    """La pared como un cascarón translúcido: se ve la forma, no tapa el interior."""
    rep = Show(surface, view)
    rep.SetRepresentationType("Surface")
    ColorBy(rep, None)
    rep.DiffuseColor = WALL_COLOR
    rep.AmbientColor = WALL_COLOR
    rep.Opacity = opacity
    rep.Specular = 0.25
    rep.SpecularPower = 40.0
    rep.Interpolation = "Gouraud"
    rep.BackfaceRepresentation = "Cull Backface"
    return rep


def apply_colors(field, cfg, lut_range):
    lut = GetColorTransferFunction(field)
    lo, hi = lut_range
    if cfg.get("lut"):
        pts = []
        for frac, rgb in cfg["lut"]:
            pts += [lo + frac * (hi - lo)] + list(rgb)
        lut.RGBPoints = pts
        lut.ColorSpace = "Lab"
    else:
        lut.ApplyPreset(cfg["preset"], True)
        lut.RescaleTransferFunction(lo, hi)
    lut.ScalarRangeInitialized = 1
    return lut


def volume_rep(resampled, view, field, lut_range, cfg, opacity_curve, diag):
    rep = Show(resampled, view)
    rep.SetRepresentationType("Volume")
    ColorBy(rep, ("POINTS", field))
    lut = apply_colors(field, cfg, lut_range)
    pwf = GetOpacityTransferFunction(field)
    lo, hi = lut_range
    pts = []
    for frac, alpha in opacity_curve:
        pts += [lo + frac * (hi - lo), alpha, 0.5, 0.0]
    pwf.Points = pts
    pwf.ScalarRangeInitialized = 1
    rep.Shade = 1
    # Opacidad acumulada por unidad de longitud: con el valor por defecto el
    # volumen queda demasiado tenue sobre el fondo oscuro.
    rep.ScalarOpacityUnitDistance = 0.5 * rep.ScalarOpacityUnitDistance
    try:
        rep.VolumeRenderingMode = "Smart"
    except Exception:
        pass
    return rep, lut


# --------------------------------------------------------------------------
# Escenas
# --------------------------------------------------------------------------
def build_scene(name, cfg, args, data, view, bounds, diag, ranges, dark):
    """Devuelve (lut, representaciones a actualizar) para la escena."""
    surface = ExtractSurface(Input=data)
    field = cfg["field"]
    lo, hi = ranges[field]
    hi = lo + cfg.get("range_scale", 1.0) * (hi - lo)   # saturar los chorros de entrada

    if name == "streamlines":
        ghost_wall(surface, view, args.wall_opacity)
        tracer = StreamTracer(Input=data, SeedType="Point Cloud")
        tracer.Vectors = ["POINTS", "velocity"]
        tracer.SeedType.Center = [0.5 * (bounds[0] + bounds[1]), 0.5 * (bounds[2] + bounds[3]),
                                  0.5 * (bounds[4] + bounds[5])]
        tracer.SeedType.Radius = 0.42 * diag
        tracer.SeedType.NumberOfPoints = args.seeds
        tracer.IntegrationDirection = "BOTH"
        tracer.MaximumStreamlineLength = 2.0 * diag
        tracer.MaximumSteps = 4000
        tracer.ComputeVorticity = 0
        tube = Tube(Input=tracer, Radius=args.tube * diag, NumberofSides=10, Capping=1)
        rep = Show(tube, view)
        ColorBy(rep, ("POINTS", field))
        lut = apply_colors(field, cfg, (lo, hi))
        rep.Specular = 0.5
        rep.SpecularPower = 30.0
        style_scalar_bar(rep, view, lut, cfg["label"], cfg["unit"], dark)
        return lut, [tracer]

    if "opacity" in cfg:
        # Render volumétrico dentro de la pared fantasma; el campo se remuestrea
        # a una imagen regular (mucho más rápido que el ray casting sobre tetraedros).
        ghost_wall(surface, view, args.wall_opacity)
        n = args.volume_dims
        resampled = ResampleToImage(Input=data, SamplingDimensions=[n, n, n], UseInputBounds=1)
        rep, lut = volume_rep(resampled, view, field, (lo, hi), cfg, cfg["opacity"], diag)
        style_scalar_bar(rep, view, lut, cfg["label"], cfg["unit"], dark)
        return lut, [resampled]

    if name == "thrombin_wall":
        # Pared opaca coloreada por la trombina de los nodos del endocardio; el
        # cero es color de tejido, así la aurícula sigue pareciendo una aurícula.
        rep = Show(surface, view)
        rep.SetRepresentationType("Surface")
        ColorBy(rep, ("POINTS", field))
        lut = apply_colors(field, cfg, (lo, hi))
        rep.Specular = 0.35
        rep.SpecularPower = 40.0
        rep.Interpolation = "Gouraud"
        style_scalar_bar(rep, view, lut, cfg["label"], cfg["unit"], dark)
        return lut, [surface]

    raise ValueError(name)


def encode(frames_dir, out_mp4, fps):
    ffmpeg = shutil.which("ffmpeg")
    if not ffmpeg:
        log("ffmpeg no está en el PATH; los fotogramas quedan en", frames_dir)
        log("para el video:  ffmpeg -framerate %d -i %s/frame_%%06d.png "
            "-c:v libx264 -pix_fmt yuv420p -crf 18 %s" % (fps, frames_dir, out_mp4))
        return
    cmd = [ffmpeg, "-y", "-loglevel", "error", "-framerate", str(fps),
           "-i", os.path.join(frames_dir, "frame_%06d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", out_mp4]
    subprocess.run(cmd, check=False)
    log("video:", out_mp4)


def main():
    args = parse_args()
    res = [int(x) for x in args.res.lower().split("x")]
    snap_res = [int(x) for x in args.snap_res.lower().split("x")]
    view_dir = vec(args.view_dir)
    view_up = vec(args.view_up)
    dark = not args.light
    scenes = list(SCENES) if args.scenes == "all" else [s.strip() for s in args.scenes.split(",")]
    for s in scenes:
        if s not in SCENES:
            sys.exit("escena desconocida: %s (opciones: %s)" % (s, ", ".join(SCENES)))
    os.makedirs(args.out, exist_ok=True)
    os.makedirs(os.path.join(args.out, "snapshots"), exist_ok=True)

    needed_sources = sorted({SCENES[s]["source"] for s in scenes})
    reader, clean = load(args.xdmf, needed_sources)
    times = list(reader.TimestepValues) or [0.0]
    log("pasos de tiempo:", len(times), " t =", times[0], "...", times[-1])

    # Rango de protrombina de referencia: su valor de entrada (máximo al inicio).
    needed = {}
    for s in scenes:
        f = SCENES[s]["field"]
        if f == "prothrombin_depletion":
            clean.UpdatePipeline(times[0])
            needed[f] = clean.PointData["prothrombin"].GetRange()[1]
        else:
            needed[f] = None
    data = add_fields(clean, needed)
    data.UpdatePipeline(times[0])
    bounds = data.GetDataInformation().GetBounds()
    log("nodos:", data.GetDataInformation().GetNumberOfPoints(), " bounds:", bounds)

    # Rangos de color fijos en el tiempo, por campo.
    ranges = {}
    for s in scenes:
        f = SCENES[s]["field"]
        if f in ranges:
            continue
        lo, hi = field_range(data, f, times, args.range_samples)
        lo = 0.0
        hi = max(hi, 1e-12)
        if f == "prothrombin_depletion":
            hi = min(hi, 100.0)
        ranges[f] = (lo, hi)
        log("rango de %s: [%g, %g]" % (f, lo, hi))

    for name in scenes:
        cfg = SCENES[name]
        # Las especies arrancan después del calentamiento: se omiten los pasos
        # en que todavía son nulas, salvo que el usuario fije --from-time.
        t0 = args.from_time
        if t0 is None and cfg["species"]:
            t0 = first_time_with_signal(data, cfg["field"], times)
            if t0 is None:
                log("la escena %s no tiene señal (campo %s nulo); se omite" % (name, cfg["field"]))
                continue
        sel = [t for t in times if (t0 is None or t >= t0 - 1e-12) and
               (args.to_time is None or t <= args.to_time + 1e-12)]
        sel = sel[::max(1, args.every)]
        if args.max_frames > 0:
            sel = sel[:args.max_frames]
        if not sel:
            log("sin pasos de tiempo para", name)
            continue

        log("escena %s: %d fotogramas" % (name, len(sel)))
        view = make_view(res, args.light)
        diag = frame_camera(view, bounds, view_dir, view_up, args.zoom)
        lut, updaters = build_scene(name, cfg, args, data, view, bounds, diag, ranges, dark)
        text_overlay(view, cfg["title"], "Upper Left Corner", 26)
        _, tdisp = text_overlay(view, "", "Upper Right Corner", 20)
        _, _ = text_overlay(view, "CFDLab USACH · Rodin", "Lower Left Corner", 14)

        frames_dir = os.path.join(args.out, name)
        os.makedirs(frames_dir, exist_ok=True)
        scene = GetAnimationScene()
        scene.UpdateAnimationUsingDataTimeSteps()
        period = None
        d_az = args.orbit / max(1, len(sel) - 1) if len(sel) > 1 else 0.0
        keyframes = {0: "first", len(sel) // 2: "mid", len(sel) - 1: "last"}
        for i, t in enumerate(sel):
            scene.AnimationTime = t
            view.ViewTime = t
            tdisp.Input.Text = "t = %.2f s" % t
            if i > 0 and d_az:
                GetActiveCamera().Azimuth(d_az)
            Render(view)
            SaveScreenshot(os.path.join(frames_dir, "frame_%06d.png" % i), view,
                           ImageResolution=res, TransparentBackground=0)
            if i in keyframes:
                SaveScreenshot(os.path.join(args.out, "snapshots",
                                            "%s_%s.png" % (name, keyframes[i])),
                               view, ImageResolution=snap_res, TransparentBackground=0)
            if i % 10 == 0 or i == len(sel) - 1:
                log("  %s: fotograma %d/%d (t = %.3f s)" % (name, i + 1, len(sel), t))
        if not args.no_video:
            encode(frames_dir, os.path.join(args.out, name + ".mp4"), args.fps)
        Delete(view)
        del view

    log("listo. Salida en", os.path.abspath(args.out))


if __name__ == "__main__":
    main()
