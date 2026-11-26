"""Shared helpers for the PTAGP validation tests (run, read, plot).

Each test script imports this module. The executable is found from
--exe, the environment variable PTAGP_EXE, or the default build path
<repo>/build/examples/Heart/Thrombosis/PTAGP/PTAGP. MPI launcher from
--mpirun (default "mpirun -n 1"; add --allow-run-as-root inside containers).
"""
import argparse
import concurrent.futures as cf
import csv
import os
import shlex
import subprocess
import sys
import xml.etree.ElementTree as ET

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
DEFAULT_EXE = os.path.join(REPO, "build", "examples", "Heart", "Thrombosis", "PTAGP", "PTAGP")
MAKE_MESH = os.path.join(HERE, "make_mesh.py")

# Values of the P-T-A-G constants used in the scripts' derived quantities.
G0 = 7.35e-3      # fibrinogen (mol/m^3)
KF = 2.0e-3       # gel threshold (mol/m^3)
T_TH = 2.48e-6    # thrombin threshold at theta = 1 (mol/m^3)
S_MAX = 2.5e16    # platelet packing (1/m^3)
N0 = 2.5e14       # platelet count (1/m^3)


def parser(description):
    ap = argparse.ArgumentParser(description=description)
    ap.add_argument("--exe", default=os.environ.get("PTAGP_EXE", DEFAULT_EXE))
    ap.add_argument("--mpirun", default=os.environ.get("PTAGP_MPIRUN", "mpirun -n 1"))
    ap.add_argument("--jobs", type=int, default=2, help="simulations run in parallel")
    ap.add_argument("--force", action="store_true", help="rerun cases that already exist")
    ap.add_argument("--plot-only", action="store_true", help="never run, only plot")
    ap.add_argument("--quick", action="store_true", help="coarser, shorter runs")
    ap.add_argument("--out", default=None, help="results folder (default: <test>/results)")
    ap.add_argument("--extra", default=os.environ.get("PTAGP_EXTRA", ""),
                    help='options added to every case, e.g. "pg_fibrin_gate=false pg_M_in_theta=false"')
    return ap


def results_dir(args, script):
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(script)), "results")
    os.makedirs(out, exist_ok=True)
    return out


def make_mesh(path, **kw):
    if os.path.exists(path):
        return path
    cmd = [sys.executable, MAKE_MESH, "-o", path]
    for k, v in kw.items():
        flag = "--" + k.replace("_", "-")
        if v is True:
            cmd.append(flag)
        else:
            cmd += [flag, str(v)]
    subprocess.run(cmd, check=True)
    return path


def run_case(args, outdir, name, options):
    """Run one simulation unless its CSV exists; returns the CSV path."""
    csvp = os.path.join(outdir, name + ".csv")
    if os.path.exists(csvp) and not args.force:
        return csvp
    if args.plot_only:
        raise FileNotFoundError(csvp)
    opts = dict(options)
    for tok in shlex.split(getattr(args, "extra", "") or ""):
        key, val = tok.split("=", 1)
        opts[key] = val
    opts.setdefault("pg_csv", csvp)
    opts.setdefault("pg_xdmf", os.path.join(outdir, name))
    cmd = shlex.split(args.mpirun) + [args.exe]
    for k, v in opts.items():
        cmd += ["-" + k, str(v).lower() if isinstance(v, bool) else str(v)]
    log = os.path.join(outdir, name + ".log")
    with open(log, "w") as f:
        f.write(" ".join(cmd) + "\n")
        f.flush()
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    if r.returncode != 0:
        raise RuntimeError(f"{name} failed, see {log}")
    return csvp


def run_many(args, outdir, cases):
    """cases: dict name -> options. Runs in parallel; returns name -> csv."""
    out = {}
    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(run_case, args, outdir, n, o): n for n, o in cases.items()}
        for fu in cf.as_completed(futs):
            out[futs[fu]] = fu.result()
            print(f"  done: {futs[fu]}", flush=True)
    return out


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def read_xdmf(path, fields):
    """Snapshots of an XDMF series: list of (time, xy, {field: values}).

    Needs h5py; returns [] (with a message) when it is not installed.
    """
    try:
        import h5py
    except ImportError:
        print("  (h5py not installed: spatial profiles skipped)")
        return []
    base = os.path.dirname(path)
    root = ET.parse(path).getroot()
    out = []
    for step in root.iter("Grid"):
        if step.get("CollectionType") != "Spatial":
            continue
        t = float(step.find("Time").get("Value"))
        grid = step.find("Grid")
        geo = grid.find("Geometry/DataItem").text.strip()
        f, ds = geo.split(":")
        with h5py.File(os.path.join(base, f), "r") as h:
            xy = np.array(h[ds])
        vals = {}
        for att in grid.findall("Attribute"):
            name = att.get("Name")
            if name in fields:
                f, ds = att.find("DataItem").text.strip().split(":")
                with h5py.File(os.path.join(base, f), "r") as h:
                    vals[name] = np.array(h[ds])
        out.append((t, xy, vals))
    return out


def first_time(t, mask):
    idx = np.flatnonzero(mask)
    return t[idx[0]] if idx.size else np.nan


def slope(t, y, t0, t1):
    m = (t >= t0) & (t <= t1)
    if m.sum() < 3:
        return np.nan
    return np.polyfit(t[m], y[m], 1)[0]


def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.dpi": 120, "savefig.dpi": 200, "font.size": 9,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "legend.frameon": False,
    })
    return plt


def write_summary(outdir, title, lines):
    path = os.path.join(outdir, "summary.md")
    with open(path, "w") as f:
        f.write(f"# {title}\n\n")
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n  -> {path}")
    return path
