#!/usr/bin/env python3
"""Post-processing of the Biwole et al. (2013) validation (run_biwole.sh).

Solver frame: x along the cavity height (upwards, angle 90 deg), y across the
cavity, y = 0 the cold wall (bottom face), y = W the hot wall (top face).
Paper frame (Fig. 4-7): s = W - y, distance from the hot (left) wall, z = x.

  melt/   front.png, front_tMMM.csv   melt front Y = 0.5 at 20, 40, ..., 160 min (Fig. 5)
          liquid.png                  mean liquid fraction vs time
  steady/ profiles.png, vz_zNNN.csv   vertical velocity vs s at z = 41, 82, 123 mm (Fig. 7.1-3)
          vs_s15.csv                  horizontal velocity (hot -> cold) vs z at s = 15 mm (Fig. 7.4)

    python3 postprocess_biwole.py biwole_runs [--exp DIR]

--exp DIR overlays digitised data: DIR/fig7_1.csv, fig7_2.csv, fig7_3.csv (s_mm, v_m/s),
fig7_4.csv (z_mm, v_m/s), fig5_MMM.csv (s_mm, z_mm of the front), if present.
"""
import argparse
import csv
import os
import xml.etree.ElementTree as ET

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402

W = 0.030   # m, cavity width
H = 0.167   # m, cavity height


def _dataitem(base, item):
    import h5py
    fname, path = item.text.strip().split(":", 1)
    with h5py.File(os.path.join(base, fname), "r") as f:
        return np.array(f[path])


def read_xdmf(path):
    base = os.path.dirname(path)
    root = ET.parse(path).getroot()
    return base, [(float(g.find("Time").get("Value")), g) for g in root.find("Domain/Grid").findall("Grid")]


def load_step(base, step, names):
    parts = []
    for g in step.findall("Grid"):
        topo = g.find("Topology")
        conn = _dataitem(base, topo.find("DataItem")).astype(np.int64).ravel()
        if topo.get("TopologyType", "").lower() == "mixed":
            tri, i = [], 0
            for _ in range(int(topo.get("NumberOfElements"))):
                typ = int(conn[i])
                nv = {4: 3, 5: 4, 1: 1, 2: 2}.get(typ, 3)
                if typ == 2:
                    nv = int(conn[i + 1]) + 1
                if typ == 4:
                    tri.append(conn[i + 1:i + 4])
                i += 1 + nv
            tri = np.array(tri)
        else:
            tri = conn.reshape(-1, 3)
        xy = _dataitem(base, g.find("Geometry/DataItem")).reshape(-1, 2)
        fields = {}
        for a in g.findall("Attribute"):
            if a.get("Name") in names:
                v = _dataitem(base, a.find("DataItem"))
                fields[a.get("Name")] = v.reshape(len(xy), -1) if a.get("Name") == "velocity" else v.ravel()
        parts.append((xy, tri, fields))
    off = np.cumsum([0] + [len(p[0]) for p in parts[:-1]])
    xy = np.vstack([p[0] for p in parts])
    tri = np.vstack([p[1] + o for p, o in zip(parts, off)])
    fields = {n: np.concatenate([p[2][n] for p in parts]) for n in names}
    return xy, tri, fields


def paper_triangulation(xy, tri):
    """(s, z) in mm: s from the hot wall, z upwards."""
    return mtri.Triangulation(1e3 * (W - xy[:, 1]), 1e3 * xy[:, 0], tri)


def read_exp(d, name):
    p = os.path.join(d, name) if d else None
    if not p or not os.path.exists(p):
        return None
    a = np.loadtxt(p, delimiter=",", comments="#")
    return a[:, 0], a[:, 1]


def melt(run, out, exp):
    xdmf = os.path.join(run, "NaturalConvectionPCM.xdmf")
    base, steps = read_xdmf(xdmf)
    minutes = [20, 40, 60, 80, 120, 140, 160]
    fig, axs = plt.subplots(1, len(minutes), figsize=(1.6 * len(minutes), 5.2), sharey=True)
    for ax, mins in zip(axs, minutes):
        t, step = min(steps, key=lambda s: abs(s[0] - 60 * mins))
        xy, tri, f = load_step(base, step, ("liquid_fraction",))
        tr = paper_triangulation(xy, tri)
        ax.tricontourf(tr, f["liquid_fraction"], levels=[-0.1, 0.5, 1.1], colors=["0.1", "white"])
        cs = ax.tricontour(tr, f["liquid_fraction"], levels=[0.5], colors="tab:red", linewidths=1.2)
        segs = [s for lvl in cs.allsegs for s in lvl]
        if segs:
            front = np.vstack(segs)
            front = front[np.argsort(front[:, 1])]
            np.savetxt(os.path.join(out, "front_t%03d.csv" % mins), front, delimiter=",",
                       header="s_mm (from hot wall),z_mm", comments="")
        e = read_exp(exp, "fig5_%03d.csv" % mins)
        if e:
            ax.plot(e[0], e[1], "o", ms=2, color="tab:blue")
        ax.set_title("%d min" % mins, fontsize=9)
        ax.set_xlim(0, 1e3 * W)
        ax.set_ylim(0, 1e3 * H)
        ax.set_aspect("equal")
        ax.set_xlabel("s (mm)")
        if abs(t - 60 * mins) > 1:
            print("warning: %d min not written, using t = %.0f s" % (mins, t))
    axs[0].set_ylabel("z (mm)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "front.png"), dpi=200)
    plt.close(fig)

    with open(os.path.join(run, "NaturalConvectionPCM.csv")) as fcsv:
        rows = list(csv.DictReader(fcsv))
    t = np.array([float(r["t"]) for r in rows]) / 60
    y = np.array([float(r["liquidFraction"]) for r in rows])
    fig, ax = plt.subplots(figsize=(5, 3.4))
    ax.plot(t, y)
    ax.set_xlabel("t (min)")
    ax.set_ylabel("mean liquid fraction (-)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "liquid.png"), dpi=200)
    plt.close(fig)


def steady(run, out, exp):
    base, steps = read_xdmf(os.path.join(run, "NaturalConvectionPCM.xdmf"))
    t, step = steps[-1]
    xy, tri, f = load_step(base, step, ("velocity",))
    tr = paper_triangulation(xy, tri)
    vz = mtri.LinearTriInterpolator(tr, f["velocity"][:, 0])    # u_x: upwards
    vs = mtri.LinearTriInterpolator(tr, -f["velocity"][:, 1])   # -u_y: hot -> cold
    s = np.linspace(0, 1e3 * W, 601)
    z = np.linspace(0, 1e3 * H, 1001)
    fig, axs = plt.subplots(2, 2, figsize=(9, 6.5))
    for k, (ax, zc) in enumerate(zip(axs.flat, (41, 82, 123))):
        v = np.asarray(vz(s, np.full_like(s, zc)).filled(np.nan))
        np.savetxt(os.path.join(out, "vz_z%03d.csv" % zc), np.c_[s, v], delimiter=",",
                   header="s_mm (from hot wall),v_z_m/s", comments="")
        ax.plot(s, v, color="tab:red", label="present")
        e = read_exp(exp, "fig7_%d.csv" % (k + 1))
        if e:
            ax.plot(e[0], e[1], "o-", ms=2, color="tab:blue", label="Biwole et al. (exp.)")
        ax.set_title("(%d) z = %d mm" % (k + 1, zc), fontsize=10)
        ax.set_xlabel("s (mm)")
        ax.set_ylabel(r"$v_z$ (m/s)")
        ax.set_ylim(-0.006, 0.006)
        ax.grid(alpha=0.3)
    v = np.asarray(vs(np.full_like(z, 15.0), z).filled(np.nan))
    np.savetxt(os.path.join(out, "vs_s15.csv"), np.c_[z, v], delimiter=",",
               header="z_mm,v_s_m/s (hot -> cold)", comments="")
    ax = axs[1, 1]
    ax.plot(z, v, color="tab:red", label="present")
    e = read_exp(exp, "fig7_4.csv")
    if e:
        ax.plot(e[0], e[1], "o-", ms=2, color="tab:blue", label="Biwole et al. (exp.)")
    ax.set_title("(4) s = 15 mm", fontsize=10)
    ax.set_xlabel("z (mm)")
    ax.set_ylabel(r"$v_s$ (m/s)")
    ax.set_ylim(-0.001, 0.001)
    ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle("steady state, t = %.0f min" % (t / 60), fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "profiles.png"), dpi=200)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("--exp", default=None)
    a = ap.parse_args()
    for case, fn in (("melt", melt), ("steady", steady)):
        run = os.path.join(a.out, case)
        if os.path.exists(os.path.join(run, "NaturalConvectionPCM.xdmf")):
            fn(run, run, a.exp)
            print("%s: written in %s" % (case, run))


if __name__ == "__main__":
    main()
