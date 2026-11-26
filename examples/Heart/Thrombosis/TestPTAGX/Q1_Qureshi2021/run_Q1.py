#!/usr/bin/env python3
"""Q1 - 2D left atrium of Qureshi et al. (CinC 2021), four appendage shapes.

Benchmark (numerical, not experimental): Qureshi A, Balmus M, Nechipurenko D,
Ataullakhanov F, Williams S, Lip G, Nordsletten D, Aslanidi O, de Vecchi A.
"Left atrial appendage morphology impacts thrombus formation risks in
multi-physics atrial models", Computing in Cardiology 48 (2021).
doi:10.22489/CinC.2021.241. Meshes: github.com/ahmedqureshi1/2D_LAA.

Their set-up, reproduced here:
  - 2D LA with four PV inlets, MV outlet, LAA; shapes CW, WS, BR, CA;
  - pulsatile PV inflow, peak 7.5 cm/s (SR) and 5.0 cm/s (AF);
  - Newtonian blood, rho 1060 kg/m^3, mu 3.5e-3 Pa s; dt 0.01 s; 60 s;
  - LAA wall moving along its normal, sinusoidal, strong in SR, weak in AF;
  - thrombin seed T = 0.09 mmol/m^3 (90 nM) at P1, P2 or P3 inside the LAA;
  - fibrinogen 4.5 mmol/m^3 (uM); Fg -> Fn at K_eff Fg Th, K_eff = 7180 m^3/(mol s);
    D_Fg = D_Fn = 2e-11, D_Th = 4.6e-11 m^2/s.
  - thrombus formation (TF): Fn >= 3.6 uM (80 % of Fg) in the LAA;
    washout: Fn < 0.5 uM in the LAA; Fn_res = Fn in the LAA at 60 s.

Their Table 1 (seed at P2):
         SR: V_mean (cm/s)  Fn_res  t_washout (s) | AF: V_mean  Fn_res  t_TF (s)
  CW         0.12            0.1     28           |    0.022     3.7     17
  WS         0.11            0.9     no washout   |    0.019     3.6     16
  BR         0.11            1.3     no washout   |    0.011     3.8     14
  CA         0.13            0.5     42           |    0.031     3.7     20
and: seed at P1 washes out in < 5 s (all shapes, SR and AF); seed at P3 gives TF
in all cases within 60 s.

Their kinetics as a limit of PTAGX: Fn = G_0 - G with D_G = D_Fn is exactly
their Fn equation; R_F = k_F T G/(K_G + G) with K_G >> G and k_F = K_eff K_G
is their K_eff Fg Th; k_a = k_on = 0 and no wall source leave T as a passive
tracer. Their R_Th (Ataullakhanov's model) is not given in the paper, so this
limit has R_Th = 0 ("passive seed"); it isolates the flow (washout vs stasis).
The other kinetics are PTAGX itself (amplification, antithrombin, endocardial
Xa, immobile fibrin, Brinkman), with the same seed: "ptagx" = the proposed
parameters (ptagx_common.MODEL: k_a of E2), "ptagx_surface" = MODEL plus the
surface amplification sigma of E5 on the endocardial layer (sensitivity),
"ptagx_legacy" = the PTAGX defaults.

Criteria (each kinetics, seed P2), from the paper: BR keeps the most fibrin in
SR and CW the least; CW and CA wash out in SR, WS and BR do not; in AF every
shape forms thrombus within 60 s, BR first and CA last; V_mean SR/AF between 4
and 10; in AF BR has the slowest LAA and CA the fastest.

Not given in the paper, and set here (see README):
  - the PV waveform shape (S, D and, in SR, A-reversal waves; period 1 s);
  - the LAA wall-velocity amplitude: calibrated on CW so that the mean LAA
    speed matches Table 1 (0.12 cm/s SR, 0.022 cm/s AF); the other three shapes
    then predict their V_mean;
  - the moving wall is a transpiration velocity on a fixed mesh (no ALE).

Stages (each skips the runs whose CSV exists):
  --stage calib   flow only, CW, wall amplitude sweep (6 runs x 4 cycles)
  --stage coag    shapes x {SR, AF} x seeds x kinetics, 60 s
  --stage report  tables and figures only
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import ptagx_common as px  # noqa: E402

pc = px.pc
MESH = os.path.join(HERE, "meshes")
WAVE = os.path.join(HERE, "waveforms")
TABLE1 = {  # shape: (SR Vmean cm/s, SR Fn_res, SR t_washout or None, AF Vmean, AF Fn_res, AF t_TF)
    "CW": (0.12, 0.1, 28.0, 0.022, 3.7, 17.0),
    "WS": (0.11, 0.9, None, 0.019, 3.6, 16.0),
    "BR": (0.11, 1.3, None, 0.011, 3.8, 14.0),
    "CA": (0.13, 0.5, 42.0, 0.031, 3.7, 20.0),
}
PEAK = {"SR": 0.075, "AF": 0.050}
# LAA contraction peaks at t = 0.9 s (atrial systole) in both rhythms
PHASE = np.pi / 2 - 2 * np.pi * 0.9 + 2 * np.pi
SEED_R, SEED_T = 2.0e-3, 9.0e-5
FN_TF, FN_WASH = 3.6e-3, 0.5e-3

QURESHI = dict(ptag_xa=False, ptag_immobile_fibrin=False, ptag_k_a=0.0, ptag_k_on=0.0,
               ptag_wall_flux=0.0, ptag_K_G=1.0, ptag_k_F=7180.0, ptag_fibrinogen=4.5e-3,
               ptag_D_G=2.0e-11, ptag_D_P=4.6e-11, ptag_brinkman=False)
PTAGX = dict(ptag_xa=True, ptag_immobile_fibrin=True)


def ptagx_options(model, **over):
    """LeftAtrium2DPTAGX options of a parameter set of ptagx_common: k_a, sigma,
    the Xa channel (kappa, epsilon), the monomer and the thrombomodulin."""
    o = dict(PTAGX, ptag_k_a=model["k_a"], ptag_sigma_s=model["sigma_s"],
             ptag_kappa_x=model.get("kappa", 3.0), ptag_eps_x=model.get("eps_x", 1.0),
             ptag_fibrin_monomer=model.get("monomer", False),
             ptag_k_tm=model.get("k_tm", 0.0))
    o.update(over)
    return o


def flow_options(shape, rhythm, wall):
    m = os.path.join(MESH, f"Q_{shape}.mesh")
    return dict(ptag_mesh=m, ptag_pv=os.path.join(WAVE, f"pv_velocity_{rhythm}.dat"),
                ptag_mv=os.path.join(WAVE, "mv_pressure_zero.dat"),
                ptag_inlet_velocity=PEAK[rhythm], ptag_mu=3.5e-3, ptag_dt=0.01,
                ptag_laa_wall_velocity=wall, ptag_laa_wall_phase=PHASE,
                ptag_laa_poly=os.path.join(MESH, f"Q_{shape}_laa.txt"),
                ptag_output_every=0, ptag_skew=True, ptag_csv_every=5)


def run(args, outdir, name, opts):
    """Run LeftAtrium2DPTAGX in its own folder (it writes <xdmf>.* there)."""
    import shlex
    import subprocess
    d = os.path.join(outdir, name)
    os.makedirs(d, exist_ok=True)
    csv = os.path.join(d, "run.csv")
    if os.path.exists(csv) and not args.force and os.path.exists(os.path.join(d, "done")):
        return csv
    if args.plot_only:
        raise FileNotFoundError(csv)
    opts = dict(opts, ptag_csv=csv, ptag_xdmf=os.path.join(d, "run"))
    cmd = shlex.split(args.mpirun) + [args.exe_x]
    for k, v in opts.items():
        if isinstance(v, (list, tuple)):
            v = ",".join(f"{x:.9g}" for x in v)
        elif isinstance(v, bool):
            v = str(v).lower()
        cmd += ["-" + k, str(v)]
    with open(os.path.join(d, "log.txt"), "w") as f:
        f.write(" ".join(cmd) + "\n")
        f.flush()
        r = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, cwd=d)
    if r.returncode != 0:
        raise RuntimeError(f"{name} failed, see {d}/log.txt")
    open(os.path.join(d, "done"), "w").close()
    return csv


def run_many(args, outdir, cases):
    import concurrent.futures as cf
    out = {}
    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        futs = {ex.submit(run, args, outdir, n, o): n for n, o in cases.items()}
        for fu in cf.as_completed(futs):
            n = futs[fu]
            try:
                out[n] = fu.result()
                print(f"  done: {n}", flush=True)
            except Exception as e:  # keep the others going
                print(f"  FAILED: {n}: {e}", flush=True)
    return out


def last_cycle_mean(csv, col="laaMeanU", period=1.0):
    d = pc.read_csv(csv)
    t = d["t"]
    m = t > t[-1] - period
    return float(np.mean(d[col][m]))


def calib(args, out):
    """Wall amplitude W per rhythm so that the CW mean LAA speed matches Table 1.

    W = 0 is run first; wall motion only adds flow, so if the still wall already
    exceeds the target, W = 0 is kept and the sweep is skipped.
    """
    sweep = {"SR": [2e-3, 5e-3], "AF": [0.5e-3, 1.5e-3]}

    def case(rh, w):
        return f"calib_CW_{rh}_W{w*1e3:g}mm", dict(
            flow_options("CW", rh, w), ptag_flow_cycles=args.calib_cycles,
            ptag_species_cycles=0, ptag_replay=False)

    first = dict(case(rh, 0.0) for rh in sweep)
    csvs = run_many(args, out, first)
    res = {}
    for rh in sweep:
        target = TABLE1["CW"][0] if rh == "SR" else TABLE1["CW"][3]
        ws = [0.0]
        v = [last_cycle_mean(csvs[f"calib_CW_{rh}_W0mm"]) * 100]
        if v[0] < target:
            more = dict(case(rh, w) for w in sweep[rh])
            c2 = run_many(args, out, more)
            ws += sweep[rh]
            v += [last_cycle_mean(c2[f"calib_CW_{rh}_W{w*1e3:g}mm"]) * 100 for w in sweep[rh]]
            w_star = float(np.interp(target, v, ws)) if v[-1] >= target else float(ws[-1] * target / v[-1])
        else:
            w_star = 0.0
        res[rh] = dict(walls=ws, vmean_cm_s=v, target=target, wall=w_star,
                       note="still wall already above the target" if v[0] >= target else "")
        print(f"  {rh}: Vmean(W) = " + ", ".join(f"{w*1e3:g} mm/s -> {x:.3f} cm/s" for w, x in zip(ws, v))
              + f"; W* = {w_star*1e3:.2f} mm/s for {target} cm/s", flush=True)
    json.dump(res, open(os.path.join(out, "calibration.json"), "w"), indent=1)
    return res


def coag(args, out):
    cal = json.load(open(os.path.join(out, "calibration.json")))
    shapes = args.shapes.split(",")
    seeds_of = {s: json.load(open(os.path.join(MESH, f"Q_{s}.json")))["seeds"] for s in shapes}
    cases = {}
    for kin in args.kinetics.split(","):
        kopt = {"qureshi": QURESHI,
                "ptagx": ptagx_options(px.MODEL),                  # PTAGX-Va, the proposal
                "ptagx_base": ptagx_options(px.BASE),              # PTAGX before the changes
                "ptagx_notm": ptagx_options(px.MODEL, ptag_k_tm=0.0),
                "ptagx_mono": ptagx_options(px.MODEL, ptag_fibrin_monomer=True),
                "ptagx_surface": ptagx_options(px.SURFACE),
                "ptagx_legacy": ptagx_options(px.LEGACY)}[kin]
        for s in shapes:
            for rh in ("SR", "AF"):
                for p in args.seeds.split(","):
                    x, y = seeds_of[s][p]
                    cases[f"{kin}_{s}_{rh}_{p}"] = dict(
                        flow_options(s, rh, cal[rh]["wall"]), **kopt,
                        ptag_seed=(x, y, SEED_R, SEED_T),
                        ptag_flow_cycles=args.warmup, ptag_species_cycles=args.cycles,
                        ptag_replay=True, ptag_species_substeps=args.substeps)
    print(f"Q1 coag: {len(cases)} runs")
    return run_many(args, out, cases)


def metrics(csv, warmup):
    d = pc.read_csv(csv)
    t = d["t"] - warmup                       # time since the seed
    m = t >= 0
    t, F, V = t[m], d["laaMaxF"][m], d["laaMeanU"][m]
    tf = pc.first_time(t, F >= FN_TF)
    peak = F.max()
    # washout: below the threshold after having exceeded it
    ipk = int(np.argmax(F))
    after = np.flatnonzero(F[ipk:] < FN_WASH)
    tw = t[ipk + after[0]] if after.size else np.nan
    return dict(Vmean=100 * float(np.mean(V)), Fn_res=1e3 * float(F[-1]), Fn_peak=1e3 * float(peak),
                t_TF=tf, t_wash=tw, laaF=float(d["laaF"][m][-1]))


def report(args, out):
    plt = pc.setup_matplotlib()
    lines = []
    cal_path = os.path.join(out, "calibration.json")
    if os.path.exists(cal_path):
        cal = json.load(open(cal_path))
        lines += ["## Calibration of the LAA wall velocity (CW)", "",
                  "| rhythm | W (mm/s) | Vmean LAA (cm/s) |", "|---|---|---|"]
        for rh, c in cal.items():
            for w, v in zip(c["walls"], c["vmean_cm_s"]):
                lines.append(f"| {rh} | {w*1e3:g} | {v:.3f} |")
            lines.append(f"| {rh} | **{c['wall']*1e3:.2f}** (chosen) | target {c['target']} |")
        lines.append("")
    runs = sorted(d for d in os.listdir(out) if os.path.isdir(os.path.join(out, d))
                  and not d.startswith("calib") and os.path.exists(os.path.join(out, d, "done")))
    if not runs:
        pc.write_summary(out, "Q1 - Qureshi 2021", lines or ["(no runs)"])
        return
    M = {r: metrics(os.path.join(out, r, "run.csv"), args.warmup * 1.0) for r in runs}
    for kin in sorted({r.rsplit("_", 3)[0] for r in runs}):
        lines += [f"## Kinetics: {kin}", "",
                  "| shape | rhythm | seed | Vmean (cm/s) [Q] | Fn peak (uM) | Fn_res 60 s (uM) [Q] | "
                  "t_TF (s) [Q] | t_washout (s) [Q] |", "|---|---|---|---|---|---|---|---|"]
        for s in ("CW", "WS", "BR", "CA"):
            for rh in ("SR", "AF"):
                for p in ("P1", "P2", "P3"):
                    r = f"{kin}_{s}_{rh}_{p}"
                    if r not in M:
                        continue
                    m = M[r]
                    q = TABLE1[s]
                    qv = q[0] if rh == "SR" else q[3]
                    qf = q[1] if rh == "SR" else q[4]
                    qtf = "-" if rh == "SR" else f"{q[5]:g}"
                    qtw = ("no" if q[2] is None else f"{q[2]:g}") if rh == "SR" else "-"
                    if p != "P2":
                        qv = qf = qtf = qtw = "-"
                    lines.append(
                        f"| {s} | {rh} | {p} | {m['Vmean']:.3f} [{qv}] | {m['Fn_peak']:.2f} | "
                        f"{m['Fn_res']:.2f} [{qf}] | {m['t_TF']:.1f} [{qtf}] | {m['t_wash']:.1f} [{qtw}] |")
        lines.append("")
        # ordering checks of the paper (seed P2)
        def get(s, rh, key):
            return M.get(f"{kin}_{s}_{rh}_P2", {}).get(key, np.nan)
        chk = []
        res_sr = {s: get(s, "SR", "Fn_res") for s in TABLE1}
        tf_af = {s: get(s, "AF", "t_TF") for s in TABLE1}
        if all(np.isfinite(list(res_sr.values()))):
            chk.append(("BR has the largest SR residual fibrin, CW the smallest",
                        max(res_sr, key=res_sr.get) == "BR" and min(res_sr, key=res_sr.get) == "CW"))
            chk.append(("SR: CW and CA wash out (< 0.5 uM at 60 s), WS and BR do not",
                        res_sr["CW"] < 0.5 and res_sr["CA"] < 0.5 and res_sr["WS"] >= 0.5 and res_sr["BR"] >= 0.5))
        if all(np.isfinite(list(tf_af.values()))):
            chk.append(("AF: TF in all shapes within 60 s", all(v <= 60 for v in tf_af.values())))
            chk.append(("AF: BR clots first, CA last",
                        min(tf_af, key=tf_af.get) == "BR" and max(tf_af, key=tf_af.get) == "CA"))
        v_sr = {s: get(s, "SR", "Vmean") for s in TABLE1}
        v_af = {s: get(s, "AF", "Vmean") for s in TABLE1}
        if all(np.isfinite(list(v_sr.values()) + list(v_af.values()))):
            chk.append(("Vmean SR / AF between 4 and 10 for every shape (Table 1: 4.2-10)",
                        all(4 <= v_sr[s] / max(v_af[s], 1e-9) <= 10 for s in TABLE1)))
            chk.append(("AF: BR slowest LAA, CA fastest", min(v_af, key=v_af.get) == "BR"
                        and max(v_af, key=v_af.get) == "CA"))
        if chk:
            lines += ["Checks against Qureshi (seed P2):", ""]
            lines += [f"- {txt}: {'PASA' if ok else 'FALLA'}" for txt, ok in chk]
            lines.append("")
            os.makedirs(os.path.join(out, kin), exist_ok=True)
            px.verdict(os.path.join(out, kin), f"Q1 Qureshi 2021 ({kin})", [(t, ok, "") for t, ok in chk])
        # figure: LAA max fibrin vs time
        fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), sharey=True)
        for j, rh in enumerate(("SR", "AF")):
            for s in ("CW", "WS", "BR", "CA"):
                r = f"{kin}_{s}_{rh}_P2"
                if r in M:
                    d = pc.read_csv(os.path.join(out, r, "run.csv"))
                    axs[j].plot(d["t"] - args.warmup, d["laaMaxF"] * 1e3, label=s)
            axs[j].axhline(FN_TF * 1e3, color="k", ls=":", lw=0.8)
            axs[j].axhline(FN_WASH * 1e3, color="k", ls="--", lw=0.8)
            axs[j].set(xlabel="t since seed (s)", title=f"{rh}, seed P2, {kin}", xlim=(0, None))
        axs[0].set_ylabel("max fibrin in LAA (uM)")
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(out, f"Q1_fibrin_{kin}.png"))
    lines += ["[Q] = Qureshi et al. 2021, Table 1 (seed P2 only).", "",
              "Reference: Qureshi A, et al. Computing in Cardiology 48 (2021). doi:10.22489/CinC.2021.241"]
    pc.write_summary(out, "Q1 - 2D left atrium of Qureshi et al. (CinC 2021)", lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=["calib", "coag", "report", "all"], default="all")
    ap.add_argument("--exe-x", default=px.EXE_X)
    ap.add_argument("--mpirun", default=os.environ.get("PTAGP_MPIRUN", "mpirun -n 1"))
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--shapes", default="CW,WS,BR,CA")
    ap.add_argument("--seeds", default="P2")
    ap.add_argument("--kinetics", default="qureshi,ptagx",
                    help="comma list of qureshi, ptagx (PTAGX-Va, proposed), ptagx_base (before), "
                         "ptagx_notm, ptagx_mono, ptagx_surface, ptagx_legacy")
    ap.add_argument("--calib-cycles", type=int, default=4)
    ap.add_argument("--warmup", type=int, default=3, help="flow-only cycles before the seed")
    ap.add_argument("--cycles", type=int, default=60, help="cycles (s) after the seed")
    ap.add_argument("--substeps", type=int, default=4)
    a = ap.parse_args()
    a.out = os.path.abspath(a.out)       # each case runs in its own folder
    os.makedirs(a.out, exist_ok=True)
    if a.stage in ("calib", "all"):
        calib(a, a.out)
    if a.stage in ("coag", "all"):
        coag(a, a.out)
    report(a, a.out)


if __name__ == "__main__":
    main()
