#!/usr/bin/env python3
"""¿Hay una dosis de factor tisular con la que la FA coagula y el ritmo sinusal no?

La auricula no usa la densidad de FT medida de Kastrup: fija j_X por la
convencion j_X = j_0 k_dX / kappa, pensada para que la iniciacion en estasis no
dependa de kappa (y lo logra: kappa X_ss = j_0 chi E~ / delta). Pero j_0 no esta
medido. Con la constante por defecto, j_X ~ 2e-15 mol m^-2 s^-1, unas 200 veces
menos que un parche trombogenico de Kastrup.

Entonces los tiempos absolutos de formacion de trombo no son una prediccion. Lo
que si lo es, y es lo que mide este script, es el contraste a dosis fija: se
barre j_0 y se busca la ventana en la que la FA enciende dentro de la ventana de
simulacion y el ritmo sinusal no. Si esa ventana existe, es el punto de
operacion del modelo y lo que hay que reportar; si no existe, el modelo no
discrimina y eso tambien es un resultado.

Por defecto corre SIN semilla de trombina, que es la pregunta de iniciacion
espontanea. Con --seed se agrega la semilla y entonces mide retencion.
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import ptagx_common as px  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_Q1 as q1  # noqa: E402

pc = px.pc
HERE = os.path.dirname(os.path.abspath(__file__))
J0 = 2.0e-12                 # j_0 por defecto de LeftAtrium2DPTAGX (mol m^-2 s^-1)
DOSES = [1.0, 3.0, 10.0, 30.0, 100.0, 300.0]
BURST = 1.0e-4               # 100 nM de trombina
KASTRUP = 4.0e-13            # dosis medida de FT, para poner j_X en escala


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--exe-x", default=px.EXE_X)
    ap.add_argument("--mpirun", default=os.environ.get("PTAGP_MPIRUN", "mpirun -n 1"))
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--plot-only", action="store_true")
    ap.add_argument("--shapes", default="CW")
    ap.add_argument("--kinetics", default="ptagx")
    ap.add_argument("--warmup", type=int, default=3)
    ap.add_argument("--cycles", type=int, default=60)
    ap.add_argument("--substeps", type=int, default=4)
    ap.add_argument("--seed", action="store_true", help="agregar la semilla de trombina")
    ap.add_argument("--doses", default=",".join(f"{d:g}" for d in DOSES))
    args = ap.parse_args()

    out = os.path.abspath(args.out)
    sub = os.path.join(out, "dose" + ("_seeded" if args.seed else ""))
    os.makedirs(sub, exist_ok=True)
    cal = json.load(open(os.path.join(out, "calibration.json")))
    doses = [float(x) for x in args.doses.split(",")]
    shapes = args.shapes.split(",")
    kopt = {"qureshi": q1.QURESHI,
            "ptagx": q1.ptagx_options(px.MODEL),
            "ptagx_mono": q1.ptagx_options(px.MODEL, ptag_fibrin_monomer=True),
            "ptagx_notm": q1.ptagx_options(px.MODEL, ptag_k_tm=0.0),
            "ptagx_base": q1.ptagx_options(px.BASE)}[args.kinetics]
    seeds_of = {s: json.load(open(os.path.join(q1.MESH, f"Q_{s}.json")))["seeds"] for s in shapes}

    cases = {}
    for s in shapes:
        for rh in ("SR", "AF"):
            for d in doses:
                o = dict(q1.flow_options(s, rh, cal[rh]["wall"]), **kopt,
                         ptag_wall_flux=d * J0,
                         ptag_flow_cycles=args.warmup, ptag_species_cycles=args.cycles,
                         ptag_replay=True, ptag_species_substeps=args.substeps)
                if args.seed:
                    x, y = seeds_of[s][args.seed if isinstance(args.seed, str) else "P2"]
                    o["ptag_seed"] = (x, y, q1.SEED_R, q1.SEED_T)
                cases[f"{s}_{rh}_d{d:g}"] = o
    print(f"dose scan: {len(cases)} corridas ({'con' if args.seed else 'sin'} semilla)", flush=True)
    csvs = q1.run_many(args, sub, cases)

    def ignite(csv):
        """(t de encendido global, ¿encendio primero en la orejuela?, fibrina final,
        fraccion de la fibrina que esta en la orejuela)."""
        dd = pc.read_csv(csv)
        t = dd["t"] - args.warmup
        m = t >= 0
        t = t[m]
        tg = pc.first_time(t, dd["maxT"][m] >= BURST)
        tl = pc.first_time(t, dd["laaMaxT"][m] >= BURST)
        intF, laaF = float(dd["intF"][m][-1]), float(dd["laaF"][m][-1])
        frac = laaF / intF if intF > 0 else np.nan
        return tg, bool(np.isfinite(tl) and (not np.isfinite(tg) or tl <= tg + 1e-9)), intF, frac

    lines = [f"Cinetica: `{args.kinetics}`. Semilla: {'si' if args.seed else 'no'}. "
             f"Ventana: {args.cycles} s. j_0 por defecto = {J0:.1e} mol m^-2 s^-1, "
             f"o sea j_X = {J0*px.K['k_dX']/px.MODEL['kappa']:.1e} "
             f"({J0*px.K['k_dX']/px.MODEL['kappa']/KASTRUP:.3f} veces la dosis medida de Kastrup).",
             "",
             "| forma | j_0/j_0def | j_X/Kastrup | RS: enciende (s) | FA: enciende (s) | "
             "FA primero en la orejuela | fibrina FA (mol) | fraccion en la orejuela |",
             "|---|---|---|---|---|---|---|---|"]
    window = {}
    for s in shapes:
        for d in doses:
            try:
                tsr = ignite(csvs[f"{s}_SR_d{d:g}"])
                taf = ignite(csvs[f"{s}_AF_d{d:g}"])
            except KeyError:
                continue
            jx = d * J0 * px.K["k_dX"] / px.MODEL["kappa"]
            good = (not np.isfinite(tsr[0])) and np.isfinite(taf[0])
            window.setdefault(s, []).append((d, good))
            lines.append(
                f"| {s} | {d:g} | {jx/KASTRUP:.3f} | "
                f"{'no' if not np.isfinite(tsr[0]) else f'{tsr[0]:.0f}'} | "
                f"{'no' if not np.isfinite(taf[0]) else f'{taf[0]:.0f}'} | "
                f"{'si' if taf[1] else 'no'} | {taf[2]:.2e} | {taf[3]:.2f} |"
                + ("  **ventana**" if good else ""))

    lines.append("")
    for s, ws in window.items():
        good = [d for d, ok in ws if ok]
        if good:
            lines.append(f"- {s}: la FA enciende y el RS no para j_0 entre {min(good):g} y "
                         f"{max(good):g} veces el valor por defecto.")
        else:
            lines.append(f"- {s}: **no hay ventana** en el barrido. O las dos encienden, o "
                         "ninguna. A dosis fija el modelo no separa FA de RS en esta forma.")
    lines += ["",
              "Como leerlo. Si hay ventana, ese es el punto de operacion del modelo y el "
              "resultado reportable es el contraste, no el tiempo. Si la ventana es muy "
              "ancha, el contraste es robusto; si es de un solo punto, el modelo esta "
              "viviendo en un filo y hay que decirlo.",
              "",
              "La columna 'primero en la orejuela' es la que dice si el modelo enciende "
              "donde deberia. Si enciende antes en el cuerpo auricular, es un falso positivo "
              "y hay que mirar la trombomodulina y E."]
    pc.write_summary(sub, "Q1 - barrido de la dosis de factor tisular", lines)
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
