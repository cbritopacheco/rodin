#!/usr/bin/env python3
"""Cross-variant comparison of Q1: which form of the model does best in the atrium.

run_Q1.py --stage report writes one verdict per kinetics. This puts all of them
side by side and ranks them, which is the thing you actually want when the
question is "which model do I keep".

Three groups of criteria, in decreasing order of how much they are worth:

  1. Checks against Qureshi et al. 2021, Table 1 (the six of run_Q1.py). They
     are a comparison against another model, not against an experiment, so a
     perfect score is not validation; it means the transport and the washout
     behave like the reference does.
  2. Discrimination. The reason the model exists is to separate a left atrial
     appendage that retains a thrombus from one that washes it out, and atrial
     fibrillation from sinus rhythm. Measured as the ratio of residual fibrin
     AF/SR, averaged over the four shapes. Larger is better; 1 is useless.
  3. Quiescence of the healthy wall. Fibrin that appears outside the appendage
     with no lesion there is a false positive. Measured as the fibrin in the
     appendage over the fibrin in the whole domain at 60 s.

Nothing here is a substitute for an experiment. It ranks the variants against
each other and against the only published 2D reference; read it as such.
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
SHAPES = ("CW", "WS", "BR", "CA")

WHAT = {
    "qureshi": "cinética del artículo (trazador pasivo, sin realimentación)",
    "ptagx": "PTAGX-c: cofactor Va + trombomodulina, sin monómero",
    "ptagx_mono": "PTAGX-c con el monómero de fibrina advectado",
    "ptagx_notm": "PTAGX-c sin trombomodulina",
    "ptagx_base": "PTAGX antes de los cambios (kappa 3, eps 1)",
    "ptagx_surface": "PTAGX con amplificacion superficial sigma",
    "ptagx_legacy": "PTAGX por defecto (k_a 0.10)",
}


def collect(out, warmup):
    runs = [d for d in os.listdir(out) if os.path.isdir(os.path.join(out, d))
            and not d.startswith("calib")
            and os.path.exists(os.path.join(out, d, "done"))]
    M = {}
    for r in runs:
        try:
            M[r] = q1.metrics(os.path.join(out, r, "run.csv"), warmup * 1.0)
        except Exception as e:                       # a run that died mid-way
            print(f"  (se ignora {r}: {e})", flush=True)
    return M


def score_table1(out, kin):
    """How many of the six Qureshi checks that run_Q1 wrote this variant passes."""
    p = os.path.join(out, kin, "verdict.md")
    if not os.path.exists(p):
        return None, None
    lines = [l for l in open(p, encoding="utf8") if l.startswith("- [")]
    ok = sum(1 for l in lines if l.startswith("- [PASA]"))
    return ok, len(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=os.path.join(HERE, "results"))
    ap.add_argument("--warmup", type=int, default=3)
    ap.add_argument("--seed", default="P2")
    args = ap.parse_args()
    out = os.path.abspath(args.out)

    M = collect(out, args.warmup)
    kins = sorted({r.rsplit("_", 3)[0] for r in M})
    if not kins:
        sys.exit(f"no hay corridas terminadas en {out}")

    rows, detail = [], []
    for kin in kins:
        def g(s, rh, key):
            return M.get(f"{kin}_{s}_{rh}_{args.seed}", {}).get(key, np.nan)

        res_sr = np.array([g(s, "SR", "Fn_res") for s in SHAPES])
        res_af = np.array([g(s, "AF", "Fn_res") for s in SHAPES])
        tf_af = np.array([g(s, "AF", "t_TF") for s in SHAPES])
        # discrimination AF/SR: fibrin left at 60 s, floor to keep it finite
        ratio = np.nanmean(res_af / np.maximum(res_sr, 1e-3))
        # confinement to the appendage
        laa = np.array([g(s, rh, "laaF") for s in SHAPES for rh in ("SR", "AF")])
        peak = np.array([g(s, rh, "Fn_peak") for s in SHAPES for rh in ("SR", "AF")])
        ok, tot = score_table1(out, kin)
        rows.append(dict(kin=kin, ok=ok, tot=tot, ratio=ratio,
                         res_sr=np.nanmean(res_sr), res_af=np.nanmean(res_af),
                         tf_af=np.nanmean(tf_af), peak=np.nanmean(peak),
                         laaF=np.nanmean(laa)))
        detail.append((kin, res_sr, res_af, tf_af))

    rows.sort(key=lambda r: (-(r["ok"] or 0), -r["ratio"]))

    lines = ["Cada fila es una forma del modelo. El orden es por checks de Qureshi y,",
             "a igualdad, por discriminacion FA/RS.", "",
             "| variante | checks Qureshi | Fn_res RS (uM) | Fn_res FA (uM) | FA/RS | t_TF FA (s) | que es |",
             "|---|---|---|---|---|---|---|"]
    for r in rows:
        sc = "-" if r["ok"] is None else f"{r['ok']}/{r['tot']}"
        lines.append(f"| `{r['kin']}` | {sc} | {r['res_sr']:.2f} | {r['res_af']:.2f} | "
                     f"{r['ratio']:.1f} | {r['tf_af']:.1f} | {WHAT.get(r['kin'], '')} |")

    lines += ["", "Fibrina residual a los 60 s por forma de orejuela (semilla "
              f"{args.seed}), en uM:", "",
              "| variante | " + " | ".join(f"{s} RS" for s in SHAPES) + " | "
              + " | ".join(f"{s} FA" for s in SHAPES) + " |",
              "|---" * (1 + 2 * len(SHAPES)) + "|"]
    for kin, rs, af, _ in detail:
        lines.append(f"| `{kin}` | " + " | ".join(f"{v:.2f}" for v in rs) + " | "
                     + " | ".join(f"{v:.2f}" for v in af) + " |")

    best = rows[0]
    lines += ["", "## Lectura", "",
              f"Mejor por este criterio: **{best['kin']}** "
              f"({WHAT.get(best['kin'], '')}), con {('-' if best['ok'] is None else str(best['ok']))} "
              f"checks de Qureshi y discriminacion FA/RS de {best['ratio']:.1f}.", "",
              "Advertencias que hay que leer antes de creerle a la tabla:", "",
              "- Los checks son contra otro modelo, no contra un experimento. Un 6/6",
              "  significa que el transporte y el lavado se comportan como la referencia.",
              "- La dosis de factor tisular de la auricula (j_0) no esta medida. Los tiempos",
              "  absolutos no son predicciones; el contraste a dosis fija si lo es.",
              "  Ver dose_scan.py.",
              "- Las corridas llevan semilla de trombina, asi que esto mide propagacion y",
              "  lavado, no iniciacion espontanea. Para iniciacion, dose_scan.py --no-seed.",
              ]
    pc.write_summary(out, "Q1 - comparacion de variantes del modelo", lines)
    print("\n".join(lines), flush=True)


if __name__ == "__main__":
    main()
