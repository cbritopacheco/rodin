#!/usr/bin/env python3
"""
Live/static diagnostics dashboard for the 0D poroelastic sphere.

This script reads the CSV produced by PoroelasticSphere_0D and provides:

1. Static diagnostics:
   - saves individual PNG figures,
   - writes summary tables,
   - writes per-cycle haemodynamics,
   - writes validation flags,
   - writes a correlation matrix.

2. Interactive dashboard:
   - grouped windows for cavity/LV dynamics,
   - wall poromechanics (porosity, multiplier, interstitial pressure),
   - coronary perfusion and the systolic impediment,
   - validation/conservation,
   - internal variables,
   - histograms,
   - cycle overlays,
   - correlation matrix with colorbar.

3. Live dashboard:
   - repeatedly reloads the CSV,
   - tolerates partially written rows,
   - refreshes dashboard windows in-place,
   - keeps the correlation-matrix window in its moved position,
   - creates the colorbar once to avoid matplotlib recursion bugs.

The CSV columns are those written by examples/Heart/PoroelasticSphere_0D.cpp:

    t, y, phi, pv, par, pd, ec, gamma, beta, w, kc, tauc,
    V, Q, pat, lambdaBar, pf, Qcor

Model parameters that are not in the CSV (geometry, storage modulus,
perfusion conductances, venous pressure) are needed for the derived
diagnostics and conservation checks. They default to the values used in
PoroelasticSphere_0D.cpp and can be overridden on the command line; keep
them in sync with the run being analyzed.

Examples
--------
Save all diagnostics:

    python PoroelasticSphere_0D_Plot.py poroelastic_sphere_0d_cycle.csv

Show static dashboard:

    python PoroelasticSphere_0D_Plot.py poroelastic_sphere_0d_cycle.csv --show

Live dashboard:

    python PoroelasticSphere_0D_Plot.py poroelastic_sphere_0d_cycle.csv --watch --interval 1

Live dashboard without saving:

    python PoroelasticSphere_0D_Plot.py poroelastic_sphere_0d_cycle.csv --watch --no-save

Periodically save tables while live-monitoring:

    python PoroelasticSphere_0D_Plot.py poroelastic_sphere_0d_cycle.csv --watch --save-every 20

Analyze a run with different parameters:

    python PoroelasticSphere_0D_Plot.py run.csv --R0 2.36e-2 --d0 1.42e-2 --gamma-ven 1.4e-9
"""

import argparse
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PERIOD = 0.85

# Systole is taken from the activation ramp to the end of the relaxation ramp
# (periodic_activation: tau = 0.13 .. 0.45); the rest of the cycle is diastole.
SYSTOLE = (0.13, 0.45)

PA_PER_MMHG = 133.322368

# numpy renamed trapz to trapezoid in 2.0; support both.
TRAPEZOID = getattr(np, "trapezoid", None) or np.trapz

# Model parameters that the CSV does not carry. Defaults match
# examples/Heart/PoroelasticSphere_0D.cpp; main() overrides them from the CLI.
PARAMS = {
    "R0": 2.4e-2,      # reference cavity radius [m]
    "d0": 1.45e-2,     # reference wall thickness [m]
    "phi0": 0.1,       # reference porosity [-]
    "KPhi": 2.0e5,     # fluid storage modulus [Pa]
    "gammaAr": 7.0e-10,   # arterial perfusion conductance [m^3/(s Pa)]
    "gammaVen": 7.0e-10,  # venous perfusion conductance [m^3/(s Pa)]
    "pSv": 1.0e3,      # venous pressure [Pa]
}


def wall_reference_volume():
    """Reference wall volume V_w0 = 4/3 pi (R_out^3 - R_in^3) [m^3]."""
    rin = PARAMS["R0"]
    rout = PARAMS["R0"] + PARAMS["d0"]
    return 4.0 / 3.0 * np.pi * (rout ** 3 - rin ** 3)


def cols_present(df, cols):
    """Return only columns present in the dataframe."""
    return [c for c in cols if c in df.columns]


def numeric_series(df, col):
    """Return one dataframe column converted robustly to numeric values."""
    return pd.to_numeric(df[col], errors="coerce")


def safe_series(df, col):
    """Return a finite numeric series for diagnostics."""
    return (
        numeric_series(df, col)
        .replace([np.inf, -np.inf], np.nan)
        .dropna()
    )


def read_csv_robust(path):
    """
    Read a CSV that may currently be written by another process.

    In live mode, pandas can catch the file while a row is incomplete, for
    example with a value such as "1e-". Loading as strings and coercing later
    avoids crashing the dashboard.
    """
    try:
        df = pd.read_csv(path, dtype=str, on_bad_lines="skip")
        if df.empty:
            return None
        return add_derived_columns(df)
    except Exception as e:
        print(f"[live] Could not read CSV yet: {e}")
        return None


def bdf2_rate(t, y):
    """
    Discrete rate on the CSV grid using the scheme of the solver.

    The stepper uses backward Euler for the first step and BDF2 afterwards,
    so reproducing that same formula here makes the conservation residuals
    below measure the solver rather than the post-processing. The first two
    samples cannot be formed (the CSV does not contain the initial state) and
    are returned as NaN.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    out = np.full(y.shape, np.nan)

    if y.size < 3:
        return out

    dt = t[1:] - t[:-1]
    with np.errstate(invalid="ignore", divide="ignore"):
        out[2:] = (1.5 * y[2:] - 2.0 * y[1:-1] + 0.5 * y[:-2]) / dt[1:]
    return out


def add_derived_columns(df):
    """
    Add derived physical and numerical diagnostics.

    Unit conversions:
    - volume: m^3 -> mL
    - cavity flow: m^3/s -> mL/s
    - perfusion flow: m^3/s -> mL/min
    - pressure: Pa -> mmHg

    Added diagnostics:
    - volume ratio J and wall/blood volumes,
    - venous outflow from the constitutive law and from the mass balance,
    - mass-conservation residual (absolute and relative),
    - cavity-flow consistency check,
    - storage-law check (p_tilde + lambda_bar against K_Phi (Phi - phi0)),
    - perfusion pressure,
    - cardiac cycle index,
    - cycle-local time.
    """
    df = df.copy()

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    vw0 = wall_reference_volume()
    df["WallReferenceVolume_mL"] = 1e6 * vw0

    # ---- Cavity ---------------------------------------------------------
    if "V" in df:
        df["CavityVolume_mL"] = 1e6 * df["V"]

    if "Q" in df:
        df["CavityFlow_mL_s"] = 1e6 * df["Q"]

    if "y" in df:
        df["EndocardialDisplacement_mm"] = 1e3 * df["y"]
        df["EndocardialRadius_mm"] = 1e3 * (PARAMS["R0"] + df["y"])

    for src, dst in (
        ("pv", "VentricularPressure_mmHg"),
        ("par", "ArterialPressure_mmHg"),
        ("pd", "DistalPressure_mmHg"),
        ("pat", "AtrialPressure_mmHg"),
        ("pf", "InterstitialPressure_mmHg"),
    ):
        if src in df:
            df[dst] = df[src] / PA_PER_MMHG

    if "lambdaBar" in df:
        df["MultiplierBar_kPa"] = 1e-3 * df["lambdaBar"]

    # ---- Wall poromechanics --------------------------------------------
    if "phi" in df:
        df["VolumeRatio_J"] = 1.0 - PARAMS["phi0"] + df["phi"]
        df["WallBloodVolume_mL"] = 1e6 * vw0 * df["phi"]
        df["WallVolume_mL"] = 1e6 * vw0 * df["VolumeRatio_J"]
        df["PorosityChange"] = df["phi"] - PARAMS["phi0"]

    # The storage law p_tilde = K_Phi (Phi - phi0) - lambda_bar must hold
    # exactly; this residual detects a parameter mismatch between the run and
    # the values assumed by this script.
    if all(c in df for c in ("pf", "phi", "lambdaBar")):
        df["StorageLawResidual"] = (
            df["pf"]
            - (PARAMS["KPhi"] * (df["phi"] - PARAMS["phi0"]) - df["lambdaBar"])
        )

    # ---- Perfusion ------------------------------------------------------
    # Arterial inflow is reported by the solver; the venous outflow follows
    # either from its own law or from the fluid mass balance
    # V_w0 dPhi/dt = Q_in - Q_out. Comparing the two is the conservation
    # check of the porosity equation.
    if "Qcor" in df:
        df["CoronaryArterialInflow_mL_min"] = 6e7 * df["Qcor"]

    if "pf" in df:
        df["CoronaryVenousOutflow"] = PARAMS["gammaVen"] * (df["pf"] - PARAMS["pSv"])
        df["CoronaryVenousOutflow_mL_min"] = 6e7 * df["CoronaryVenousOutflow"]

    if all(c in df for c in ("t", "phi")):
        df["PorosityRate"] = bdf2_rate(df["t"].to_numpy(), df["phi"].to_numpy())
        df["WallStorageFlux"] = vw0 * df["PorosityRate"]
        df["WallStorageFlux_mL_min"] = 6e7 * df["WallStorageFlux"]

        if "Qcor" in df:
            df["CoronaryVenousOutflowFromBalance"] = df["Qcor"] - df["WallStorageFlux"]
            df["CoronaryVenousOutflowFromBalance_mL_min"] = (
                6e7 * df["CoronaryVenousOutflowFromBalance"]
            )

    if all(c in df for c in ("CoronaryVenousOutflowFromBalance", "CoronaryVenousOutflow")):
        df["MassBalanceResidual"] = (
            df["CoronaryVenousOutflowFromBalance"] - df["CoronaryVenousOutflow"]
        )
        df["AbsMassBalanceResidual"] = df["MassBalanceResidual"].abs()
        denom = df["Qcor"].abs().replace(0.0, np.nan)
        df["RelativeMassBalanceResidual"] = df["AbsMassBalanceResidual"] / denom

    if all(c in df for c in ("par", "pf")):
        df["PerfusionPressure"] = df["par"] - df["pf"]
        df["PerfusionPressure_mmHg"] = df["PerfusionPressure"] / PA_PER_MMHG

    # Interstitial pressure relative to the cavity: when p_tilde exceeds p_v
    # the wall squeezes its own vasculature harder than the cavity does.
    if all(c in df for c in ("pf", "pv")):
        df["InterstitialMinusCavity"] = df["pf"] - df["pv"]

    # ---- Cavity flow consistency ---------------------------------------
    # The CSV Q column is a first-order difference of V written by the
    # example, so it differs from the solver's own BDF2 rate of V by exactly
    # 0.5 dt Vddot. The mismatch is therefore a measure of the volume
    # acceleration, largest at the activation switch, and not a solver error:
    # for the shipped parameters its median is ~0.1% of the peak flow and it
    # spikes only over the few samples of the upstroke.
    if all(c in df for c in ("t", "V", "Q")):
        df["CavityFlowFromVolume"] = bdf2_rate(df["t"].to_numpy(), df["V"].to_numpy())
        df["CavityFlowMismatch"] = df["Q"] - df["CavityFlowFromVolume"]

    if "t" in df:
        df["cycle"] = np.floor(df["t"] / PERIOD).astype("Int64")
        df["tau"] = df["t"] - PERIOD * df["cycle"].astype(float)

    return df


def phase_split(df, col):
    """Mean of a column over systole and over diastole, plus peak/mean ratio."""
    if col not in df or "tau" not in df:
        return None
    s = numeric_series(df, col)
    tau = numeric_series(df, "tau")
    sys_mask = (tau >= SYSTOLE[0]) & (tau < SYSTOLE[1])
    mean_all = s.mean()
    out = {
        "systolic_mean": s[sys_mask].mean(),
        "diastolic_mean": s[~sys_mask].mean(),
        "mean": mean_all,
        "peak": s.max(),
    }
    out["peak_over_mean"] = (
        out["peak"] / mean_all if mean_all not in (0.0, np.nan) else np.nan
    )
    return out


def last_full_cycle(df):
    """Return the rows of the last complete cardiac cycle, or None."""
    if "cycle" not in df or df["cycle"].dropna().empty:
        return None
    cycles = df["cycle"].dropna().unique()
    if len(cycles) < 2:
        return None
    last = sorted(cycles)[-1]
    grp = df[df["cycle"] == last]
    # The final cycle may be truncated; fall back to the previous one.
    if "tau" in grp and grp["tau"].max() < 0.95 * PERIOD and len(cycles) >= 3:
        last = sorted(cycles)[-2]
        grp = df[df["cycle"] == last]
    return grp


def shade_systole(ax, df):
    """Shade the systolic part of every cycle."""
    if "t" not in df:
        return
    t = numeric_series(df, "t")
    if t.empty:
        return
    tmin, tmax = float(t.min()), float(t.max())
    n = int(np.floor(tmax / PERIOD)) + 1
    for k in range(n):
        lo = k * PERIOD + SYSTOLE[0]
        hi = k * PERIOD + SYSTOLE[1]
        # Clip to the data range so no band is drawn past the last sample.
        lo, hi = max(lo, tmin), min(hi, tmax)
        if hi > lo:
            ax.axvspan(lo, hi, color="0.85", zorder=0, lw=0)


def phased_flows(ax, df):
    """Coronary arterial inflow against venous outflow, with systole shaded.

    Coronary arterial inflow is diastolic dominant and venous outflow is
    systolic dominant. In this model that phase opposition is not imposed: it
    emerges from the reduction, because contraction drives the averaged
    multiplier negative, which raises the interstitial pressure and therefore
    both impedes the inflow and drives the outflow. A physiological result
    shows the two curves peaking in opposite halves of the shaded bands.
    """
    cols = cols_present(
        df,
        [
            "CoronaryArterialInflow_mL_min",
            "CoronaryVenousOutflow_mL_min",
        ],
    )
    if not cols or "t" not in df:
        ax.set_axis_off()
        return
    shade_systole(ax, df)
    ts(ax, df, cols, "Coronary inflow vs venous outflow", "mL/min")

    parts = []
    for c, label in (
        ("CoronaryArterialInflow_mL_min", "art"),
        ("CoronaryVenousOutflow_mL_min", "ven"),
    ):
        st = phase_split(df, c)
        if st is None or not np.isfinite(st["mean"]):
            continue
        parts.append(
            f"{label}: sys {st['systolic_mean']:.0f}, "
            f"dia {st['diastolic_mean']:.0f}, pk/mean {st['peak_over_mean']:.1f}"
        )
    if parts:
        ax.set_title(
            "Coronary inflow vs venous outflow [mL/min]\n" + " | ".join(parts),
            fontsize=8,
        )


def impediment_panel(ax, df):
    """Averaged multiplier and interstitial pressure, with systole shaded.

    The systolic impediment of coronary flow is the chain
    contraction -> lambda_bar < 0 -> p_tilde up -> inflow down, so these two
    curves are the mechanism itself rather than a consequence of it.
    """
    cols = cols_present(df, ["MultiplierBar_kPa"])
    if not cols or "t" not in df:
        ax.set_axis_off()
        return
    shade_systole(ax, df)
    t = numeric_series(df, "t")
    y = numeric_series(df, "MultiplierBar_kPa")
    valid = t.notna() & y.notna()
    ax.plot(t[valid], y[valid], color="C3", label="lambda_bar [kPa]")
    ax.axhline(0.0, color="0.4", lw=0.8)
    ax.set_xlabel("t [s]")
    ax.set_ylabel("lambda_bar [kPa]")
    ax.grid(True, alpha=0.3)

    handles, labels = ax.get_legend_handles_labels()

    if "InterstitialPressure_mmHg" in df:
        ax2 = ax.twinx()
        yp = numeric_series(df, "InterstitialPressure_mmHg")
        valid = t.notna() & yp.notna()
        ax2.plot(t[valid], yp[valid], color="C0", label="p_tilde [mmHg]")
        ax2.set_ylabel("p_tilde [mmHg]")
        h2, l2 = ax2.get_legend_handles_labels()
        handles += h2
        labels += l2

    ax.set_title("Multiplier and interstitial pressure")
    ax.legend(handles, labels, fontsize=7, loc="lower left")


def cycle_overlay(ax, df, col, title, ylabel=""):
    """Overlay every cardiac cycle against the in-cycle time tau."""
    if col not in df or "tau" not in df or "cycle" not in df:
        ax.set_axis_off()
        return
    shade_systole_tau(ax)
    for cyc, grp in df.groupby("cycle"):
        if grp.empty:
            continue
        ax.plot(grp["tau"], numeric_series(grp, col), lw=1.2, label=f"cycle {cyc}")
    ax.set_title(title)
    ax.set_xlabel("tau [s]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)


def shade_systole_tau(ax):
    """Shade systole on an in-cycle (tau) axis."""
    ax.axvspan(SYSTOLE[0], SYSTOLE[1], color="0.85", zorder=0, lw=0)


def savefig(path):
    """Save current figure and close it."""
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def ts(ax, df, cols, title, ylabel=""):
    """Plot one or more time series."""
    cols = cols_present(df, cols)

    if "t" not in df:
        ax.set_visible(False)
        return

    t = numeric_series(df, "t")

    for c in cols:
        y = numeric_series(df, c)
        valid = t.notna() & y.notna()
        ax.plot(t[valid], y[valid], label=c)

    ax.set_title(title)
    ax.set_xlabel("t [s]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)

    if cols:
        ax.legend(fontsize=7)


def scatter(ax, df, x, y, title, xlabel=None, ylabel=None):
    """Scatter plot for input-output or validation relationships."""
    if x not in df or y not in df:
        ax.set_visible(False)
        return

    xs = numeric_series(df, x)
    ys = numeric_series(df, y)
    valid = xs.notna() & ys.notna()

    ax.scatter(xs[valid], ys[valid], s=8, alpha=0.75)
    ax.set_title(title)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    ax.grid(True, alpha=0.3)


def line_xy(ax, df, x, y, title, xlabel=None, ylabel=None):
    """Plot y(x), useful for PV loops and phase portraits."""
    if x not in df or y not in df:
        ax.set_visible(False)
        return

    xs = numeric_series(df, x)
    ys = numeric_series(df, y)
    valid = xs.notna() & ys.notna()

    ax.plot(xs[valid], ys[valid])
    ax.set_title(title)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    ax.grid(True, alpha=0.3)


def pv_loop(ax, df):
    """Pressure-volume loop, with the last complete cycle highlighted."""
    if "CavityVolume_mL" not in df or "VentricularPressure_mmHg" not in df:
        ax.set_visible(False)
        return

    line_xy(
        ax,
        df,
        "CavityVolume_mL",
        "VentricularPressure_mmHg",
        "Cavity pressure-volume loop",
        "cavity volume [mL]",
        "p_v [mmHg]",
    )

    grp = last_full_cycle(df)
    if grp is not None and not grp.empty:
        ax.plot(
            numeric_series(grp, "CavityVolume_mL"),
            numeric_series(grp, "VentricularPressure_mmHg"),
            color="C3",
            lw=2.0,
            label="last cycle",
        )
        stats = cycle_haemodynamics(grp)
        if stats is not None:
            ax.set_title(
                "Cavity pressure-volume loop\n"
                f"EDV {stats['EDV_mL']:.1f}, ESV {stats['ESV_mL']:.1f}, "
                f"SV {stats['SV_mL']:.1f} mL, EF {100 * stats['EF']:.1f} %",
                fontsize=9,
            )
        ax.legend(fontsize=7)


def hist(ax, df, col, title=None):
    """Histogram of a diagnostic variable."""
    if col not in df:
        ax.set_visible(False)
        return

    values = safe_series(df, col)

    if values.empty:
        ax.set_title(title or f"Histogram: {col}")
        ax.text(0.5, 0.5, "waiting for numeric data", ha="center", va="center")
        return

    ax.hist(values, bins=60)
    ax.set_title(title or f"Histogram: {col}")
    ax.set_xlabel(col)
    ax.set_ylabel("count")
    ax.grid(True, alpha=0.3)


def corr_matrix(ax, df, cols, title="Correlation matrix"):
    """Plot a correlation matrix for selected diagnostics."""
    cols = cols_present(df, cols)

    if not cols:
        ax.set_visible(False)
        return None

    clean = df[cols].apply(pd.to_numeric, errors="coerce")
    clean = clean.replace([np.inf, -np.inf], np.nan)
    clean = clean.dropna(axis=1, how="all")

    if clean.shape[0] < 2 or clean.shape[1] < 2:
        ax.set_title(title)
        ax.text(
            0.5,
            0.5,
            "waiting for enough numeric data",
            ha="center",
            va="center",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        return None

    corr = clean.corr()

    im = ax.imshow(corr, aspect="auto", vmin=-1.0, vmax=1.0)
    ax.set_title(title)
    ax.set_xticks(range(len(corr.columns)))
    ax.set_yticks(range(len(corr.index)))
    ax.set_xticklabels(corr.columns, rotation=90, fontsize=6)
    ax.set_yticklabels(corr.index, fontsize=6)

    return im


def cycle_haemodynamics(grp):
    """End-diastolic/systolic volumes, stroke volume, ejection fraction, work."""
    if "CavityVolume_mL" not in grp or grp["CavityVolume_mL"].dropna().empty:
        return None

    v = numeric_series(grp, "CavityVolume_mL").dropna()
    edv = float(v.max())
    esv = float(v.min())
    out = {
        "EDV_mL": edv,
        "ESV_mL": esv,
        "SV_mL": edv - esv,
        "EF": (edv - esv) / edv if edv else np.nan,
    }

    if "V" in grp and "pv" in grp:
        vol = numeric_series(grp, "V").to_numpy()
        pres = numeric_series(grp, "pv").to_numpy()
        ok = np.isfinite(vol) & np.isfinite(pres)
        if ok.sum() > 2:
            # Stroke work is the (clockwise) area of the PV loop, -oint p dV
            # with the trapezoidal rule on the closed contour.
            out["StrokeWork_J"] = -float(TRAPEZOID(pres[ok], vol[ok]))

    for col, key in (
        ("VentricularPressure_mmHg", "pv"),
        ("ArterialPressure_mmHg", "par"),
        ("InterstitialPressure_mmHg", "pf"),
    ):
        if col in grp:
            s = safe_series(grp, col)
            if not s.empty:
                out[f"{key}_min_mmHg"] = float(s.min())
                out[f"{key}_max_mmHg"] = float(s.max())

    if "phi" in grp:
        s = safe_series(grp, "phi")
        if not s.empty:
            out["phi_min"] = float(s.min())
            out["phi_max"] = float(s.max())

    if "MultiplierBar_kPa" in grp:
        s = safe_series(grp, "MultiplierBar_kPa")
        if not s.empty:
            out["lambdaBar_min_kPa"] = float(s.min())

    for col, key in (
        ("CoronaryArterialInflow_mL_min", "coronary_inflow"),
        ("CoronaryVenousOutflow_mL_min", "venous_outflow"),
    ):
        if col in grp:
            s = safe_series(grp, col)
            if not s.empty:
                out[f"{key}_mean_mL_min"] = float(s.mean())

    return out


def dashboard_specs():
    """Return the grouped dashboard windows."""
    return [
        (
            "Cavity / LV dynamics",
            [
                lambda ax, df: ts(
                    ax,
                    df,
                    [
                        "AtrialPressure_mmHg",
                        "VentricularPressure_mmHg",
                        "ArterialPressure_mmHg",
                        "DistalPressure_mmHg",
                    ],
                    "Pressures",
                    "mmHg",
                ),
                lambda ax, df: ts(ax, df, ["CavityVolume_mL"], "Cavity volume", "mL"),
                lambda ax, df: ts(ax, df, ["CavityFlow_mL_s"], "Cavity flow", "mL/s"),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["EndocardialDisplacement_mm"],
                    "Endocardial displacement",
                    "mm",
                ),
                lambda ax, df: pv_loop(ax, df),
                lambda ax, df: cycle_overlay(
                    ax, df, "VentricularPressure_mmHg",
                    "Cavity pressure, cycles overlaid", "mmHg",
                ),
            ],
        ),
        (
            "Wall poromechanics",
            [
                lambda ax, df: ts(ax, df, ["phi"], "Lagrangian porosity", "-"),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["VolumeRatio_J"],
                    "Volume ratio  J = 1 - phi0 + Phi",
                    "-",
                ),
                lambda ax, df: impediment_panel(ax, df),
                lambda ax, df: ts(
                    ax,
                    df,
                    [
                        "VentricularPressure_mmHg",
                        "ArterialPressure_mmHg",
                        "InterstitialPressure_mmHg",
                    ],
                    "Cavity, arterial and interstitial pressure",
                    "mmHg",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["WallBloodVolume_mL", "WallVolume_mL"],
                    "Wall blood volume and wall volume",
                    "mL",
                ),
                lambda ax, df: cycle_overlay(
                    ax, df, "MultiplierBar_kPa",
                    "Averaged multiplier, cycles overlaid", "kPa",
                ),
            ],
        ),
        (
            "Coronary perfusion",
            [
                lambda ax, df: phased_flows(ax, df),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["PerfusionPressure_mmHg"],
                    "Perfusion pressure  p_ar - p_tilde",
                    "mmHg",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["InterstitialMinusCavity"],
                    "Interstitial minus cavity pressure  p_tilde - p_v",
                    "Pa",
                ),
                lambda ax, df: cycle_overlay(
                    ax, df, "CoronaryArterialInflow_mL_min",
                    "Coronary inflow, cycles overlaid", "mL/min",
                ),
                lambda ax, df: cycle_overlay(
                    ax, df, "CoronaryVenousOutflow_mL_min",
                    "Venous outflow, cycles overlaid", "mL/min",
                ),
                lambda ax, df: scatter(
                    ax,
                    df,
                    "PerfusionPressure",
                    "CoronaryArterialInflow_mL_min",
                    "Perfusion pressure vs coronary inflow",
                    "p_ar - p_tilde [Pa]",
                    "inflow [mL/min]",
                ),
            ],
        ),
        (
            "Validation / conservation",
            [
                lambda ax, df: ts(
                    ax,
                    df,
                    ["MassBalanceResidual"],
                    "Fluid mass-balance residual",
                    "m³/s",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["RelativeMassBalanceResidual"],
                    "Relative mass-balance residual",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [
                        "CoronaryVenousOutflow_mL_min",
                        "CoronaryVenousOutflowFromBalance_mL_min",
                        "WallStorageFlux_mL_min",
                    ],
                    "Venous outflow: law vs mass balance",
                    "mL/min",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["StorageLawResidual"],
                    "Storage-law residual  p_tilde - [K_Phi (Phi - phi0) - lambda_bar]",
                    "Pa",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["CavityFlowMismatch"],
                    "Cavity flow: reported minus BDF2 rate  (= 0.5 dt Vddot)",
                    "m³/s",
                ),
                lambda ax, df: scatter(
                    ax,
                    df,
                    "PorosityChange",
                    "StorageLawResidual",
                    "Storage-law residual vs porosity change",
                    "Phi - phi0",
                    "residual [Pa]",
                ),
            ],
        ),
        (
            "Internal / active variables",
            [
                lambda ax, df: ts(ax, df, ["ec"], "ec"),
                lambda ax, df: ts(ax, df, ["gamma"], "gamma"),
                lambda ax, df: ts(ax, df, ["beta"], "beta"),
                lambda ax, df: ts(ax, df, ["kc"], "kc"),
                lambda ax, df: ts(ax, df, ["tauc"], "tauc"),
                lambda ax, df: ts(ax, df, ["w"], "w  (load-dependent relaxation)"),
            ],
        ),
        (
            "Histograms",
            [
                lambda ax, df: hist(ax, df, "VentricularPressure_mmHg"),
                lambda ax, df: hist(ax, df, "ArterialPressure_mmHg"),
                lambda ax, df: hist(ax, df, "CavityVolume_mL"),
                lambda ax, df: hist(ax, df, "phi"),
                lambda ax, df: hist(ax, df, "MultiplierBar_kPa"),
                lambda ax, df: hist(ax, df, "CoronaryArterialInflow_mL_min"),
            ],
        ),
    ]


def diagnostic_columns(df):
    """Columns used for summaries and correlations."""
    return cols_present(
        df,
        [
            "AtrialPressure_mmHg",
            "VentricularPressure_mmHg",
            "ArterialPressure_mmHg",
            "DistalPressure_mmHg",
            "InterstitialPressure_mmHg",
            "PerfusionPressure_mmHg",
            "EndocardialDisplacement_mm",
            "CavityVolume_mL",
            "CavityFlow_mL_s",
            "phi",
            "VolumeRatio_J",
            "MultiplierBar_kPa",
            "WallBloodVolume_mL",
            "CoronaryArterialInflow_mL_min",
            "CoronaryVenousOutflow_mL_min",
            "RelativeMassBalanceResidual",
            "ec",
            "gamma",
            "beta",
            "kc",
            "tauc",
            "w",
        ],
    )


def create_dashboard_windows():
    """Create all dashboard windows."""
    windows = []

    for title, funcs in dashboard_specs():
        fig, axs = plt.subplots(2, 3, figsize=(17, 9))

        try:
            fig.canvas.manager.set_window_title(title)
        except Exception:
            pass

        fig.suptitle(title, fontsize=14)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        windows.append((title, fig, axs.ravel(), funcs))

    return windows


def update_dashboard_windows(windows, df):
    """Refresh all dashboard windows with a new dataframe."""
    for _, fig, axs, funcs in windows:
        if not plt.fignum_exists(fig.number):
            continue

        for ax in axs:
            ax.clear()
            ax.set_visible(True)

        for ax, func in zip(axs, funcs):
            func(ax, df)

        fig.canvas.draw_idle()
        fig.canvas.flush_events()


def show_cycle_dashboard(df):
    """Show cycle overlays for periodic validation."""
    if "cycle" not in df or "tau" not in df:
        return

    if df["cycle"].dropna().empty or df["cycle"].max() < 1:
        return

    cycle_cols = cols_present(
        df,
        [
            "VentricularPressure_mmHg",
            "ArterialPressure_mmHg",
            "CavityVolume_mL",
            "phi",
            "MultiplierBar_kPa",
            "CoronaryArterialInflow_mL_min",
        ],
    )

    if not cycle_cols:
        return

    fig, axs = plt.subplots(2, 3, figsize=(17, 9))

    try:
        fig.canvas.manager.set_window_title("Cycle overlays")
    except Exception:
        pass

    fig.suptitle("Cycle overlays", fontsize=14)
    axs = axs.ravel()

    for ax, col in zip(axs, cycle_cols):
        shade_systole_tau(ax)
        for cycle, g in df.groupby("cycle"):
            ax.plot(g["tau"], g[col], label=f"cycle {cycle}", alpha=0.8)

        ax.set_title(col)
        ax.set_xlabel("cycle time [s]")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=7)

    fig.tight_layout(rect=(0, 0, 1, 0.96))


def show_correlation_dashboard(df):
    """Show static correlation matrix with colorbar."""
    cols = diagnostic_columns(df)

    if not cols:
        return

    fig, ax = plt.subplots(figsize=(12.5, 9))
    fig.subplots_adjust(left=0.22, right=0.78, bottom=0.30, top=0.90)
    cax = fig.add_axes([0.84, 0.30, 0.025, 0.55])

    try:
        fig.canvas.manager.set_window_title("Correlation matrix")
    except Exception:
        pass

    im = corr_matrix(ax, df, cols, "Diagnostic correlation matrix")

    if im is not None:
        fig.colorbar(im, cax=cax, label="correlation")


def grouped_dashboard(df):
    """Show static grouped dashboards."""
    windows = create_dashboard_windows()
    update_dashboard_windows(windows, df)
    show_cycle_dashboard(df)
    show_correlation_dashboard(df)
    plt.show()


def live_dashboard(csv_path, interval=1.0, save_every=None, outdir=None):
    """
    Monitor a growing CSV and update the full dashboard set live.

    The colorbar is created once using a ScalarMappable and is never recreated
    inside the refresh loop. This avoids matplotlib colorbar recursion on macOS.
    """
    plt.ion()

    windows = create_dashboard_windows()

    corr_fig, corr_ax = plt.subplots(figsize=(12.5, 9))
    corr_fig.subplots_adjust(left=0.22, right=0.78, bottom=0.30, top=0.90)
    corr_cax = corr_fig.add_axes([0.84, 0.30, 0.025, 0.55])

    corr_norm = mpl.colors.Normalize(vmin=-1.0, vmax=1.0)
    corr_sm = mpl.cm.ScalarMappable(norm=corr_norm, cmap=plt.get_cmap())
    corr_sm.set_array([])
    corr_fig.colorbar(corr_sm, cax=corr_cax, label="correlation")

    try:
        corr_fig.canvas.manager.set_window_title("Live: correlation matrix")
    except Exception:
        pass

    refresh = 0

    while any(plt.fignum_exists(fig.number) for _, fig, _, _ in windows):
        df = read_csv_robust(csv_path)

        if df is None:
            time.sleep(interval)
            continue

        update_dashboard_windows(windows, df)

        if plt.fignum_exists(corr_fig.number):
            corr_ax.clear()

            corr_matrix(
                corr_ax,
                df,
                diagnostic_columns(df),
                "Live diagnostic correlation matrix",
            )

            corr_fig.canvas.draw_idle()
            corr_fig.canvas.flush_events()

        if save_every is not None and outdir is not None and refresh % save_every == 0:
            outdir.mkdir(parents=True, exist_ok=True)
            save_live_snapshot(df, outdir)

        refresh += 1
        time.sleep(interval)

    if corr_fig is not None and plt.fignum_exists(corr_fig.number):
        plt.close(corr_fig)

    plt.ioff()


def save_live_snapshot(df, outdir):
    """Save lightweight live snapshot tables."""
    write_tables(df, outdir)


def save_all_figures(df, outdir):
    """Save all individual diagnostic figures."""
    outdir.mkdir(parents=True, exist_ok=True)

    def make_ts(name, cols, title, ylabel=""):
        plt.figure(figsize=(11, 5))
        ts(plt.gca(), df, cols, title, ylabel)
        savefig(outdir / name)

    # ---- Cavity / circulation ------------------------------------------
    make_ts(
        "01_pressures_0d.png",
        [
            "AtrialPressure_mmHg",
            "VentricularPressure_mmHg",
            "ArterialPressure_mmHg",
            "DistalPressure_mmHg",
        ],
        "0D pressure dynamics",
        "mmHg",
    )
    make_ts("02_cavity_volume.png", ["CavityVolume_mL"], "Cavity volume", "mL")
    make_ts("03_cavity_flow.png", ["CavityFlow_mL_s"], "Cavity flow", "mL/s")
    make_ts(
        "04_endocardial_kinematics.png",
        ["EndocardialDisplacement_mm"],
        "Endocardial displacement",
        "mm",
    )

    plt.figure(figsize=(6.5, 6))
    pv_loop(plt.gca(), df)
    savefig(outdir / "05_pv_loop.png")

    # ---- Wall poromechanics --------------------------------------------
    make_ts("06_porosity.png", ["phi"], "Lagrangian porosity", "-")
    make_ts(
        "07_volume_ratio.png",
        ["VolumeRatio_J"],
        "Volume ratio  J = 1 - phi0 + Phi",
        "-",
    )
    make_ts(
        "08_wall_volumes.png",
        ["WallBloodVolume_mL", "WallVolume_mL"],
        "Wall blood volume and wall volume",
        "mL",
    )

    fig, ax = plt.subplots(figsize=(11, 5))
    impediment_panel(ax, df)
    savefig(outdir / "09_multiplier_interstitial.png")

    make_ts(
        "10_pressures_wall.png",
        [
            "VentricularPressure_mmHg",
            "ArterialPressure_mmHg",
            "InterstitialPressure_mmHg",
        ],
        "Cavity, arterial and interstitial pressure",
        "mmHg",
    )

    # ---- Coronary perfusion --------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 5))
    phased_flows(ax, df)
    savefig(outdir / "11_coronary_inflow_vs_venous.png")

    make_ts(
        "12_perfusion_pressure.png",
        ["PerfusionPressure_mmHg"],
        "Perfusion pressure  p_ar - p_tilde",
        "mmHg",
    )

    for name, col, ylab in (
        ("13_coronary_inflow_cycles.png", "CoronaryArterialInflow_mL_min", "mL/min"),
        ("14_venous_outflow_cycles.png", "CoronaryVenousOutflow_mL_min", "mL/min"),
        ("15_porosity_cycles.png", "phi", "-"),
        ("16_multiplier_cycles.png", "MultiplierBar_kPa", "kPa"),
    ):
        if col in df:
            fig, ax = plt.subplots(figsize=(8, 4.5))
            cycle_overlay(ax, df, col, col, ylab)
            savefig(outdir / name)

    # ---- Validation ----------------------------------------------------
    make_ts(
        "17_mass_balance.png",
        ["MassBalanceResidual"],
        "Fluid mass-balance residual",
        "m³/s",
    )
    make_ts(
        "18_relative_mass_balance.png",
        ["RelativeMassBalanceResidual"],
        "Relative mass-balance residual",
    )
    make_ts(
        "19_venous_law_vs_balance.png",
        [
            "CoronaryVenousOutflow_mL_min",
            "CoronaryVenousOutflowFromBalance_mL_min",
            "WallStorageFlux_mL_min",
        ],
        "Venous outflow: law vs mass balance",
        "mL/min",
    )
    make_ts(
        "20_storage_law_residual.png",
        ["StorageLawResidual"],
        "Storage-law residual",
        "Pa",
    )
    make_ts(
        "21_cavity_flow_mismatch.png",
        ["CavityFlowMismatch"],
        "Cavity flow: reported minus BDF2 rate  (= 0.5 dt Vddot)",
        "m³/s",
    )

    plt.figure(figsize=(6, 6))
    scatter(
        plt.gca(),
        df,
        "PerfusionPressure",
        "CoronaryArterialInflow_mL_min",
        "Perfusion pressure vs coronary inflow",
        "p_ar - p_tilde [Pa]",
        "inflow [mL/min]",
    )
    savefig(outdir / "22_perfusion_pressure_vs_inflow.png")

    plt.figure(figsize=(6, 6))
    scatter(
        plt.gca(),
        df,
        "MultiplierBar_kPa",
        "CoronaryArterialInflow_mL_min",
        "Multiplier vs coronary inflow",
        "lambda_bar [kPa]",
        "inflow [mL/min]",
    )
    savefig(outdir / "23_multiplier_vs_inflow.png")

    # ---- Active variables ----------------------------------------------
    make_ts(
        "24_internal_variables.png",
        ["ec", "gamma", "beta", "kc", "tauc", "w"],
        "Internal variables",
    )

    hist_dir = outdir / "histograms"
    hist_dir.mkdir(exist_ok=True)

    for col in [
        "VentricularPressure_mmHg",
        "ArterialPressure_mmHg",
        "DistalPressure_mmHg",
        "InterstitialPressure_mmHg",
        "CavityVolume_mL",
        "CavityFlow_mL_s",
        "phi",
        "MultiplierBar_kPa",
        "CoronaryArterialInflow_mL_min",
        "CoronaryVenousOutflow_mL_min",
        "RelativeMassBalanceResidual",
    ]:
        if col in df:
            plt.figure(figsize=(7, 5))
            hist(plt.gca(), df, col)
            savefig(hist_dir / f"hist_{col}.png")


def write_tables(df, outdir):
    """Write summary, per-cycle, correlation, and validation tables."""
    outdir.mkdir(parents=True, exist_ok=True)

    diagnostic_cols = diagnostic_columns(df)

    if diagnostic_cols:
        clean = df[diagnostic_cols].apply(pd.to_numeric, errors="coerce")
        clean = clean.replace([np.inf, -np.inf], np.nan)

        summary = clean.agg(["min", "max", "mean", "std", "median"]).T
        summary["abs_max"] = clean.abs().max()
        summary.to_csv(outdir / "summary.csv")
        print(summary)

        corr = clean.corr()
        corr.to_csv(outdir / "correlation_matrix.csv")

        fig, ax = plt.subplots(figsize=(12.5, 9))
        fig.subplots_adjust(left=0.22, right=0.78, bottom=0.30, top=0.90)
        cax = fig.add_axes([0.84, 0.30, 0.025, 0.55])

        im = corr_matrix(ax, df, diagnostic_cols, "Diagnostic correlation matrix")

        if im is not None:
            fig.colorbar(im, cax=cax, label="correlation")

        plt.savefig(outdir / "correlation_matrix.png", dpi=200)
        plt.close(fig)

    # ---- Per-cycle haemodynamics ---------------------------------------
    if "cycle" in df and not df["cycle"].dropna().empty:
        rows = []
        for cycle, grp in df.groupby("cycle"):
            stats = cycle_haemodynamics(grp)
            if stats is None:
                continue
            stats["cycle"] = int(cycle)
            rows.append(stats)
        if rows:
            cycles = pd.DataFrame(rows).set_index("cycle").sort_index()
            cycles.to_csv(outdir / "cycles.csv")
            print(cycles)

    flags = []

    def add_flag(name, value):
        flags.append({"diagnostic": name, "value": value})

    # Haemodynamics of the last complete cycle: the periodic regime is what
    # should be compared against physiology, not the initial transient.
    grp = last_full_cycle(df)
    if grp is not None:
        stats = cycle_haemodynamics(grp)
        if stats is not None:
            for k, v in stats.items():
                add_flag(f"last_cycle_{k}", v)

    # Perfusion phasing. Coronary inflow should be diastolic dominant and
    # venous outflow systolic dominant; a systolic/diastolic ratio below 1 for
    # the inflow and above 1 for the outflow is the systolic impediment.
    for col, label in (
        ("CoronaryArterialInflow_mL_min", "coronary_inflow"),
        ("CoronaryVenousOutflow_mL_min", "venous_outflow"),
    ):
        st = phase_split(df, col)
        if st is None:
            continue
        add_flag(f"{label}_mean_mL_min", st["mean"])
        add_flag(f"{label}_systolic_mean_mL_min", st["systolic_mean"])
        add_flag(f"{label}_diastolic_mean_mL_min", st["diastolic_mean"])
        add_flag(f"{label}_peak_over_mean", st["peak_over_mean"])
        if np.isfinite(st["diastolic_mean"]) and st["diastolic_mean"] != 0.0:
            add_flag(
                f"{label}_systolic_over_diastolic",
                st["systolic_mean"] / st["diastolic_mean"],
            )

    # Conservation and admissibility.
    if "RelativeMassBalanceResidual" in df:
        add_flag(
            "max_relative_mass_balance_residual",
            safe_series(df, "RelativeMassBalanceResidual").max(),
        )

    if "MassBalanceResidual" in df:
        add_flag(
            "max_abs_mass_balance_residual",
            safe_series(df, "MassBalanceResidual").abs().max(),
        )

    if "StorageLawResidual" in df:
        add_flag(
            "max_abs_storage_law_residual",
            safe_series(df, "StorageLawResidual").abs().max(),
        )

    if "CavityFlowMismatch" in df:
        add_flag(
            "max_abs_cavity_flow_mismatch",
            safe_series(df, "CavityFlowMismatch").abs().max(),
        )

    if "phi" in df:
        s = safe_series(df, "phi")
        add_flag("min_porosity", s.min())
        add_flag("max_porosity", s.max())

    if "VolumeRatio_J" in df:
        s = safe_series(df, "VolumeRatio_J")
        add_flag("min_volume_ratio_J", s.min())
        add_flag("max_volume_ratio_J", s.max())

    if "lambdaBar" in df:
        s = safe_series(df, "lambdaBar")
        add_flag("min_lambda_bar_Pa", s.min())
        add_flag("max_lambda_bar_Pa", s.max())

    if "pf" in df:
        s = safe_series(df, "pf")
        add_flag("min_interstitial_pressure_Pa", s.min())
        add_flag("max_interstitial_pressure_Pa", s.max())

    if "pv" in df:
        s = safe_series(df, "pv")
        add_flag("min_cavity_pressure_Pa", s.min())
        add_flag("max_cavity_pressure_Pa", s.max())

    if "par" in df:
        s = safe_series(df, "par")
        add_flag("min_arterial_pressure_Pa", s.min())
        add_flag("max_arterial_pressure_Pa", s.max())

    pd.DataFrame(flags).to_csv(outdir / "validation_flags.csv", index=False)


def main():
    global PERIOD, SYSTOLE

    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("poroelastic_diagnostics"))
    parser.add_argument("--show", action="store_true", help="show grouped dashboard windows")
    parser.add_argument("--watch", action="store_true", help="monitor CSV in realtime")
    parser.add_argument("--interval", type=float, default=1.0, help="live refresh interval in seconds")
    parser.add_argument("--no-save", action="store_true", help="do not save PNG/CSV diagnostics")
    parser.add_argument("--save-every", type=int, default=None, help="in live mode, save tables every N refreshes")

    parser.add_argument("--period", type=float, default=PERIOD, help="cardiac period [s]")
    parser.add_argument("--systole", type=float, nargs=2, default=list(SYSTOLE),
                        metavar=("START", "END"), help="systolic window within a cycle [s]")

    # Model parameters not carried by the CSV; keep in sync with the run.
    parser.add_argument("--R0", type=float, default=PARAMS["R0"], help="reference cavity radius [m]")
    parser.add_argument("--d0", type=float, default=PARAMS["d0"], help="reference wall thickness [m]")
    parser.add_argument("--phi0", type=float, default=PARAMS["phi0"], help="reference porosity")
    parser.add_argument("--KPhi", type=float, default=PARAMS["KPhi"], help="fluid storage modulus [Pa]")
    parser.add_argument("--gamma-ar", type=float, default=PARAMS["gammaAr"],
                        help="arterial perfusion conductance [m^3/(s Pa)]")
    parser.add_argument("--gamma-ven", type=float, default=PARAMS["gammaVen"],
                        help="venous perfusion conductance [m^3/(s Pa)]")
    parser.add_argument("--psv", type=float, default=PARAMS["pSv"], help="venous pressure [Pa]")

    args = parser.parse_args()

    PERIOD = args.period
    SYSTOLE = (args.systole[0], args.systole[1])
    PARAMS.update(
        {
            "R0": args.R0,
            "d0": args.d0,
            "phi0": args.phi0,
            "KPhi": args.KPhi,
            "gammaAr": args.gamma_ar,
            "gammaVen": args.gamma_ven,
            "pSv": args.psv,
        }
    )

    if args.watch:
        live_dashboard(
            args.csv,
            interval=args.interval,
            save_every=args.save_every,
            outdir=None if args.no_save else args.outdir,
        )
        return

    df = pd.read_csv(args.csv, dtype=str)
    df = add_derived_columns(df)

    if not args.no_save:
        args.outdir.mkdir(parents=True, exist_ok=True)
        save_all_figures(df, args.outdir)
        write_tables(df, args.outdir)
        print(f"Wrote diagnostics to {args.outdir}")

    if args.show:
        grouped_dashboard(df)


if __name__ == "__main__":
    main()
