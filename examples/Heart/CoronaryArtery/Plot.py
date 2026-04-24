#!/usr/bin/env python3
"""
Live/static diagnostics dashboard for the coupled 0D--3D coronary simulation.

This script reads the CSV produced by CoupledLV0DCoronary3D and provides:

1. Static diagnostics:
   - saves individual PNG figures,
   - writes summary tables,
   - writes validation flags,
   - writes a correlation matrix.

2. Interactive dashboard:
   - grouped windows for 0D/LV dynamics,
   - coronary/RCR dynamics,
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

Examples
--------
Save all diagnostics:

    python Plot.py CoronaryArtery.csv

Show static dashboard:

    python Plot.py CoronaryArtery.csv --show

Live dashboard:

    python Plot.py CoronaryArtery.csv --watch --interval 1

Live dashboard without saving:

    python Plot.py CoronaryArtery.csv --watch --no-save

Periodically save tables while live-monitoring:

    python Plot.py CoronaryArtery.csv --watch --save-every 20
"""

import argparse
import time
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


OUTLET_IDS = [4, 5, 6, 7, 8, 9]
PERIOD = 0.85


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


def add_derived_columns(df):
    """
    Add derived physical and numerical diagnostics.

    Unit conversions:
    - volume: m^3 -> mL
    - flow: m^3/s -> mL/s

    Added diagnostics:
    - outlet flux sum check,
    - outlet flux mismatch,
    - absolute flow balance,
    - relative flow balance,
    - mean/spread outlet pressures,
    - mean/spread capacitor pressures,
    - cardiac cycle index,
    - cycle-local time.
    """
    df = df.copy()

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    outlet_fluxes = [f"CoronaryOutlet{i}Flux" for i in OUTLET_IDS]
    outlet_pressures = [f"CoronaryOutlet{i}Pressure" for i in OUTLET_IDS]
    cap_pressures = [f"CoronaryOutlet{i}CapPressure" for i in OUTLET_IDS]

    if "LeftVentricleVolume" in df:
        df["LeftVentricleVolume_mL"] = 1e6 * df["LeftVentricleVolume"]

    if "LeftVentricleFlow" in df:
        df["LeftVentricleFlow_mL_s"] = 1e6 * df["LeftVentricleFlow"]

    if "CoronaryInletFlux" in df:
        df["CoronaryInletFlux_mL_s"] = 1e6 * df["CoronaryInletFlux"]
        df["AbsCoronaryInletFlux"] = df["CoronaryInletFlux"].abs()

    if "CoronaryOutletFluxTotal" in df:
        df["CoronaryOutletFluxTotal_mL_s"] = 1e6 * df["CoronaryOutletFluxTotal"]
        df["AbsCoronaryOutletFluxTotal"] = df["CoronaryOutletFluxTotal"].abs()

    if all(c in df for c in outlet_fluxes) and "CoronaryOutletFluxTotal" in df:
        df["OutletFluxSumCheck"] = df[outlet_fluxes].sum(axis=1)
        df["OutletFluxMismatch"] = (
            df["OutletFluxSumCheck"] - df["CoronaryOutletFluxTotal"]
        )

        total = df["CoronaryOutletFluxTotal"].replace(0.0, np.nan)

        for c in outlet_fluxes:
            df[c + "_mL_s"] = 1e6 * df[c]
            df[c + "Fraction"] = df[c] / total

    if "FlowBalance" in df:
        df["AbsFlowBalance"] = df["FlowBalance"].abs()

    if "AbsFlowBalance" in df and "AbsCoronaryInletFlux" in df:
        denom = df["AbsCoronaryInletFlux"].replace(0.0, np.nan)
        df["RelativeFlowBalance"] = df["AbsFlowBalance"] / denom

    if all(c in df for c in outlet_pressures):
        df["MeanOutletPressure"] = df[outlet_pressures].mean(axis=1)
        df["OutletPressureSpread"] = (
            df[outlet_pressures].max(axis=1)
            - df[outlet_pressures].min(axis=1)
        )

    if all(c in df for c in cap_pressures):
        df["MeanCapPressure"] = df[cap_pressures].mean(axis=1)
        df["CapPressureSpread"] = (
            df[cap_pressures].max(axis=1)
            - df[cap_pressures].min(axis=1)
        )

    if "t" in df:
        df["cycle"] = np.floor(df["t"] / PERIOD).astype("Int64")
        df["tau"] = df["t"] - PERIOD * df["cycle"].astype(float)

    return df


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


def dashboard_specs():
    """Return the original grouped dashboard windows."""
    return [
        (
            "0D / LV dynamics",
            [
                lambda ax, df: ts(
                    ax,
                    df,
                    [
                        "LeftAtriumPressure",
                        "LeftVentriclePressure",
                        "AortaPressure",
                        "DistalPressure",
                    ],
                    "Pressures",
                    "Pa",
                ),
                lambda ax, df: ts(ax, df, ["LeftVentricleVolume_mL"], "LV volume", "mL"),
                lambda ax, df: ts(ax, df, ["LeftVentricleFlow_mL_s"], "LV flow", "mL/s"),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["LeftVentricleDisplacement", "LeftVentricleVelocity"],
                    "LV displacement / velocity",
                ),
                lambda ax, df: line_xy(
                    ax,
                    df,
                    "LeftVentricleVolume_mL",
                    "LeftVentriclePressure",
                    "LV pressure-volume loop",
                    "LV volume [mL]",
                    "LV pressure [Pa]",
                ),
                lambda ax, df: line_xy(
                    ax,
                    df,
                    "LeftVentricleDisplacement",
                    "LeftVentricleVelocity",
                    "LV phase portrait",
                    "displacement",
                    "velocity",
                ),
            ],
        ),
        (
            "Coronary / RCR dynamics",
            [
                lambda ax, df: ts(
                    ax,
                    df,
                    ["CoronaryInletFlux_mL_s", "CoronaryOutletFluxTotal_mL_s"],
                    "Coronary total fluxes",
                    "mL/s",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [f"CoronaryOutlet{i}Flux_mL_s" for i in OUTLET_IDS],
                    "Individual outlet fluxes",
                    "mL/s",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [f"CoronaryOutlet{i}FluxFraction" for i in OUTLET_IDS],
                    "Outlet flux fractions",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [f"CoronaryOutlet{i}Pressure" for i in OUTLET_IDS],
                    "Outlet pressures",
                    "Pa",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [f"CoronaryOutlet{i}CapPressure" for i in OUTLET_IDS],
                    "Capacitor pressures",
                    "Pa",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [
                        "MeanOutletPressure",
                        "MeanCapPressure",
                        "OutletPressureSpread",
                        "CapPressureSpread",
                    ],
                    "RCR pressure summary",
                    "Pa",
                ),
            ],
        ),
        (
            "Validation / conservation",
            [
                lambda ax, df: ts(
                    ax,
                    df,
                    ["FlowBalance", "OutletFluxMismatch"],
                    "Absolute flux-balance errors",
                    "m³/s",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    ["RelativeFlowBalance"],
                    "Relative flow-balance error",
                ),
                lambda ax, df: scatter(
                    ax,
                    df,
                    "CoronaryInletFlux_mL_s",
                    "CoronaryOutletFluxTotal_mL_s",
                    "Inlet vs total outlet flux",
                    "inlet [mL/s]",
                    "outlet [mL/s]",
                ),
                lambda ax, df: scatter(
                    ax,
                    df,
                    "AortaPressure",
                    "CoronaryInletFlux_mL_s",
                    "Aorta pressure vs coronary inlet flux",
                    "Aorta pressure [Pa]",
                    "inlet flux [mL/s]",
                ),
                lambda ax, df: scatter(
                    ax,
                    df,
                    "MeanOutletPressure",
                    "CoronaryOutletFluxTotal_mL_s",
                    "Mean outlet pressure vs outlet flux",
                    "Mean outlet pressure [Pa]",
                    "outlet flux [mL/s]",
                ),
                lambda ax, df: ts(
                    ax,
                    df,
                    [
                        "AbsFlowBalance",
                        "AbsCoronaryInletFlux",
                        "AbsCoronaryOutletFluxTotal",
                    ],
                    "Absolute flux magnitudes",
                    "m³/s",
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
                lambda ax, df: ts(
                    ax,
                    df,
                    ["ec", "gamma", "beta", "kc", "tauc"],
                    "All internal variables",
                ),
            ],
        ),
        (
            "Histograms",
            [
                lambda ax, df: hist(ax, df, "LeftVentriclePressure"),
                lambda ax, df: hist(ax, df, "AortaPressure"),
                lambda ax, df: hist(ax, df, "LeftVentricleVolume_mL"),
                lambda ax, df: hist(ax, df, "CoronaryInletFlux_mL_s"),
                lambda ax, df: hist(ax, df, "CoronaryOutletFluxTotal_mL_s"),
                lambda ax, df: hist(ax, df, "RelativeFlowBalance"),
            ],
        ),
    ]


def diagnostic_columns(df):
    """Columns used for summaries and correlations."""
    return cols_present(
        df,
        [
            "LeftAtriumPressure",
            "LeftVentriclePressure",
            "AortaPressure",
            "DistalPressure",
            "LeftVentricleDisplacement",
            "LeftVentricleVelocity",
            "LeftVentricleVolume_mL",
            "LeftVentricleFlow_mL_s",
            "CoronaryInletFlux_mL_s",
            "CoronaryOutletFluxTotal_mL_s",
            "FlowBalance",
            "RelativeFlowBalance",
            "MeanOutletPressure",
            "MeanCapPressure",
            "OutletPressureSpread",
            "CapPressureSpread",
            "ec",
            "gamma",
            "beta",
            "kc",
            "tauc",
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
            "LeftVentriclePressure",
            "AortaPressure",
            "LeftVentricleVolume_mL",
            "CoronaryInletFlux_mL_s",
            "CoronaryOutletFluxTotal_mL_s",
            "FlowBalance",
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
    """Save all original individual diagnostic figures."""
    outdir.mkdir(parents=True, exist_ok=True)

    def make_ts(name, cols, title, ylabel=""):
        plt.figure(figsize=(11, 5))
        ts(plt.gca(), df, cols, title, ylabel)
        savefig(outdir / name)

    make_ts(
        "01_pressures_0d.png",
        ["LeftAtriumPressure", "LeftVentriclePressure", "AortaPressure", "DistalPressure"],
        "0D pressure dynamics",
        "Pa",
    )
    make_ts("02_lv_volume.png", ["LeftVentricleVolume_mL"], "LV volume", "mL")
    make_ts("03_lv_flow.png", ["LeftVentricleFlow_mL_s"], "LV flow", "mL/s")
    make_ts(
        "04_lv_kinematics.png",
        ["LeftVentricleDisplacement", "LeftVentricleVelocity"],
        "LV kinematics",
    )
    make_ts(
        "05_coronary_total_fluxes.png",
        ["CoronaryInletFlux_mL_s", "CoronaryOutletFluxTotal_mL_s"],
        "Coronary total fluxes",
        "mL/s",
    )
    make_ts(
        "06_outlet_fluxes.png",
        [f"CoronaryOutlet{i}Flux_mL_s" for i in OUTLET_IDS],
        "Individual outlet fluxes",
        "mL/s",
    )
    make_ts(
        "07_outlet_flux_fractions.png",
        [f"CoronaryOutlet{i}FluxFraction" for i in OUTLET_IDS],
        "Outlet flux fractions",
    )
    make_ts(
        "08_outlet_pressures.png",
        [f"CoronaryOutlet{i}Pressure" for i in OUTLET_IDS],
        "RCR outlet pressures",
        "Pa",
    )
    make_ts(
        "09_capacitor_pressures.png",
        [f"CoronaryOutlet{i}CapPressure" for i in OUTLET_IDS],
        "RCR capacitor pressures",
        "Pa",
    )
    make_ts(
        "10_flow_balance.png",
        ["FlowBalance", "OutletFluxMismatch"],
        "Flux balance",
        "m³/s",
    )
    make_ts(
        "11_relative_balance.png",
        ["RelativeFlowBalance"],
        "Relative flow-balance error",
    )
    make_ts(
        "12_internal_variables.png",
        ["ec", "gamma", "beta", "kc", "tauc"],
        "Internal variables",
    )

    plt.figure(figsize=(6, 6))
    line_xy(
        plt.gca(),
        df,
        "LeftVentricleVolume_mL",
        "LeftVentriclePressure",
        "LV pressure-volume loop",
        "LV volume [mL]",
        "LV pressure [Pa]",
    )
    savefig(outdir / "13_lv_pv_loop.png")

    plt.figure(figsize=(6, 6))
    line_xy(
        plt.gca(),
        df,
        "LeftVentricleDisplacement",
        "LeftVentricleVelocity",
        "LV phase portrait",
        "displacement",
        "velocity",
    )
    savefig(outdir / "14_lv_phase_portrait.png")

    plt.figure(figsize=(6, 6))
    scatter(
        plt.gca(),
        df,
        "AortaPressure",
        "CoronaryInletFlux_mL_s",
        "Aorta pressure vs coronary inlet flux",
        "Aorta pressure [Pa]",
        "flux [mL/s]",
    )
    savefig(outdir / "15_aorta_pressure_vs_inlet_flux.png")

    plt.figure(figsize=(6, 6))
    scatter(
        plt.gca(),
        df,
        "CoronaryInletFlux_mL_s",
        "CoronaryOutletFluxTotal_mL_s",
        "Inlet vs outlet flux",
        "inlet [mL/s]",
        "outlet [mL/s]",
    )
    savefig(outdir / "16_inlet_vs_outlet_flux.png")

    hist_dir = outdir / "histograms"
    hist_dir.mkdir(exist_ok=True)

    for col in [
        "LeftVentriclePressure",
        "AortaPressure",
        "DistalPressure",
        "LeftVentricleVolume_mL",
        "LeftVentricleFlow_mL_s",
        "CoronaryInletFlux_mL_s",
        "CoronaryOutletFluxTotal_mL_s",
        "FlowBalance",
        "RelativeFlowBalance",
        "MeanOutletPressure",
        "OutletPressureSpread",
    ]:
        if col in df:
            plt.figure(figsize=(7, 5))
            hist(plt.gca(), df, col)
            savefig(hist_dir / f"hist_{col}.png")


def write_tables(df, outdir):
    """Write summary, correlation, and validation tables."""
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

    flags = []

    def add_flag(name, value):
        flags.append({"diagnostic": name, "value": value})

    if "RelativeFlowBalance" in df:
        add_flag("max_relative_flow_balance", safe_series(df, "RelativeFlowBalance").max())

    if "FlowBalance" in df:
        add_flag("max_abs_flow_balance", safe_series(df, "FlowBalance").abs().max())

    if "LeftVentriclePressure" in df:
        add_flag("min_left_ventricle_pressure", safe_series(df, "LeftVentriclePressure").min())
        add_flag("max_left_ventricle_pressure", safe_series(df, "LeftVentriclePressure").max())

    if "AortaPressure" in df:
        add_flag("min_aorta_pressure", safe_series(df, "AortaPressure").min())
        add_flag("max_aorta_pressure", safe_series(df, "AortaPressure").max())

    if "CoronaryOutletFluxTotal" in df:
        add_flag("min_total_outlet_flux", safe_series(df, "CoronaryOutletFluxTotal").min())
        add_flag("max_total_outlet_flux", safe_series(df, "CoronaryOutletFluxTotal").max())

    pd.DataFrame(flags).to_csv(outdir / "validation_flags.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("coronary_diagnostics"))
    parser.add_argument("--show", action="store_true", help="show grouped dashboard windows")
    parser.add_argument("--watch", action="store_true", help="monitor CSV in realtime")
    parser.add_argument("--interval", type=float, default=1.0, help="live refresh interval in seconds")
    parser.add_argument("--no-save", action="store_true", help="do not save PNG/CSV diagnostics")
    parser.add_argument("--save-every", type=int, default=None, help="in live mode, save tables every N refreshes")

    args = parser.parse_args()

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
