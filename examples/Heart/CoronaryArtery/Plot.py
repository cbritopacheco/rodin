#!/usr/bin/env python3
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


OUTLET_IDS = [4, 5, 6, 7, 8, 9]
PERIOD = 0.85


def savefig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def cols_present(df, cols):
    return [c for c in cols if c in df.columns]


def add_derived_columns(df):
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

    if all(c in df for c in outlet_fluxes):
        df["OutletFluxSumCheck"] = df[outlet_fluxes].sum(axis=1)
        df["OutletFluxMismatch"] = df["OutletFluxSumCheck"] - df["CoronaryOutletFluxTotal"]

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
        df["OutletPressureSpread"] = df[outlet_pressures].max(axis=1) - df[outlet_pressures].min(axis=1)

    if all(c in df for c in cap_pressures):
        df["MeanCapPressure"] = df[cap_pressures].mean(axis=1)
        df["CapPressureSpread"] = df[cap_pressures].max(axis=1) - df[cap_pressures].min(axis=1)

    if "t" in df:
        df["cycle"] = np.floor(df["t"] / PERIOD).astype(int)
        df["tau"] = df["t"] - PERIOD * df["cycle"]

    return df


def ts(ax, df, cols, title, ylabel=""):
    cols = cols_present(df, cols)
    for c in cols:
        ax.plot(df["t"], df[c], label=c)
    ax.set_title(title)
    ax.set_xlabel("t [s]")
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    if cols:
        ax.legend(fontsize=7)


def scatter(ax, df, x, y, title, xlabel=None, ylabel=None):
    if x not in df or y not in df:
        ax.set_visible(False)
        return

    ax.scatter(df[x], df[y], s=8, alpha=0.75)
    ax.set_title(title)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    ax.grid(True, alpha=0.3)


def line_xy(ax, df, x, y, title, xlabel=None, ylabel=None):
    if x not in df or y not in df:
        ax.set_visible(False)
        return

    ax.plot(df[x], df[y])
    ax.set_title(title)
    ax.set_xlabel(xlabel or x)
    ax.set_ylabel(ylabel or y)
    ax.grid(True, alpha=0.3)


def hist(ax, df, col, title=None):
    if col not in df:
        ax.set_visible(False)
        return

    ax.hist(df[col].replace([np.inf, -np.inf], np.nan).dropna(), bins=60)
    ax.set_title(title or f"Histogram: {col}")
    ax.set_xlabel(col)
    ax.set_ylabel("count")
    ax.grid(True, alpha=0.3)


def grouped_dashboard(df):
    windows = [
        (
            "0D / LV dynamics",
            [
                lambda ax: ts(ax, df, ["LeftAtriumPressure", "LeftVentriclePressure", "AortaPressure", "DistalPressure"], "Pressures", "Pa"),
                lambda ax: ts(ax, df, ["LeftVentricleVolume_mL"], "LV volume", "mL"),
                lambda ax: ts(ax, df, ["LeftVentricleFlow_mL_s"], "LV flow", "mL/s"),
                lambda ax: ts(ax, df, ["LeftVentricleDisplacement", "LeftVentricleVelocity"], "LV displacement / velocity"),
                lambda ax: line_xy(ax, df, "LeftVentricleVolume_mL", "LeftVentriclePressure", "LV pressure-volume loop", "LV volume [mL]", "LV pressure [Pa]"),
                lambda ax: line_xy(ax, df, "LeftVentricleDisplacement", "LeftVentricleVelocity", "LV phase portrait", "displacement", "velocity"),
            ],
        ),
        (
            "Coronary / RCR dynamics",
            [
                lambda ax: ts(ax, df, ["CoronaryInletFlux_mL_s", "CoronaryOutletFluxTotal_mL_s"], "Coronary total fluxes", "mL/s"),
                lambda ax: ts(ax, df, [f"CoronaryOutlet{i}Flux_mL_s" for i in OUTLET_IDS], "Individual outlet fluxes", "mL/s"),
                lambda ax: ts(ax, df, [f"CoronaryOutlet{i}FluxFraction" for i in OUTLET_IDS], "Outlet flux fractions"),
                lambda ax: ts(ax, df, [f"CoronaryOutlet{i}Pressure" for i in OUTLET_IDS], "Outlet pressures", "Pa"),
                lambda ax: ts(ax, df, [f"CoronaryOutlet{i}CapPressure" for i in OUTLET_IDS], "Capacitor pressures", "Pa"),
                lambda ax: ts(ax, df, ["MeanOutletPressure", "MeanCapPressure", "OutletPressureSpread", "CapPressureSpread"], "RCR pressure summary", "Pa"),
            ],
        ),
        (
            "Validation / conservation",
            [
                lambda ax: ts(ax, df, ["FlowBalance", "OutletFluxMismatch"], "Absolute flux-balance errors", "m³/s"),
                lambda ax: ts(ax, df, ["RelativeFlowBalance"], "Relative flow-balance error"),
                lambda ax: scatter(ax, df, "CoronaryInletFlux_mL_s", "CoronaryOutletFluxTotal_mL_s", "Inlet vs total outlet flux", "inlet [mL/s]", "outlet [mL/s]"),
                lambda ax: scatter(ax, df, "AortaPressure", "CoronaryInletFlux_mL_s", "Aorta pressure vs coronary inlet flux", "Aorta pressure [Pa]", "inlet flux [mL/s]"),
                lambda ax: scatter(ax, df, "MeanOutletPressure", "CoronaryOutletFluxTotal_mL_s", "Mean outlet pressure vs outlet flux", "Mean outlet pressure [Pa]", "outlet flux [mL/s]"),
                lambda ax: ts(ax, df, ["AbsFlowBalance", "AbsCoronaryInletFlux", "AbsCoronaryOutletFluxTotal"], "Absolute flux magnitudes", "m³/s"),
            ],
        ),
        (
            "Internal / active variables",
            [
                lambda ax: ts(ax, df, ["ec"], "ec"),
                lambda ax: ts(ax, df, ["gamma"], "gamma"),
                lambda ax: ts(ax, df, ["beta"], "beta"),
                lambda ax: ts(ax, df, ["kc"], "kc"),
                lambda ax: ts(ax, df, ["tauc"], "tauc"),
                lambda ax: ts(ax, df, ["ec", "gamma", "beta", "kc", "tauc"], "All internal variables"),
            ],
        ),
        (
            "Histograms",
            [
                lambda ax: hist(ax, df, "LeftVentriclePressure"),
                lambda ax: hist(ax, df, "AortaPressure"),
                lambda ax: hist(ax, df, "LeftVentricleVolume_mL"),
                lambda ax: hist(ax, df, "CoronaryInletFlux_mL_s"),
                lambda ax: hist(ax, df, "CoronaryOutletFluxTotal_mL_s"),
                lambda ax: hist(ax, df, "RelativeFlowBalance"),
            ],
        ),
    ]

    for title, funcs in windows:
        fig, axs = plt.subplots(2, 3, figsize=(17, 9))
        try:
            fig.canvas.manager.set_window_title(title)
        except Exception:
            pass
        fig.suptitle(title, fontsize=14)
        axs = axs.ravel()

        for ax, func in zip(axs, funcs):
            func(ax)

        fig.tight_layout(rect=(0, 0, 1, 0.96))

    if "cycle" in df and df["cycle"].max() >= 1:
        cycle_cols = [
            "LeftVentriclePressure",
            "AortaPressure",
            "LeftVentricleVolume_mL",
            "CoronaryInletFlux_mL_s",
            "CoronaryOutletFluxTotal_mL_s",
            "FlowBalance",
        ]
        cycle_cols = cols_present(df, cycle_cols)

        if cycle_cols:
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

    plt.show()


def save_all_figures(df, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    def make_ts(name, cols, title, ylabel=""):
        plt.figure(figsize=(11, 5))
        ax = plt.gca()
        ts(ax, df, cols, title, ylabel)
        savefig(outdir / name)

    make_ts("01_pressures_0d.png", ["LeftAtriumPressure", "LeftVentriclePressure", "AortaPressure", "DistalPressure"], "0D pressure dynamics", "Pa")
    make_ts("02_lv_volume.png", ["LeftVentricleVolume_mL"], "LV volume", "mL")
    make_ts("03_lv_flow.png", ["LeftVentricleFlow_mL_s"], "LV flow", "mL/s")
    make_ts("04_lv_kinematics.png", ["LeftVentricleDisplacement", "LeftVentricleVelocity"], "LV kinematics")
    make_ts("05_coronary_total_fluxes.png", ["CoronaryInletFlux_mL_s", "CoronaryOutletFluxTotal_mL_s"], "Coronary total fluxes", "mL/s")
    make_ts("06_outlet_fluxes.png", [f"CoronaryOutlet{i}Flux_mL_s" for i in OUTLET_IDS], "Individual outlet fluxes", "mL/s")
    make_ts("07_outlet_flux_fractions.png", [f"CoronaryOutlet{i}FluxFraction" for i in OUTLET_IDS], "Outlet flux fractions")
    make_ts("08_outlet_pressures.png", [f"CoronaryOutlet{i}Pressure" for i in OUTLET_IDS], "RCR outlet pressures", "Pa")
    make_ts("09_capacitor_pressures.png", [f"CoronaryOutlet{i}CapPressure" for i in OUTLET_IDS], "RCR capacitor pressures", "Pa")
    make_ts("10_flow_balance.png", ["FlowBalance", "OutletFluxMismatch"], "Flux balance", "m³/s")
    make_ts("11_relative_balance.png", ["RelativeFlowBalance"], "Relative flow-balance error")
    make_ts("12_internal_variables.png", ["ec", "gamma", "beta", "kc", "tauc"], "Internal variables")

    plt.figure(figsize=(6, 6))
    line_xy(plt.gca(), df, "LeftVentricleVolume_mL", "LeftVentriclePressure", "LV pressure-volume loop", "LV volume [mL]", "LV pressure [Pa]")
    savefig(outdir / "13_lv_pv_loop.png")

    plt.figure(figsize=(6, 6))
    line_xy(plt.gca(), df, "LeftVentricleDisplacement", "LeftVentricleVelocity", "LV phase portrait", "displacement", "velocity")
    savefig(outdir / "14_lv_phase_portrait.png")

    plt.figure(figsize=(6, 6))
    scatter(plt.gca(), df, "AortaPressure", "CoronaryInletFlux_mL_s", "Aorta pressure vs coronary inlet flux", "Aorta pressure [Pa]", "flux [mL/s]")
    savefig(outdir / "15_aorta_pressure_vs_inlet_flux.png")

    plt.figure(figsize=(6, 6))
    scatter(plt.gca(), df, "CoronaryInletFlux_mL_s", "CoronaryOutletFluxTotal_mL_s", "Inlet vs outlet flux", "inlet [mL/s]", "outlet [mL/s]")
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
    diagnostic_cols = cols_present(df, [
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
    ])

    if diagnostic_cols:
        summary = df[diagnostic_cols].replace([np.inf, -np.inf], np.nan).agg(
            ["min", "max", "mean", "std", "median"]
        ).T
        summary["abs_max"] = df[diagnostic_cols].replace([np.inf, -np.inf], np.nan).abs().max()
        summary.to_csv(outdir / "summary.csv")
        print(summary)

        corr = df[diagnostic_cols].replace([np.inf, -np.inf], np.nan).corr()
        corr.to_csv(outdir / "correlation_matrix.csv")

        plt.figure(figsize=(11, 9))
        plt.imshow(corr, aspect="auto")
        plt.colorbar(label="correlation")
        plt.xticks(range(len(corr.columns)), corr.columns, rotation=90, fontsize=7)
        plt.yticks(range(len(corr.index)), corr.index, fontsize=7)
        plt.title("Diagnostic correlation matrix")
        savefig(outdir / "correlation_matrix.png")

    flags = []

    def add_flag(name, value):
        flags.append({"diagnostic": name, "value": value})

    if "RelativeFlowBalance" in df:
        add_flag("max_relative_flow_balance", df["RelativeFlowBalance"].replace([np.inf, -np.inf], np.nan).max())

    if "FlowBalance" in df:
        add_flag("max_abs_flow_balance", df["FlowBalance"].abs().max())

    if "LeftVentriclePressure" in df:
        add_flag("min_left_ventricle_pressure", df["LeftVentriclePressure"].min())
        add_flag("max_left_ventricle_pressure", df["LeftVentriclePressure"].max())

    if "AortaPressure" in df:
        add_flag("min_aorta_pressure", df["AortaPressure"].min())
        add_flag("max_aorta_pressure", df["AortaPressure"].max())

    if "CoronaryOutletFluxTotal" in df:
        add_flag("min_total_outlet_flux", df["CoronaryOutletFluxTotal"].min())
        add_flag("max_total_outlet_flux", df["CoronaryOutletFluxTotal"].max())

    pd.DataFrame(flags).to_csv(outdir / "validation_flags.csv", index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("coronary_diagnostics"))
    parser.add_argument("--show", action="store_true", help="show grouped dashboard windows")
    parser.add_argument("--no-save", action="store_true", help="do not save PNG/CSV diagnostics")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
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
