"""Compare coronary outlet laws with prescribed pressures in one 0D bed."""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from revision_data import PERIOD, Q_REF, T0, T1, hydraulic_curve, trap


P_RA = 700.0  # Pa
COMPLIANCE = 4e-10  # m^3/Pa
DISTAL_FRACTION = 0.13
REFERENCE_RATES = np.array([800.0, 400.0])  # proximal, distal, in 1/s
CLOSURES = ("R_mu", "R_0", "R_0_replay_fit")
COLORS = ("#0072B2", "#4D4D4D", "#D55E00")
LINE_STYLES = ("-", "--", ":")
LABELS = (r"$\mathcal{R}_\mu$", r"$\mathcal{R}_0$", r"$\mathcal{R}_0^*$")


def pressure_replay(raw):
    """Repeat cycle 3, removing linear endpoint drift about the cycle mean."""
    time = np.linspace(0, PERIOD, 851)
    source_time = np.linspace(T0, T1, time.size)
    pressures = np.array([
        np.interp(source_time, raw["t"], raw[name])
        for name in ("AortaPressure", "IntramyoPressure")
    ])
    drift = pressures[:, -1] - pressures[:, 0]
    pressures -= drift[:, None] * (time / PERIOD - 0.5)
    pressures[:, -1] = pressures[:, 0]
    return CubicSpline(time, pressures, axis=1, bc_type="periodic")


def cycle(forcing, head_ratio, branch_flow, max_step):
    """Return one periodic cycle of the transmural pressure and branch flows."""
    def pressure_drops(t, transmural):
        aortic, myocardial = forcing(t)
        waterfall = np.maximum(myocardial, P_RA)
        inlet = waterfall + head_ratio * (aortic - waterfall)
        compartment = transmural + myocardial
        return np.array([inlet - compartment, compartment - waterfall])

    def balance(t, y):
        proximal, distal = branch_flow(pressure_drops(t, y[0]))
        return [(proximal - distal) / COMPLIANCE]

    def integrate(initial, dense_output=False):
        solution = solve_ivp(
            balance, (0, PERIOD), [initial], method="DOP853", max_step=max_step,
            rtol=5e-10, atol=1e-7, dense_output=dense_output,
        )
        if not solution.success:
            raise RuntimeError(solution.message)
        return solution

    def period_residual(initial):
        return integrate(initial).y[0, -1] - initial

    time = np.linspace(0, PERIOD, 2001)
    boundary_pressures = pressure_drops(time, 0.0) * np.array([[1.0], [-1.0]])
    # Outside the boundary pressure range both resistors drive the state back in.
    initial = brentq(period_residual, boundary_pressures.min() - 1, boundary_pressures.max() + 1, xtol=1e-7)
    solution = integrate(initial, dense_output=True)
    if abs(solution.y[0, -1] - initial) > 1e-4:
        raise RuntimeError("Terminal pressure did not return to its initial state")
    transmural = solution.sol(time)[0]
    return time, transmural, branch_flow(pressure_drops(time, transmural))


def mean(time, values):
    return trap(values, x=time) / PERIOD


def plot_results(path, head_ratios, results, pressure_curve, branch_budgets):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.7), layout="constrained")
    flows = np.geomspace(0.02, 2.0, 400)
    variable_drops = pressure_curve(flows) / branch_budgets[:, None] / flows
    curve_labels = (r"Proximal ($800\,\mathrm{s}^{-1}$)", r"Distal ($400\,\mathrm{s}^{-1}$)")
    for drop, color, label in zip(variable_drops, ("#0072B2", "#5F3A87"), curve_labels):
        axes[0].plot(flows, drop, color=color, lw=2, label=label)
    axes[0].axhline(1.0, c=COLORS[1], ls="--", lw=1.4)
    axes[0].set(xscale="log", xlabel=r"Branch flow $|Q|/Q_0$",
                ylabel=r"Resistance $\mathcal{R}_\mu/\mathcal{R}_0$",
                title="A   Branch resistance")
    axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.13), frameon=False, fontsize=8)

    for axis, ratio, letter in zip(axes[1:], head_ratios, ("B", "C")):
        for (time, _, flow), color, ls, label in zip(results[ratio], COLORS, LINE_STYLES, LABELS):
            axis.plot(time / PERIOD, 1e6 * flow[0], color=color, ls=ls, lw=2, label=label)
        axis.axhline(0, color="0.75", lw=0.8, zorder=0)
        name = "Pressure replay" if ratio == 1.0 else f"Reduced head ({ratio:.0%})"
        axis.set(xlabel="Fraction of cardiac cycle", ylabel="Inflow (mL/s)", title=f"{letter}   {name}")
    axes[1].legend(loc="lower center", bbox_to_anchor=(1.04, 1.13), ncols=3, frameon=False)
    for axis in axes:
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(alpha=0.20)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--low-head-ratio", type=float, default=0.15)
    args = parser.parse_args()
    if not 0.05 <= args.low_head_ratio < 1.0:
        parser.error("--low-head-ratio must be at least 0.05 and below 1")

    raw = np.genfromtxt(Path(__file__).with_name("replay_pressures.csv"), delimiter=",", names=True)
    forcing = pressure_replay(raw)
    # The archive starts at 1 ms. Its first aortic sample approximates p_ar(0).
    head = raw["AortaPressure"][0] - P_RA
    budgets = head * np.array([1.0 - DISTAL_FRACTION, DISTAL_FRACTION])
    resistance = budgets / Q_REF
    stress, rate, _ = hydraulic_curve()
    stress, rate = np.r_[0.0, stress], np.r_[0.0, rate]
    reference_stress = np.interp(REFERENCE_RATES, rate, stress)

    def variable_flow(drops):
        # Wall stress scales with pressure drop at fixed tube geometry.
        # Read the apparent rate at that stress, then scale by the reference flow.
        return np.array([
            np.sign(drop) * (Q_REF / ref_rate) * np.interp(abs(drop) * ref_stress / budget, stress, rate)
            for drop, ref_rate, ref_stress, budget in zip(drops, REFERENCE_RATES, reference_stress, budgets)
        ])

    def pressure_curve(relative_flow):
        return np.array([
            budget / ref_stress * np.interp(ref_rate * relative_flow, rate, stress)
            for budget, ref_stress, ref_rate in zip(budgets, reference_stress, REFERENCE_RATES)
        ])

    np.testing.assert_allclose(variable_flow(budgets), Q_REF, rtol=1e-12)
    step = PERIOD / 200

    def constant(scale):
        def flow(drops):
            return (drops.T / (scale * resistance)).T
        return flow

    def fit_constant(ratio, variable):
        time, _, flows = variable
        target = mean(time, abs(flows[0]))

        def calibration_error(scale):
            time, _, flows = cycle(forcing, ratio, constant(scale), step)
            return (mean(time, abs(flows[0])) - target) / target

        return brentq(calibration_error, 0.25, 4.0, xtol=1e-6)

    head_ratios = (1.0, args.low_head_ratio)
    variable_runs = {ratio: cycle(forcing, ratio, variable_flow, step) for ratio in head_ratios}
    matched_scale = fit_constant(1.0, variable_runs[1.0])
    refit_scale = fit_constant(args.low_head_ratio, variable_runs[args.low_head_ratio])
    fixed = dict(zip(CLOSURES[1:], (1.0, matched_scale)))
    results = {}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    max_refinement = 0.0
    for ratio, variable_run in variable_runs.items():
        fits = fixed.copy()
        if ratio != 1.0:
            fits["R_0_low_head_fit"] = refit_scale
        cases = [("R_mu", None, variable_flow)] + [(name, scale, constant(scale)) for name, scale in fits.items()]
        runs = [variable_run] + [cycle(forcing, ratio, law, step) for _, _, law in cases[1:]]
        results[ratio] = runs[:len(CLOSURES)]
        time, _, variable = variable_run
        mean_absolute_variable = mean(time, abs(variable[0]))
        for (name, scale, law), (_, pressure, (proximal, distal)) in zip(cases, runs):
            _, _, refined = cycle(forcing, ratio, law, step / 2)
            refinement = np.max(abs(refined - [proximal, distal])) / np.max(abs(refined))
            max_refinement = max(max_refinement, refinement)
            if abs(trap(proximal - distal, x=time)) > 1e-5 * trap(abs(proximal), x=time):
                raise RuntimeError("Periodic mass balance check failed")
            distal_rate = REFERENCE_RATES[1] * abs(distal) / Q_REF
            rms_difference = np.sqrt(mean(time, (proximal - variable[0])**2))
            rows.append({
                "head_ratio": ratio,
                "closure": name,
                "resistance_multiplier": scale,
                "mean_net_inflow_ml_s": 1e6 * mean(time, proximal),
                "mean_absolute_inflow_ml_s": 1e6 * mean(time, abs(proximal)),
                "reverse_volume_ml_cycle": 1e6 * trap(np.maximum(-proximal, 0), x=time),
                "normalized_inflow_rms_difference_pct": 100 * rms_difference / mean_absolute_variable,
                "min_transmural_pressure_pa": pressure.min(),
                "max_transmural_pressure_pa": pressure.max(),
                "min_distal_apparent_shear_rate_s-1": distal_rate.min(),
                "max_distal_apparent_shear_rate_s-1": distal_rate.max(),
            })

    if max_refinement > 1e-4:
        raise RuntimeError("Time-step refinement changed a flow waveform by more than 0.01% of its peak")

    with (args.output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    plot_results(args.output_dir / "terminal_comparison.png", head_ratios, results, pressure_curve, budgets)
    print(f"Reference head: {head:.1f} Pa")
    print(f"Resistance multipliers: replay {matched_scale:.5f}, reduced head {refit_scale:.5f}")
    print(f"Largest refinement difference: {100 * max_refinement:.2g}% of peak branch flow")
    print(f"Wrote {args.output_dir / 'summary.csv'} and {args.output_dir / 'terminal_comparison.png'}")


if __name__ == "__main__":
    main()
