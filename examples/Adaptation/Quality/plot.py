"""Plot the cellwise Rodin samples; retain discontinuities at cell boundaries."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Adaptation_Quality output directory")
    args = parser.parse_args()
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Latin Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": 12,
        "svg.fonttype": "none",
        "pdf.fonttype": 3,
    })
    data = [np.genfromtxt(args.directory / f"P{p}.csv", delimiter=",", names=True)
            for p in (1, 2)]
    j_norm = Normalize(0.8, 1.2)
    q_norm = Normalize(1, 1.51)
    q_map = LinearSegmentedColormap.from_list(
        "distortion", ["#e4f1e7", "#8bbba4", "#ebca85", "#bb553f"])
    fig, axes = plt.subplots(2, 3, figsize=(10, 7),
                             gridspec_kw={"hspace": 0.32, "wspace": 0.12})
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.04, top=0.79)
    for row, samples in enumerate(data):
        degree = row + 1
        assert len(np.unique(samples["cell"])) == 8
        for cell in np.unique(samples["cell"]):
            s = samples[samples["cell"] == cell]
            n = int(s["i"].max())
            lookup = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(s["i"], s["k"]))}
            triangles = []
            for i in range(n):
                for k in range(n - i):
                    triangles.append([lookup[i, k], lookup[i + 1, k], lookup[i, k + 1]])
                    if i + k < n - 1:
                        triangles.append([lookup[i + 1, k], lookup[i + 1, k + 1],
                                          lookup[i, k + 1]])

            # Independent analytic checks of interpolation, geometry and diagnostics.
            if degree == 1:
                slope = 0.2 if s["x0"].mean() < 0.5 else -0.2
                bow = np.interp(s["x0"], [0, 0.5, 1], [0, 0.1, 0])
            else:
                slope = 0.4 * (1 - 2 * s["x0"])
                bow = 0.4 * s["x0"] * (1 - s["x0"])
            expected_j = 1 - 0.5 * slope
            expected_q = (2.25 + slope**2) / (2 * expected_j)
            np.testing.assert_allclose(s["x"], s["x0"] + 0.5 * s["y0"], atol=1e-11)
            np.testing.assert_allclose(s["y"], s["y0"] + bow, atol=1e-11)
            np.testing.assert_allclose(s["j"], expected_j, atol=1e-11)
            np.testing.assert_allclose(s["Q"], expected_q, atol=1e-11)

            triangulation = mtri.Triangulation(s["x"], s["y"], triangles)
            axes[row, 1].tripcolor(triangulation, s["j"], shading="flat",
                                   cmap="Blues", norm=j_norm, rasterized=False,
                                   edgecolors="face", linewidth=0.15)
            axes[row, 2].tripcolor(triangulation, s["Q"], shading="flat",
                                   cmap=q_map, norm=q_norm, rasterized=False,
                                   edgecolors="face", linewidth=0.15)
            for mask, order in ((s["i"] == 0, "k"), (s["k"] == 0, "i"),
                                (s["i"] + s["k"] == n, "i")):
                edge = np.sort(s[mask], order=order)
                axes[row, 0].plot(edge["x0"], edge["y0"], color="#333333", lw=0.8)
                for col in (1, 2):
                    axes[row, col].plot(edge["x"], edge["y"], color="#333333", lw=0.8)
            vertices = s[((s["i"] == 0) & (s["k"] == 0)) |
                         (s["i"] == n) | (s["k"] == n)]
            axes[row, 0].scatter(vertices["x0"], vertices["y0"], s=14, c="#222222", zorder=4)
            for col in (1, 2):
                axes[row, col].scatter(vertices["x"], vertices["y"], s=14, c="#222222", zorder=4)

        print(f"P{degree}: j=[{samples['j'].min():.6g}, {samples['j'].max():.6g}], "
              f"Q=[{samples['Q'].min():.6g}, {samples['Q'].max():.6g}]")
        axes[row, 0].set_ylabel(rf"$P_{degree}$", rotation=0, labelpad=18,
                                fontsize=18, y=0.46, va="center")
        for ax in axes[row]:
            ax.set_aspect("equal")
            ax.set_xlim((-0.33, 1.33) if ax is axes[row, 0] else (-0.08, 1.58))
            ax.set_ylim(-0.08, 1.18)
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
    reference_box = axes[0, 0].get_position()
    fig.text(reference_box.x0 + reference_box.width / 2, 0.965,
             "Background mesh", ha="center", va="center")
    for column, norm, cmap, label in ((1, j_norm, "Blues", r"$j=\det F$"),
                                    (2, q_norm, q_map, r"$Q=\|F\|_F^2/(2j)$")):
        box = axes[0, column].get_position()
        center = box.x0 + box.width / 2
        color_axis = fig.add_axes([box.x0 + 0.05 * box.width, 0.87,
                                  0.9 * box.width, 0.016])
        bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                           cax=color_axis, orientation="horizontal")
        bar.ax.xaxis.set_ticks_position("top")
        bar.ax.tick_params(labelsize=10, pad=6)
        fig.text(center, 0.965, label, ha="center", va="center")
    fig.savefig(args.directory / "quality.png", dpi=220)
    fig.savefig(args.directory / "quality.svg")
    fig.savefig(args.directory / "quality.pdf")
    plt.close(fig)
    print("Analytic checks passed for both degrees.")


if __name__ == "__main__":
    main()
