"""Export matched oblique views of the refined Y for the bifurcation benchmark schematic."""

import argparse
from pathlib import Path

from report_refinement import draw_geometry, plt


def draw(mesh_directory, output):
    # Fixed limits and figure dimensions register both projected surfaces with
    # the arrows and boundary conditions in y_benchmark_setup.tex.
    xmin, xmax, ymin, ymax = -4.5, 6., -5.8, 4.25
    for domain in ("full", "cut"):
        lines = (mesh_directory / f"{domain}.mesh").read_text().splitlines()
        fig = plt.figure(figsize=((xmax - xmin) * .405 / 2.54, (ymax - ymin) * .405 / 2.54))
        axis = fig.add_axes((0, 0, 1, 1))
        draw_geometry(axis, lines, schematic=True)
        axis.set(xlim=(xmin, xmax), ylim=(ymin, ymax))
        fig.savefig(output / f"y_{domain}_oblique.pdf", transparent=True)
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_directory", type=Path, help="Directory with the refined full.mesh and cut.mesh")
    parser.add_argument("--output", type=Path, default=Path("figures/revision"))
    args = parser.parse_args()
    draw(args.mesh_directory, args.output)
