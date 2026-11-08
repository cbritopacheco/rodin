#!/usr/bin/env python3
"""Plot the SurfaceMesh refinement CSV files generated at 6, 11 and 21 nodes."""
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, turbo_colors
apply_style()
COLORS = turbo_colors(2)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('root', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
rows = [np.genfromtxt(args.root / f'refinement-{n}' / 'SurfaceMesh.csv',
                      delimiter=',', names=True) for n in [6, 11, 21]]
h = np.array([2 / (float(row['nodes']) - 1) for row in rows])
error = np.array([float(row['relative_l2_error']) for row in rows])
area = np.array([abs(float(row['area']) - 4 * np.pi) / (4 * np.pi) for row in rows])
fig, ax = plt.subplots(figsize=(6.3, 4.2), layout='constrained')
ax.loglog(h, error, 'o-', color=COLORS[0], label='Relative field error')
ax.loglog(h, area, 's-', color=COLORS[1], label='Relative area error')
ax.loglog(h, error[-1] * (h / h[-1]) ** 2, '--', color=COLORS[0], label=r'$h_0^2$ reference')
ax.set(xlabel=r'Background spacing $h_0$', ylabel='Relative error',
       title='Surface field and geometric refinement')
ax.set_xticks(h, labels=['0.4', '0.2', '0.1'])
ax.minorticks_off()
ax.grid(alpha=.2)
ax.legend()
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=180, transparent=True)
