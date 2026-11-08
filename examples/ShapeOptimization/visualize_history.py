#!/usr/bin/env python3
"""Plot the computed cantilever objective history from its exported CSV."""
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style

apply_style()
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('history', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
data = np.genfromtxt(args.history, delimiter=',', names=True)
fig, ax = plt.subplots(figsize=(6, 4.3), layout='constrained')
ax.plot(data['iteration'], data['objective'], 'o-', markersize=3)
ax.set(xlabel='Iteration', ylabel=r'$J(\Omega)$', title='Cantilever objective history')
ax.spines[['top', 'right']].set_visible(False)
args.output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.output, dpi=180, transparent=True)
plt.close(fig)
