#!/usr/bin/env python3
"""Plot Before.csv and After.csv exported by an MMG example."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('root', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
import sys as _plot_sys
from pathlib import Path as _PlotPath
_plot_examples = next(parent for parent in _PlotPath(__file__).resolve().parents
                      if (parent / 'rodin.mplstyle').is_file())
_plot_sys.path.insert(0, str(_plot_examples))
from plot_style import apply_style
apply_style()
fig, axes = plt.subplots(1, 2, figsize=(9, 4.5), layout='constrained')
for ax, stem in zip(axes, ['Before', 'After']):
    data = np.genfromtxt(args.root / (stem+'.csv'), names=True, delimiter=',')
    triangles = np.stack([data['x0'],data['y0'],data['x1'],data['y1'],
                          data['x2'],data['y2']], axis=1).reshape(-1,3,2)
    ax.add_collection(PolyCollection(triangles, facecolors='#16669933',
                                     edgecolors='#718096', linewidths=.5))
    ax.set(xlim=(-.03,1.03), ylim=(-.03,1.03))
    ax.set_aspect('equal'); ax.set_axis_off()
    ax.set_title(f'{stem}: {len(data)} triangles', pad=12)
fig.savefig(args.output, dpi=180, transparent=True)
