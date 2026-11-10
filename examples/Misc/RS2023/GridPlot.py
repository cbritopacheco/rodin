from scipy.interpolate import griddata
from matplotlib import pyplot as plt
from pathlib import Path as _PlotPath
import sys as _plot_sys
from pathlib import Path as _PlotPath
_plot_examples = next(parent for parent in _PlotPath(__file__).resolve().parents
                      if (parent / 'rodin.mplstyle').is_file())
_plot_sys.path.insert(0, str(_plot_examples))
from plot_style import apply_style
apply_style()
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str,
                        help='Name of file to plot.')
    args = parser.parse_args()
    X = np.loadtxt(args.filename, delimiter=',')
    plt.hexbin(X[:, 0], X[:, 1], C=X[:, 2], gridsize=50, edgecolors='face',
               linewidths=0.5)
    plt.xlabel('$m$')
    plt.ylabel('$\epsilon$')
    plt.title('$|| u_\epsilon - u_0 ||_{L^2 (Q)}$ when $k = 2$')
    plt.colorbar()
    plt.show()


