from scipy.interpolate import griddata
from matplotlib import pyplot as plt
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


