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
    grid_m, grid_eps = np.mgrid[X[:, 0].min():X[:, 0].max():500j, X[:,
                                                                     1].min():X[:,
                                                                                1].max():500j]
    grid_err = griddata((X[:, 0], X[:, 1]), X[:, 2], (grid_m, grid_eps))
    plt.contourf(grid_m, grid_eps, grid_err)
    # plt.scatter(X[:, 0], X[:, 1], c=X[:, 2], edgecolors='black')
    plt.xlabel('$m$')
    plt.ylabel('$\epsilon$')
    plt.title('$|| u_\epsilon - u_0 ||_{L^2 (Q)}$ when $k = 2$')
    plt.colorbar()
    plt.show()


