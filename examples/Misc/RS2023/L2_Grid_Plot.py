from scipy.interpolate import griddata
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import argparse
import math
from scipy import stats
from sklearn import preprocessing
plt.style.use("bmh")

conductivities = [2.0]
waveNumbers = [1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 20]
angles = [ math.pi / 4 ]


x_column = 'm'
y_column = 'epsilon'

latex_label = {
    'm' : '$m$',
    'epsilon' : '$\epsilon$',
    'waveNumber' : '$k$',
    'angle': '$\theta$'
}

file_label = {
    'm' : 'Period',
    'epsilon' : 'Epsilon',
    'waveNumber' : 'Wavenumber',
    'angle' : 'Angle',
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', type=str,
                        help='Name of file to plot.')
    args = parser.parse_args()
    df = pd.read_csv(args.filename)
    print('Read dataset !')

    thresh = math.exp(df['error'].median())

    df.loc[df['error'] > thresh, 'error'] = thresh

    i = 0
    for waveNumber in waveNumbers:
        print(waveNumber)
        for angle in angles:
            for conductivity in conductivities:
                plt.figure()
                sdf = df[
                    (df['conductivity'] > conductivity - 0.001) &
                    (df['conductivity'] < conductivity + 0.001) &
                    (df['angle'] > angle - 0.001) &
                    (df['angle'] < angle + 0.001) &
                    (df['waveNumber'] > waveNumber - 0.01) &
                    (df['waveNumber'] < waveNumber + 0.01)
                ]

                plt.hexbin(
                    sdf.loc[:, x_column],
                    sdf.loc[:, y_column],
                    gridsize=50,
                    reduce_C_function=np.mean,
                    cmap='inferno',
                    C=sdf.loc[:, 'error'], edgecolors='face', linewidths=0.5)

                # plt.contourf(
                #     sdf.loc[:, 'm'],
                #     sdf.loc[:, 'epsilon'],
                #     sdf.loc[:, 'error'])
                plt.colorbar()

                plt.xlabel(latex_label[x_column])
                plt.ylabel(latex_label[y_column])
                plt.title(
                    '$|| u_\epsilon - u_0 ||_{L^2 (Q / B(x, 0.25))} \
                    \ \mathrm{with} \ \gamma |_{B(x, e)} = %.2E, \
                    \ \\theta = %d^\circ, \ k=%.2E$' % (
                        conductivity, (180.0 / math.pi) * angle, waveNumber), pad=20)
                out=("%s_VS_%s_Gamma=%.2E_Angle=%.2E_Wavenumber=%.2E") % (
                    file_label[x_column], file_label[y_column], conductivity,
                    (180.0 / math.pi) * angle, waveNumber)
                plt.savefig(str(i) + '_' + out + '.svg')
                plt.savefig(str(i) + '_' + out + '.png')
                # plt.show()
                i += 1



