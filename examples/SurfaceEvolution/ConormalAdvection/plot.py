import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import make_interp_spline
plt.style.use('grayscale')
# plt.style.use('seaborn-paper')

prefix = 'L2ErrorPhysical_'
suffix = '.csv'
expranges = {
    (1, 32): 0.1,
    (33, 64): 0.2,
    (65, 96): 0.3,
    (97, 128): 0.4,
    (129, 160): 0.5,
    (161, 192): 0.6,
    (193, 224): 0.7,
    (225, 256): 0.8,
    (257, 288): 0.9,
    (289, 320): 1.0,

}

meshsizes = np.linspace(0.02, 1.0, 32)
error = np.zeros(320)
for filename in glob.glob('%s*%s' % (prefix, suffix)):
    expid = int(filename[len(prefix):len(filename) - len(suffix)])
    rang = None
    for r in expranges:
        if r[0] <= expid and expid <= r[1]:
            rang = r
            break
    if (rang is None):
        pass
    with open(filename) as f:
        for line in f:
            pass
        last_line = line
        split = last_line.split(',')
        error[expid - 1] = float(split[1])

for r in expranges:
    spline = make_interp_spline(meshsizes, error[r[0] - 1:r[1]])
    x = np.linspace(meshsizes.min(), meshsizes.max(), 500)
    y = abs(spline(x))
    plt.plot(x, y, c=cm.RdYlBu_r(expranges[r]))
    # plt.plot(meshsizes, error[r[0] - 1:r[1]], c=cm.RdYlBu_r(expranges[r]))
    if (r == (289, 320)):
        break

cmap = plt.get_cmap('RdYlBu_r', 10)
sm = plt.cm.ScalarMappable(cmap=cmap)
test = np.mean(error.reshape((10, 32)), axis=0)
spline = make_interp_spline(meshsizes, test)
x = np.linspace(meshsizes.min(), meshsizes.max(), 500)
y = abs(spline(x))
cb = plt.colorbar(sm, ticks=np.linspace(0, 1, 10))
cb.ax.set_ylabel('$c$')
plt.plot(x, y, 'k--', linewidth=2,
    label='$\mathbb{E}[\mathcal{E}(t)]$')
plt.title('$ \delta t = c h $')
plt.ylabel('$\mathcal{E}(T)$', fontsize=14)
plt.xlabel('$h$', fontsize=14)
plt.legend(loc='upper left')
plt.grid(alpha=0.4, aa=True)
plt.savefig('Error.svg')
plt.show()


