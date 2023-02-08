import glob
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('grayscale')

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
    plt.plot(meshsizes, error[r[0] - 1:r[1]], label=('c=%f' % expranges[r]))
    if (r == (161, 192)):
        break
plt.legend()
plt.grid(alpha=0.1, aa=True)
plt.show()


