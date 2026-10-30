#!/usr/bin/env python3
"""Plot energy and stress evaluated by ConstitutiveResponse.cpp."""
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style
apply_style()
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('root',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
data=np.genfromtxt(args.root/'ConstitutiveResponse.csv',delimiter=',',names=True)
for name in ['energy', 'stress']:
    fig, ax = plt.subplots(figsize=(4.8, 4.3), layout='constrained')
    if name == 'energy':
        ax.plot(data['stretch'], data['energy'], color='#166699', linewidth=2)
        ax.set_title('Stored energy')
        ax.set_ylabel(r'$W$')
    else:
        ax.plot(data['stretch'], data['stress'], color='#166699', linewidth=2, label=r'$P_{11}$')
        ax.plot(data['stretch'][::7], data['finite_difference'][::7], 'o', color='#13b6d0',
                markersize=4, label=r'$dW/ds$ (finite difference)')
        ax.set_title('Stress and energy derivative')
        ax.set_ylabel(r'$P_{11}$')
        ax.legend(fontsize=9)
    ax.set_xlabel(r'$s$')
    ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(args.output.with_name(args.output.stem + '-' + name + args.output.suffix),
                dpi=180, transparent=True)
    plt.close(fig)
