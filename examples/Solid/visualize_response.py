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
fig,axes=plt.subplots(1,2,figsize=(9,3.8),layout='constrained')
axes[0].plot(data['stretch'],data['energy'],color='#166699',linewidth=2)
axes[0].set_title('Stored energy');axes[0].set_ylabel(r'$W$')
axes[1].plot(data['stretch'],data['stress'],color='#166699',linewidth=2,label=r'$P_{11}$')
axes[1].plot(data['stretch'][::7],data['finite_difference'][::7],'o',color='#13b6d0',
             markersize=4,label=r'$dW/ds$ (finite difference)')
axes[1].set_title('Stress and energy derivative');axes[1].set_ylabel(r'$P_{11}$')
axes[1].legend(fontsize=9)
for ax in axes:
    ax.set_xlabel(r'$s$');ax.spines[['top','right']].set_visible(False)
fig.savefig(args.output,dpi=180,transparent=True)
