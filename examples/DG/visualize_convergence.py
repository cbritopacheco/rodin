#!/usr/bin/env python3
"""Plot the measured SIPG refinement errors exported by Poisson.cpp."""
import argparse
from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, turbo_colors
apply_style()
COLORS = turbo_colors(2)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('root',type=Path)
parser.add_argument('output',type=Path)
args=parser.parse_args()
data=np.genfromtxt(args.root/'DGConvergence.csv',names=True,delimiter=',')
data=data[data['patch']==0]
h=1/data['intervals']
fig,ax=plt.subplots(figsize=(6,4.3),layout='constrained')
ax.loglog(h,data['l2'],'o-',color=COLORS[0],label=r'$\|u-u_h\|_{L^2(\Omega)}$')
ax.loglog(h,data['h1'],'s-',color=COLORS[1],label=r'$|u-u_h|_{H^1(\mathcal{T}_h)}$')
ax.loglog(h,data['l2'][-1]*(h/h[-1])**2,'--',color=COLORS[0],label=r'$O(h^2)$')
ax.loglog(h,data['h1'][-1]*(h/h[-1]),':',color=COLORS[1],label=r'$O(h)$')
ax.set_xlabel(r'$h=1/n$');ax.set_ylabel('Error')
ax.set_title('SIPG refinement study')
ax.spines[['top','right']].set_visible(False)
ax.legend(fontsize=10)
fig.savefig(args.output,dpi=180,transparent=True)
