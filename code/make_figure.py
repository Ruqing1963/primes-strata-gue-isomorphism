#!/usr/bin/env python3
"""Regenerate the Q-Q + CDF figure. Usage: cd code && python make_figure.py"""
import numpy as np, csv
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.signal import find_peaks, savgol_filter
from scipy.interpolate import interp1d
plt.rcParams.update({'figure.facecolor':'white','font.family':'serif','font.size':11})

def load_csv(path,col):
    v=[]
    with open(path) as f:
        for row in csv.reader(f):
            if not row or row[0].startswith('#') or not row[0][0].lstrip('-').isdigit(): continue
            v.append(float(row[col]))
    return np.array(v)

gamma=np.sort(load_csv('../data/riemann_zeros_real.csv',1))
t=gamma; N=(t/(2*np.pi))*np.log(t/(2*np.pi))-t/(2*np.pi)+7/8
s_riem=np.diff(N); s_riem=s_riem[s_riem>0]; s_riem/=s_riem.mean()
depth=load_csv('../data/gubbio_real.csv',0); chi=load_csv('../data/gubbio_real.csv',1)
idx=np.argsort(depth); depth,chi=depth[idx],chi[idx]
dz=0.02; grid=np.arange(depth.min(),depth.max(),dz)
chi_sm=savgol_filter(np.log(interp1d(depth,chi)(grid)),11,3)
peaks,_=find_peaks(chi_sm,distance=int(0.30/dz),prominence=0.1)
pos=np.sort(grid[peaks]); s_gub=np.diff(pos); s_gub=s_gub[s_gub>0]; s_gub/=s_gub.mean()

fig,axes=plt.subplots(1,2,figsize=(13,5.5))
ax=axes[0]
q=np.linspace(1,99,50); qx=np.percentile(s_riem,q); qy=np.percentile(s_gub,q)
ax.scatter(qx,qy,c='#c0392b',s=28,alpha=0.7,edgecolor='#7b241c',linewidth=0.4)
lim=[0,max(qx.max(),qy.max())*1.05]
ax.plot(lim,lim,'k--',lw=1.5,alpha=0.6,label='y = x (perfect isomorphism)')
slope,inter=np.polyfit(qx,qy,1)
ax.plot(np.array(lim),slope*np.array(lim)+inter,color='#2166ac',lw=1.8,label=f'best fit (slope={slope:.2f})')
ax.set_xlabel(f'Real Riemann zero spacing quantiles (n={len(s_riem)})')
ax.set_ylabel(f'Real Gubbio spacing quantiles (n={len(s_gub)})')
ax.set_title('Q-Q: Number Theory vs Stratigraphy (both real)')
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal'); ax.legend(fontsize=8.5); ax.grid(True,alpha=0.2)
ax=axes[1]
xr=np.sort(s_riem); yr=np.arange(1,len(xr)+1)/len(xr)
xg=np.sort(s_gub); yg=np.arange(1,len(xg)+1)/len(xg)
ax.plot(xr,yr,color='#c0392b',lw=2,label=f'Real Riemann zeros (n={len(xr)})')
ax.step(xg,yg,where='post',color='#1b7837',lw=2,label=f'Real Gubbio (n={len(xg)})')
def gue(s):return (32/np.pi**2)*s**2*np.exp(-4*s**2/np.pi)
sg=np.linspace(0,4,2000); mg=trapezoid(sg*gue(sg),sg)/trapezoid(gue(sg),sg)
cdf=cumulative_trapezoid(gue(sg),sg,initial=0); cdf/=cdf[-1]
ax.plot(sg/mg,cdf,'--',color='#8e44ad',lw=1.6,label='GUE theory')
ax.set_xlabel('Normalized spacing s'); ax.set_ylabel('Cumulative probability')
ax.set_title('Real zeros track GUE; Gubbio more rigid')
ax.set_xlim(0,3); ax.legend(fontsize=8.5); ax.grid(True,alpha=0.2)
fig.tight_layout()
fig.savefig('../figures/qq_isomorphism.pdf',dpi=300,bbox_inches='tight')
print("figure saved to figures/qq_isomorphism.pdf")
