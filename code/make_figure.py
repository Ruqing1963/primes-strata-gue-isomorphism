#!/usr/bin/env python3
"""Regenerate the two-section Q-Q + CDF figure. Usage: cd code && python make_figure.py"""
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
def strat(path,prom):
    h=load_csv(path,0); c=load_csv(path,1); idx=np.argsort(h); h,c=h[idx],c[idx]
    dz=0.02; grid=np.arange(h.min(),h.max()-1e-9,dz)
    val=np.log(interp1d(h,c)(grid)) if c.min()>0 else interp1d(h,c)(grid)
    pk,_=find_peaks(savgol_filter(val,11,3),distance=int(0.20/dz),prominence=prom)
    pos=np.sort(grid[pk]); s=np.diff(pos); s=s[s>0]; return s/s.mean()

gamma=np.sort(load_csv('../data/riemann_zeros_real.csv',1)); t=gamma
N=(t/(2*np.pi))*np.log(t/(2*np.pi))-t/(2*np.pi)+7/8
s_riem=np.diff(N); s_riem=s_riem[s_riem>0]; s_riem/=s_riem.mean()
sg=strat('../data/gubbio_real.csv',0.1); sb=strat('../data/bluelias_full.csv',0.002)

fig,axes=plt.subplots(1,2,figsize=(13.5,5.6))
ax=axes[0]; q=np.linspace(1,99,50); qx=np.percentile(s_riem,q)
for s,name,col in [(sg,'Gubbio','#c0392b'),(sb,'Blue Lias','#1b7837')]:
    ax.scatter(qx,np.percentile(s,q),s=24,alpha=0.7,color=col,label=f'{name} (n={len(s)})')
ax.plot([0,2.1],[0,2.1],'k--',lw=1.5,alpha=0.6,label='y = x (identity)')
ax.set_xlabel('Real Riemann zero spacing quantiles (n=4519)')
ax.set_ylabel('Real strata spacing quantiles')
ax.set_title('Q-Q: two real sections vs Riemann zeros')
ax.set_xlim(0,2.1); ax.set_ylim(0,2.1); ax.set_aspect('equal'); ax.legend(fontsize=8.5); ax.grid(True,alpha=0.2)
ax=axes[1]
for s,name,col in [(s_riem,f'Riemann zeros (n={len(s_riem)})','#c0392b'),
                   (sg,f'Gubbio (n={len(sg)})','#e08214'),(sb,f'Blue Lias (n={len(sb)})','#1b7837')]:
    xs=np.sort(s); ys=np.arange(1,len(xs)+1)/len(xs)
    (ax.plot if 'Riemann' in name else lambda *a,**k: ax.step(*a,where='post',**k))(xs,ys,color=col,lw=2,label=name)
def gue(s):return (32/np.pi**2)*s**2*np.exp(-4*s**2/np.pi)
sgx=np.linspace(0,4,2000); mg=trapezoid(sgx*gue(sgx),sgx)/trapezoid(gue(sgx),sgx)
cdf=cumulative_trapezoid(gue(sgx),sgx,initial=0); cdf/=cdf[-1]
ax.plot(sgx/mg,cdf,'--',color='#8e44ad',lw=1.5,label='GUE theory')
ax.set_xlabel('Normalized spacing s'); ax.set_ylabel('Cumulative probability')
ax.set_title('Both strata more rigid than the zeros'); ax.set_xlim(0,3); ax.legend(fontsize=8.5); ax.grid(True,alpha=0.2)
fig.tight_layout(); fig.savefig('../figures/qq_isomorphism.pdf',dpi=300,bbox_inches='tight')
print("figure saved")
