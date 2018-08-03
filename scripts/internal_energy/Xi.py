from pylab import *

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

rhoi = 917.0
rhow = 1000.0
ki   = 2.1
kw   = 0.561
ci   = 2009.0
cw   = 4217.6
a    = 60*60*24*365
     
h    = 1000.0
w    = logspace(log(1e-5), log(1.0), 1000, base=e)
#u    = linspace(1e-3, 50.0, 1000)
u    = logspace(log(2.3e-3), log(1.0), 1000, base=e)

ks   = (1-w)*ki   + w*kw
cs   = (1-w)*ci   + w*cw
rhos = (1-w)*rhoi + w*rhow 

kap  = a * ks/(cs*rhos)   

U,K = meshgrid(u,kap)

A   = h*U / (2*K)
xi  = 1/tanh(A) - 1/A

lnK = log(K)
lnU = log(U)

fig = figure(figsize=(8,5))
ax  = fig.add_subplot(111)

Amin = A.min()
Amax = A.max()

ximin = xi.min()
ximax = xi.max()

#xi[xi < ximin] = ximin + 1e-12
#xi[xi > ximax] = ximax - 1e-12

from matplotlib.ticker import LogFormatter
from matplotlib import colors

formatter   = LogFormatter(10, labelOnlyBase=False)
norm        = colors.LogNorm()
xilevels    = linspace(ximin, ximax, 12)

ltXi = [(-4.7328, 2.57940),
        (-3.7734, 3.21494),
        (-4.2663, 2.0604),
        (-3.2683, 2.91657),
        (-3.2682, 2.44957),
        (-2.3568, 3.11765),
        (-3.3693, 1.9177),
        (-1.6517, 3.1241),
        (-2.2342, 2.16418),
        (-1.3486, 2.5469)]

cm = get_cmap('RdGy')
cs = ax.contourf(lnU, lnK, xi, cmap=cm, levels=xilevels)
cs = ax.contour(lnU, lnK, xi, linewidths=2.0, colors='k', levels=xilevels)
ax.clabel(cs, inline=1, colors='k', manual=ltXi, fmt='%1.2f')
#colorbar(cs, ax=ax)
ax.set_xlabel(r'$\Vert \mathbf{u} \Vert$')
ax.set_ylabel(r'$\Xi_s(W)$')
ax.set_xticklabels(np.round(exp(ax.get_xticks()), decimals=2))
ax.set_yticklabels(np.round(exp(ax.get_yticks()), decimals=2))
ax.grid()

tight_layout()
savefig('../../images/internal_energy/alpha_small_u.pdf')
show()



