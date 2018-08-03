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
u    = logspace(log(5e-3), log(100000.0),1000, base=e)

ks   = (1-w)*ki   + w*kw
cs   = (1-w)*ci   + w*cw
rhos = (1-w)*rhoi + w*rhow 

kap  = a * ks/(cs*rhos)   

U,K = meshgrid(u,kap)

A   = h*U / (2*K)
xi  = 1/tanh(A) - 1/A

lnK = log(K)
lnU = log(U)

fig = figure(figsize=(8,10))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

Amin = A.min()
Amax = A.max()

ximin = xi.min()
ximax = xi.max()

#A[A < Amin] = Amin + 1e-12
#A[A > Amax] = Amax - 1e-12

#xi[xi < ximin] = ximin + 1e-12
#xi[xi > ximax] = ximax - 1e-12

from matplotlib.ticker import LogFormatter
from matplotlib import colors

formatter   = LogFormatter(10, labelOnlyBase=False)
norm        = colors.LogNorm()
Alevels     = logspace(log(Amin), log(Amax), 12, base=e)
xilevels    = linspace(ximin, ximax, 12)

ltA = [(-5.579, 2.930),
       (-3.816, 2.858),
       (-2.247, 2.774),
       (-0.336, 2.722),
       (1.526, 2.638),
       (3.339, 2.560),
       (5.103, 2.534),
       (7.014, 2.469),
       (8.778, 2.417),
       (10.543, 2.294)]

ltXi = [(-3.7734, 3.21494),
        (-3.2683, 2.91657),
        (-3.2682, 2.44957),
        (-3.3693, 1.9177),
        (-1.6517, 3.1241),
        (-1.3486, 2.5469)]

cm  = get_cmap('RdGy')
cs1 = ax1.contourf(lnU, lnK, A, cmap=cm, levels=Alevels, norm=norm)
cs1 = ax1.contour(lnU, lnK, A, linewidths=2.0, colors='k',
                  levels=Alevels, norm=norm)
ax1.clabel(cs1, inline=1, colors='k', manual=ltA, fmt='%1.2e')
#colorbar(cs1, ax=ax1)
#ax1.set_xlabel(r'$\Vert \mathbf{u} \Vert$')
ax1.set_ylabel(r'$\Xi_s(W)$')
ax1.set_xticklabels(np.round(exp(ax1.get_xticks()), decimals=2))
ax1.set_yticklabels(np.round(exp(ax1.get_yticks()), decimals=2))
ax1.grid()

cs2 = ax2.contourf(lnU, lnK, xi, cmap=cm, levels=xilevels)#, norm=norm)
cs2 = ax2.contour(lnU, lnK, xi, linewidths=2.0, colors='k',
                  levels=xilevels)#, norm=norm)
ax2.clabel(cs2, inline=1, colors='k', manual=ltXi, fmt='%1.2f')
#colorbar(cs2, ax=ax2)
ax2.set_xlabel(r'$\Vert \mathbf{u} \Vert$')
ax2.set_ylabel(r'$\Xi_s(W)$')
ax2.set_xticklabels(np.round(exp(ax2.get_xticks()), decimals=2))
ax2.set_yticklabels(np.round(exp(ax2.get_yticks()), decimals=2))
ax2.grid()

tight_layout()
savefig('../../images/internal_energy/alpha.pdf')
show()



