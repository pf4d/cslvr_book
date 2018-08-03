from pylab import *

mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['legend.fontsize'] = 'medium'

rhoi  = 917.0
rhow  = 1000.0
g     = 9.81
gamma = 9.8e-8 
Tw    = 273.15
ki    = 2.1
kw    = 0.561
ci    = 2009.0
cw    = 4217.6
a     = 60*60*24*365

H     = linspace(0, 5000, 1000)
w     = linspace(0,1,1000)
t     = Tw - gamma*rhoi*g*H
t     = linspace(250, 273, 1000)

T,W = meshgrid(t,w)

def ki_d(T):
  return 9.828 * exp(-0.0057*T)

def ci_d(T):
  return 146.3 + 7.253*T

ks   = (1-W)*ki   + W*kw
cs   = (1-W)*ci   + W*cw
rhos = (1-W)*rhoi + W*rhow

kd   = (1-W)*ki_d(T) + W*kw
cd   = (1-W)*ci_d(T) + W*cw

fig = figure(figsize=(8,10))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.plot(w, a*ks[:,0]/(cs[:,0] * rhos[:,0]), 'k', lw=2.0, label=r'$\kappa_d$')
ax1.set_xlabel(r'$W$')
ax1.set_ylabel(r'$\Xi_s(W)$')
ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

#cs1 = ax1.contourf(T,W, ci/ki*ks/cs)
#colorbar(cs1, ax=ax1)
#ax1.set_xlabel(r'$T$')
#ax1.set_ylabel(r'$W$')
ax1.grid()


cm = get_cmap('RdGy')
Z  = cs/ks*kd/cd
lv = linspace(Z.min(), Z.max(), 12)
lt = [(251.3,  0.0857),
      (252.67, 0.1914),
      (254.17, 0.3219),
      (256.22, 0.4244),
      (258.52, 0.5487),
      (261.37, 0.6605),
      (265.58, 0.6854),
      (269.43, 0.5456),
      (270.60, 0.2846),
      (271.97, 0.0950)]
#im  = ax2.imshow(Z, extent=(T.min(), T.max(), W.min(), W.max()))
cs1 = ax2.contourf(T,W,Z, cmap=cm, levels=lv, vmin=Z.min(), vmax=Z.max())
cs2 = ax2.contour(T,W, Z, linewidths=2.0, colors='k', levels=lv)
ax2.clabel(cs2, inline=1, colors='k', manual=lt, fmt='%1.2f')
#colorbar(cs2, ax=ax2)
ax2.set_xlabel(r'$T$')
ax2.set_ylabel(r'$W$')
ax2.grid()

tight_layout()
savefig('../../images/internal_energy/kappa_d.pdf')
show()


